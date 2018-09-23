#include "BranchAllocationSystem.hpp"
#include "BranchArray.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "GeneIIDMultiDiscrete.hpp"
#include "GeneIIDMultiGamma.hpp"
#include "LogNormalMixArray.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDBeta.hpp"
#include "IIDNormal.hpp"
#include "IIDGamma.hpp"
#include "IIDDirichlet.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Sum.hpp"
#include "SimpleSubMatrixSelector.hpp"
#include "Tree.hpp"

/**
 * \brief A Branch-mixture selector
 */

template <class T>
class BranchMixtureSelector : public BranchSelector<T> {
  public:
    //! Constructor takes the array of components and the allocation vector
    BranchMixtureSelector(const Selector<T> &incomponents, const BranchAllocationSystem &inalloc)
        : components(incomponents), alloc(inalloc) {}
    ~BranchMixtureSelector() {}

    const Tree &GetTree() const override { return alloc.GetTree(); }
    T &operator[](int index) override { return components[alloc.GetBranchAlloc(index)]; }
    const T &GetVal(int index) const override {
        return components.GetVal(alloc.GetBranchAlloc(index));
    }

  private:
    const Selector<T> &components;
    const BranchAllocationSystem &alloc;
};

/**
 * \brief A site-homogeneous and branch-heterogeneous Muse and Gaut omega-codon
 * model
 *
 * The model has the following structure:
 * - branch lengths iid Exponential of rate lambda
 * - nucleotide relative exchangeabilities and stationaries are uniform
 * Dirichlet
 * - over branch j, omega_j = w * v_j, with v_j ~ Gamma(branchmean,
 * branchinvshape)
 * - w ~ Gamma(sitemean,siteinvshape)
 *
 * when model used in isolation: w = 1, and the v_j's are estimated, along with
 * their hyperparameters branchmean and branchinvshape. On the other hand, in a
 * multi gene context, the v_j's are shared across genes, and w_i for gene i is
 * allowed to vary (and then all 4 hyperparameters are estimated across genes).
 */

class SparseConditionOmegaModel : public ProbModel {
    // tree and data
    Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;
    int Ncond;
    int Nlevel;

    int blmode;
    int nucmode;

    // Branch lengths

    double lambda;
    BranchIIDGamma *blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise *branchlength;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    // suff stats branch lengths, as a function of their hyper parameter lambda
    // (bl are iid gamma, of scale parameter lambda)
    GammaSuffStat hyperlengthsuffstat;

    // Nucleotide rates

    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    vector<double> nucrelrate;
    vector<double> nucstat;
    GTRSubMatrix *nucmatrix;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // all this is fixed, being controlled at the inter-gene level
    double condvhypermean;
    double condvhypervar;
    IIDNormal *condv;

    double genew;

    Sum *meanlogomegaarray;

    vector<double> picenter;
    double piconcentration;
    IIDDirichlet* pi;
    IIDGamma* meanpos;
    IIDGamma* invshapepos;
    IIDGamma* meanneg;
    IIDGamma* invshapeneg;

    MultiDiscrete* alloc;
    MultiGamma* devpos;
    MultiGamma* devneg;

    // condition-specific omega's
    LogNormalMixArray *condomegaarray;

    // which branch is under which condition
    BranchAllocationSystem *branchalloc;

    // a codon matrix (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixArray *codonmatrixarray;
    MGOmegaCodonSubMatrix *rootcodonmatrix;
    SimpleSubMatrixSelector *submatrixtree;
    SimpleMGOmegaCodonSubMatrixSelector *codonmatrixtree;

    PhyloProcess *phyloprocess;

    // suffstats

    // suff stats for substitution paths
    PathSuffStatNodeArray *pathsuffstatarray;

    // or, alternatively, collected as Poisson suff stat across all nodes
    // (branches + root), as a function of omega
    OmegaPathSuffStatArray *omegapathsuffstatarray;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    SparseConditionOmegaModel(string datafile, string treefile, int inNcond, int inNlevel) {

        blmode = 0;
        nucmode = 0;

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();
        Ncond = inNcond;
        Nlevel = inNlevel;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        // specifies which condition for which branch
        branchalloc = new BranchAllocationSystem(*tree, Ncond);
        std::cerr << "-- conditions over branches ok\n";
    }

    //! model allocation
    void Allocate() {

        // Branch lengths

        lambda = 10.0;
        blhypermean = new BranchIIDGamma(*tree, 1.0, lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree, *blhypermean, 1.0 / blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);

        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        condvhypermean = 1.0;
        condvhypervar = 1.0;
        condv = new IIDNormal(Ncond, condvhypermean, condvhypervar);
        genew = 0.0;
        meanlogomegaarray = new Sum(*condv, genew);

        picenter.assign(3,0);
        picenter[0] = 0.05;
        picenter[2] = 0.05;
        picenter[1] = 0.90;
        piconcentration = 20;
        pi = new IIDDirichlet(Ncond,picenter,piconcentration);

        alloc = new MultiDiscrete(*pi);

        meanpos = new IIDGamma(Ncond,1.0,1.0);
        invshapepos = new IIDGamma(Ncond,1.0,1.0);
        meanneg = new IIDGamma(Ncond,1.0,1.0);
        invshapeneg = new IIDGamma(Ncond,1.0,1.0);

        devpos = new MultiGamma(*meanpos,*invshapepos);
        devneg = new MultiGamma(*meanneg,*invshapeneg);

        condomegaarray = new LogNormalMixArray(*meanlogomegaarray, *alloc, *devpos, *devneg);

        codonmatrixarray =
            new MGOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, condomegaarray);
        codonmatrixtree = new SimpleMGOmegaCodonSubMatrixSelector(*codonmatrixarray, *branchalloc);
        submatrixtree = new SimpleSubMatrixSelector(*codonmatrixarray, *branchalloc);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, 1.0);
        phyloprocess =
            new PhyloProcess(tree, codondata, branchlength, 0, submatrixtree, rootcodonmatrix);

        pathsuffstatarray = new PathSuffStatNodeArray(*tree);
        omegapathsuffstatarray = new OmegaPathSuffStatArray(Ncond);

        phyloprocess->Unfold();
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    LogNormalMixArray *GetOmegaArray() const { return condomegaarray; }

    void SetOmegaTree(const Selector<double> &from) { condomegaarray->Copy(from); }

    /*
    void SetOmegaHyperParameters(const Selector<double>& inpipos, const Selector<double>& inmeanpos, const Selector<double>& ininvshapepos, const Selector<double>& inpineg, const Selector<double>& inmeanneg, const Selector<double>& ininvshapeneg) {
        pipos->Copy(inpipos);
        meanpos->Copy(inmeanpos);
        invshapepos->Copy(ininvshapepos);
        pineg->Copy(inpineg);
        meanneg->Copy(inmeanneg);
        invshapeneg->Copy(ininvshapeneg);
    }
    */

    //-------------------
    // Setting and updating
    // ------------------

    //! \brief set estimation method for branch lengths and nuc rates
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode) {
        blmode = inblmode;
        nucmode = innucmode;
    }

    // Branch lengths

    //! whether branch lengths are fixed externally (e.g. when branch lengths are
    //! shared across genes in a multi-gene context)
    bool FixedBranchLengths() const { return blmode == 2; }

    //! set branch lengths to a new value (multi-gene analyses)
    void SetBranchLengths(const BranchSelector<double> &inbranchlength) {
        branchlength->Copy(inbranchlength);
    }

    //! get a copy of branch lengths into array given as argument
    void GetBranchLengths(BranchArray<double> &inbranchlength) const {
        inbranchlength.Copy(*branchlength);
    }

    //! set branch lengths hyperparameters to a new value (multi-gene analyses)
    void SetBranchLengthsHyperParameters(const BranchSelector<double> &inblmean,
                                                        double inblinvshape) {
        blhypermean->Copy(inblmean);
        blhyperinvshape = inblinvshape;
        branchlength->SetShape(1.0 / blhyperinvshape);
    }

    // Nucleotide rates

    //! whether nuc rates are fixed externally (e.g. when nuc rates are shared
    //! across genes in a multi-gene context)
    bool FixedNucRates() const { return nucmode == 2; }

    //! set nucleotide rates (relative exchangeabilities and eq. frequencies) to a
    //! new value (multi-gene analyses)
    void SetNucRates(const std::vector<double> &innucrelrate,
                                    const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        TouchMatrices();
    }

    //! get a copy of nucleotide rates into arrays given as arguments
    void GetNucRates(std::vector<double> &innucrelrate,
                                    std::vector<double> &innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double> &innucrelratehypercenter,
                                                   double innucrelratehyperinvconc,
                                                   const std::vector<double> &innucstathypercenter,
                                                   double innucstathyperinvconc) {
        nucrelratehypercenter = innucrelratehypercenter;
        nucrelratehyperinvconc = innucrelratehyperinvconc;
        nucstathypercenter = innucstathypercenter;
        nucstathyperinvconc = innucstathyperinvconc;
    }

    void SetGeneW(double inw) {
        genew = inw;
        meanlogomegaarray->SetGeneVal(genew);
    }

    void SetCondV(const Selector<double> &incondv) { condv->Copy(incondv); }

    void SetCondVHyperParams(double incondvhypermean, double incondvhypervar) {
        condvhypermean = incondvhypermean;
        condvhypervar = incondvhypervar;
        condv->SetMean(condvhypermean);
        condv->SetVar(condvhypervar);
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices() {
        TouchNucMatrix();
        codonmatrixarray->UpdateCodonMatrices();
        rootcodonmatrix->CorruptMatrix();
    }

    void Update() override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        condv->SetMean(condvhypermean);
        condv->SetVar(condvhypervar);
        meanlogomegaarray->SetGeneVal(genew);
        condomegaarray->Update();
        TouchMatrices();
        ResampleSub(1.0);
    }

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name) {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        condv->SetMean(condvhypermean);
        condv->SetVar(condvhypervar);
        meanlogomegaarray->SetGeneVal(genew);
        condomegaarray->Update();
        TouchMatrices();
        phyloprocess->PostPredSample(name);
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        if (!FixedBranchLengths()) {
            total += BranchLengthsLogPrior();
        }
        if (!FixedNucRates()) {
            total += NucRatesLogPrior();
        }
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const {
        double total = 0;
        if (blmode == 0) {
            total += LambdaHyperLogPrior();
        }
        total += branchlength->GetLogProb();
        return total;
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here,
    //! lambda ~ exponential of rate 10)
    double LambdaHyperLogPrior() const { return -lambda / 10; }

    // Nucleotide rates

    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq.
    //! freqs. (nucstat) -- uniform Dirichlet in both cases
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate, nucrelratehypercenter,
                                             1.0 / nucrelratehyperinvconc);
        total += Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    // Branch lengths

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current
    //! substitution mapping
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! \brief return log prob of current substitution mapping, as a function of
    //! branch lengths
    //!
    //! Calculated using the lengthpathsuffstat
    //! (which summarizes all information about how the prob of the substitution
    //! mapping depends on branch lengths). lengthpathsuffstat is assumed to be
    //! updated.
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    // Nucleotide rates

    //! \brief const acess to nuc-pathsuffstat
    //!
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate)
    //! and equilibrium frequencies (nucstat) conditional on the current
    //! substitution mapping.
    const NucPathSuffStat &GetNucPathSuffStat() const { return nucpathsuffstat; }

    //! \brief return log prob of current substitution mapping, as a function of
    //! nucleotide parameters (nucrelrate and nucstat)
    //!
    //! Calculated using nucpathsuffstat
    //! (which summarizes all information about how the probability of the
    //! substitution mapping depends on nucleotide mutation rates) and the
    //! nucmatrix. Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    const OmegaPathSuffStatArray *GetOmegaPathSuffStatArray() const {
        return omegapathsuffstatarray;
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrices. Both pathsuffstats and
    //! codonmatrices are assumed to be updated.
    double PathSuffStatLogProb() const {
        return pathsuffstatarray->GetLogProb(*submatrixtree, *rootcodonmatrix);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // Branch lengths

    //! \brief log prob factor to be recomputed when moving branch lengths
    //! hyperparameters (here, lambda)
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    // Nucleotide rates

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        TouchMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            if (!FixedBranchLengths()) {
                MoveBranchLengths();
            }

            CollectPathSuffStat();
            CollectOmegaSuffStat();

            if (!FixedNucRates()) {
                TouchMatrices();
                MoveNucRates();
            }
        }
    }

    // Branch lengths

    //! overall schedule branch length updatdes
    void MoveBranchLengths() {
        ResampleBranchLengths();
        if (blmode == 0) {
            MoveLambda();
        }
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &SparseConditionOmegaModel::LambdaHyperLogProb, &SparseConditionOmegaModel::NoUpdate,
                    this);
        ScalingMove(lambda, 0.3, 10, &SparseConditionOmegaModel::LambdaHyperLogProb, &SparseConditionOmegaModel::NoUpdate,
                    this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    // Nucleotide rates

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &SparseConditionOmegaModel::NucRatesLogProb,
                    &SparseConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &SparseConditionOmegaModel::NucRatesLogProb,
                    &SparseConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &SparseConditionOmegaModel::NucRatesLogProb,
                    &SparseConditionOmegaModel::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &SparseConditionOmegaModel::NucRatesLogProb,
                    &SparseConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &SparseConditionOmegaModel::NucRatesLogProb,
                    &SparseConditionOmegaModel::TouchNucMatrix, this);

        TouchMatrices();
    }

    // Omega

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    void CollectOmegaSuffStat() {
        omegapathsuffstatarray->Clear();
        omegapathsuffstatarray->AddSuffStat(*codonmatrixarray, *pathsuffstatarray, *branchalloc);
    }

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixtree, *rootcodonmatrix, *pathsuffstatarray);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "meanom\tvarom\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << condomegaarray->GetMean() << '\t' << condomegaarray->GetVar() << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const {}

    void ToStream(ostream &os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << *condomegaarray << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

    void FromStream(istream &is) {
        is >> lambda;
        is >> *branchlength;
        is >> *condomegaarray;
        is >> nucrelrate;
        is >> nucstat;
    }
};
