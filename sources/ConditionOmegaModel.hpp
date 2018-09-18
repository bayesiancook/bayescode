#include "BranchArray.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "ConditionSpecificMeanGammaArray.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Product.hpp"
#include "SimpleSubMatrixSelector.hpp"
#include "Tree.hpp"

/**
 * \brief A Branch-mixture selector
 */

template <class T>
class BranchMixtureSelector : public BranchSelector<T> {
  public:
    //! Constructor takes the array of components and the allocation vector
    BranchMixtureSelector(const Selector<T> &incomponents, const BranchSelector<int> &inalloc)
        : components(incomponents), alloc(inalloc) {}
    ~BranchMixtureSelector() {}

    const Tree &GetTree() const override { return alloc.GetTree(); }
    T &operator[](int index) override { return components[alloc.GetVal(index)]; }
    const T &GetVal(int index) const override { return components.GetVal(alloc.GetVal(index)); }

  private:
    const Selector<T> &components;
    const BranchSelector<int> &alloc;
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

class ConditionOmegaModel : public ProbModel {
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

    // branch lengths
    /*
        double lambda;
        BranchIIDGamma* blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise* branchlength;
        PoissonSuffStatBranchArray* lengthpathsuffstatarray;
        GammaSuffStat hyperlengthsuffstat;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    */

    double lambda;
    BranchIIDGamma *branchlength;
    GammaSuffStat hyperlengthsuffstat;

    // nucleotide exchange rates and equilibrium frequencies (stationary
    // probabilities)
    std::vector<double> nucrelrate;
    std::vector<double> nucstat;

    // all this is fixed, being controlled at the inter-gene level
    double condvhypermean;
    double condvhyperinvshape;
    IIDGamma *condv;

    double genew;

    Product *meanomegaarray;
    double omegainvshape;

    // condition-specific omega's
    ConditionSpecificMeanGammaArray *condomegaarray;

    // which branch is under which condition
    SimpleBranchArray<int> *branchalloc;

    // a nucleotide matrix (parameterized by nucrelrate and nucstat)
    GTRSubMatrix *nucmatrix;
    // a codon matrix (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixArray *codonmatrixarray;
    MGOmegaCodonSubMatrix *rootcodonmatrix;
    SimpleSubMatrixSelector *submatrixtree;
    SimpleMGOmegaCodonSubMatrixSelector *codonmatrixtree;

    PhyloProcess *phyloprocess;

    // suffstats

    // suff stats for substitution paths
    PathSuffStatNodeArray *pathsuffstatarray;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as Poisson suff stat across all nodes
    // (branches + root), as a function of omega
    OmegaPathSuffStatArray *omegapathsuffstatarray;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    ConditionOmegaModel(string datafile, string treefile, int inNcond, int inNlevel) {
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
        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);

        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));

        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        condvhypermean = 1.0;
        condvhyperinvshape = 1.0;
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condv = new IIDGamma(Ncond, alpha, beta);
        genew = 1.0;
        meanomegaarray = new Product(*condv, genew);
        omegainvshape = 1.0;
        condomegaarray = new ConditionSpecificMeanGammaArray(*meanomegaarray, omegainvshape);

        codonmatrixarray =
            new MGOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, condomegaarray);
        codonmatrixtree = new SimpleMGOmegaCodonSubMatrixSelector(*codonmatrixarray, *branchalloc);
        submatrixtree = new SimpleSubMatrixSelector(*codonmatrixarray, *branchalloc);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, 1.0);
        phyloprocess =
            new PhyloProcess(tree, codondata, branchlength, 0, submatrixtree, rootcodonmatrix);

        pathsuffstatarray = new PathSuffStatNodeArray(*tree);
        omegapathsuffstatarray = new OmegaPathSuffStatArray(Ncond);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        phyloprocess->Unfold();
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    ConditionSpecificMeanGammaArray *GetOmegaArray() const { return condomegaarray; }

    void SetOmegaTree(const Selector<double> &from) { condomegaarray->Copy(from); }

    //-------------------
    // Setting and updating
    // ------------------

    void SetGeneW(double inw) {
        genew = inw;
        meanomegaarray->SetMulVal(genew);
    }

    void SetCondV(const Selector<double> &incondv) { condv->Copy(incondv); }

    void SetCondVHyperParams(double incondvhypermean, double incondvhyperinvshape) {
        condvhypermean = incondvhypermean;
        condvhyperinvshape = incondvhyperinvshape;
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condv->SetShape(alpha);
        condv->SetScale(beta);
    }

    void SetOmegaHyperInvShape(double inomegainvshape) {
        omegainvshape = inomegainvshape;
        condomegaarray->SetInvShape(omegainvshape);
    }

    //! \brief set branch lengths to a new value
    //!
    //! Used in a multigene context.
    void SetBranchLengths(const BranchSelector<double> &inbranchlength) {
        branchlength->Copy(inbranchlength);
    }

    //! \brief set nucleotide relative exchangeabilities (nucrelrate) and
    //! equilibrium frequencies (nucstat) to a new value
    //!
    //! Notifies corruption to the nucleotide and the codon matrices
    void SetNucRates(
        const std::vector<double> &innucrelrate, const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        TouchMatrices();
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
        branchlength->SetScale(lambda);
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condv->SetShape(alpha);
        condv->SetScale(beta);
        meanomegaarray->SetMulVal(genew);
        condomegaarray->SetInvShape(omegainvshape);
        TouchMatrices();
        ResampleSub(1.0);
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
        /*
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        */
        total += OmegaLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    //! \brief log prior over hyperparameters of prior over branch lengths (here,
    //! lambda ~ exponential of rate 10)
    double BranchLengthsHyperLogPrior() const {
        // exponential of mean 10
        return -lambda / 10;
    }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const { return branchlength->GetLogProb(); }

    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq.
    //! freqs. (nucstat) -- uniform Dirichlet in both cases
    double NucRatesLogPrior() const { return 0; }

    //! log prior over omega
    double OmegaLogPrior() const { return condomegaarray->GetLogProb(); }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current
    //! substitution mapping
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! \brief const acess to nuc-pathsuffstat
    //!
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate)
    //! and equilibrium frequencies (nucstat) conditional on the current
    //! substitution mapping.
    const NucPathSuffStat &GetNucPathSuffStat() const { return nucpathsuffstat; }

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

    //! \brief return log prob of current substitution mapping, as a function of
    //! branch lengths
    //!
    //! Calculated using the lengthpathsuffstat
    //! (which summarizes all information about how the prob of the substitution
    //! mapping depends on branch lengths). lengthpathsuffstat is assumed to be
    //! updated.
    double BranchLengthsHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

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

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! \brief log prob factor to be recomputed when moving branch lengths
    //! hyperparameters (here, lambda)
    //!
    //! simply: BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb()
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    //!
    //! simply: NucRatesLogPrior() + NucRatesSuffStatLogProb();
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
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();
            CollectPathSuffStat();
            CollectOmegaSuffStat();
            ResampleOmega();
            MoveNucRates();
        }
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveBranchLengthsHyperParameter() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ;
        ScalingMove(lambda, 1.0, 10, &ConditionOmegaModel::BranchLengthsHyperLogProb,
            &ConditionOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &ConditionOmegaModel::BranchLengthsHyperLogProb,
            &ConditionOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    void CollectOmegaSuffStat() {
        omegapathsuffstatarray->Clear();
        omegapathsuffstatarray->AddSuffStat(*codonmatrixarray, *pathsuffstatarray, *branchalloc);
    }

    //! Gibbs resample omega (based on sufficient statistics of current
    //! substitution mapping)
    void ResampleOmega() {
        condomegaarray->GibbsResample(*omegapathsuffstatarray);
        codonmatrixarray->UpdateCodonMatrices();
    }

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixtree, *rootcodonmatrix, *pathsuffstatarray);
    }

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &ConditionOmegaModel::NucRatesLogProb,
            &ConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &ConditionOmegaModel::NucRatesLogProb,
            &ConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &ConditionOmegaModel::NucRatesLogProb,
            &ConditionOmegaModel::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &ConditionOmegaModel::NucRatesLogProb,
            &ConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &ConditionOmegaModel::NucRatesLogProb,
            &ConditionOmegaModel::TouchNucMatrix, this);

        TouchMatrices();
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
