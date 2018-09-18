#include "BranchArray.hpp"
#include "BranchProduct.hpp"
#include "BranchSpecificMeanGammaTree.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixBranchArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

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

class BranchOmegaModel : public ProbModel {
    // tree and data
    Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    // branch lengths iid expo (gamma of shape 1 and scale lambda)
    // where lambda is a hyperparameter
    double lambda;
    BranchIIDGamma *branchlength;

    // nucleotide exchange rates and equilibrium frequencies (stationary
    // probabilities)
    std::vector<double> nucrelrate;
    std::vector<double> nucstat;

    // all this is fixed, being controlled at the inter-gene level
    double branchvhypermean;
    double branchvhyperinvshape;
    BranchIIDGamma *branchv;
    double genew;
    BranchProduct *meanomegatree;
    double omegainvshape;

    // branch-specific omega's
    BranchSpecificMeanGammaTree *omegatree;

    // a nucleotide matrix (parameterized by nucrelrate and nucstat)
    GTRSubMatrix *nucmatrix;
    // a codon matrix (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixBranchArray *codonmatrixtree;
    MGOmegaCodonSubMatrix *rootcodonmatrix;

    PhyloProcess *phyloprocess;

    // suffstats

    // suff stats for substitution paths
    PathSuffStatNodeArray *pathsuffstatarray;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as Poisson suff stat across all nodes
    // (branches + root), as a function of omega
    OmegaPathSuffStatBranchArray *omegapathsuffstatarray;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    // suff stats branch lengths, as a function of their hyper parameter lambda
    // (bl are iid gamma, of scale parameter lambda)
    GammaSuffStat hyperlengthsuffstat;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    BranchOmegaModel(string datafile, string treefile) {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();
    }

    //! model allocation
    void Allocate() {
        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));

        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        branchvhypermean = 1.0;
        branchvhyperinvshape = 1.0;
        double alpha = 1.0 / branchvhyperinvshape;
        double beta = alpha / branchvhypermean;
        branchv = new BranchIIDGamma(*tree, alpha, beta);
        genew = 1.0;
        meanomegatree = new BranchProduct(*branchv, genew);
        omegainvshape = 1.0;
        omegatree = new BranchSpecificMeanGammaTree(*meanomegatree, omegainvshape);

        codonmatrixtree =
            new MGOmegaCodonSubMatrixBranchArray(GetCodonStateSpace(), nucmatrix, omegatree);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, 1.0);
        phyloprocess =
            new PhyloProcess(tree, codondata, branchlength, 0, codonmatrixtree, rootcodonmatrix);

        pathsuffstatarray = new PathSuffStatNodeArray(*tree);
        omegapathsuffstatarray = new OmegaPathSuffStatBranchArray(*tree);

        phyloprocess->Unfold();
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    BranchSpecificMeanGammaTree *GetOmegaTree() const { return omegatree; }

    void SetOmegaTree(const BranchSelector<double> &from) { omegatree->Copy(from); }

    //-------------------
    // Setting and updating
    // ------------------

    void SetGeneW(double inw) {
        genew = inw;
        meanomegatree->SetMulVal(genew);
    }

    void SetBranchV(const BranchSelector<double> &inbranchv) { branchv->Copy(inbranchv); }

    void SetBranchVHyperParams(double inbranchvhypermean, double inbranchvhyperinvshape) {
        branchvhypermean = inbranchvhypermean;
        branchvhyperinvshape = inbranchvhyperinvshape;
        double alpha = 1.0 / branchvhyperinvshape;
        double beta = alpha / branchvhypermean;
        branchv->SetShape(alpha);
        branchv->SetScale(beta);
    }

    void SetOmegaHyperInvShape(double inomegainvshape) {
        omegainvshape = inomegainvshape;
        omegatree->SetInvShape(omegainvshape);
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
        codonmatrixtree->UpdateCodonMatrices();
        rootcodonmatrix->CorruptMatrix();
    }

    void Update() override {
        branchlength->SetScale(lambda);
        double alpha = 1.0 / branchvhyperinvshape;
        double beta = alpha / branchvhypermean;
        branchv->SetShape(alpha);
        branchv->SetScale(beta);
        meanomegatree->SetMulVal(genew);
        omegatree->SetInvShape(omegainvshape);
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
    double OmegaLogPrior() const { return omegatree->GetLogProb(); }

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

    const OmegaPathSuffStatBranchArray *GetOmegaPathSuffStatArray() const {
        return omegapathsuffstatarray;
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrices. Both pathsuffstats and
    //! codonmatrices are assumed to be updated.
    double PathSuffStatLogProb() const {
        return pathsuffstatarray->GetLogProb(*codonmatrixtree, *rootcodonmatrix);
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
        ScalingMove(lambda, 1.0, 10, &BranchOmegaModel::BranchLengthsHyperLogProb,
            &BranchOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &BranchOmegaModel::BranchLengthsHyperLogProb,
            &BranchOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    void CollectOmegaSuffStat() {
        omegapathsuffstatarray->Clear();
        omegapathsuffstatarray->AddSuffStat(*codonmatrixtree, *rootcodonmatrix, *pathsuffstatarray);
        // cerr << "gene: " << omegapathsuffstatarray->GetMeanCount() << '\t' <<
        // omegapathsuffstatarray->GetMeanBeta() << '\t'; cerr <<
        // omegapathsuffstatarray->GetMeanCount() /
        // omegapathsuffstatarray->GetMeanBeta() << '\n';
    }

    //! Gibbs resample omega (based on sufficient statistics of current
    //! substitution mapping)
    void ResampleOmega() {
        omegatree->GibbsResample(*omegapathsuffstatarray);
        codonmatrixtree->UpdateCodonMatrices();
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

        ProfileMove(nucrelrate, 0.1, 1, 3, &BranchOmegaModel::NucRatesLogProb,
            &BranchOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &BranchOmegaModel::NucRatesLogProb,
            &BranchOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &BranchOmegaModel::NucRatesLogProb,
            &BranchOmegaModel::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &BranchOmegaModel::NucRatesLogProb,
            &BranchOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &BranchOmegaModel::NucRatesLogProb,
            &BranchOmegaModel::TouchNucMatrix, this);

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
        os << omegatree->GetMean() << '\t' << omegatree->GetVar() << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const {}

    void ToStream(ostream &os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << *omegatree << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

    void FromStream(istream &is) {
        is >> lambda;
        is >> *branchlength;
        is >> *omegatree;
        is >> nucrelrate;
        is >> nucstat;
    }
};
