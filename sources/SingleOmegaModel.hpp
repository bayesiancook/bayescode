#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "component_defs.hpp"

using namespace tc;
/**
 * \brief A standard site- and branch-homogeneous Muse and Gaut omega-codon model
 *
 * The model has the following structure:
 * - branch lengths iid Exponential of rate lambda
 * - nucleotide relative exchangeabilities and stationaries are uniform Dirichlet
 * - there is one single omega=dN/dS for all sites and across all branches
 * - prior over omega ~ Gamma(omegahypermean,omegahyperinvshape).
 *
 * The 2 hyperparameters omegahypermean and hyperinvshape are fixed when this model is used in isolation.
 * On the other hand, this model can be used in a multigene context (see MultiGeneSingeOmegaModel),
 * in which case omegahypermean and hyperinvshape are estimated across genes.
 */

class SingleOmegaModel : public ProbModel {
    Model model;                    // component version of graphical model
    unique_ptr<Assembly> assembly;  // pointer to assembly (FIXME allow empty assembly for ease of use)

    // tree and data
    Tree* tree;
    FileSequenceAlignment* data;
    const TaxonSet* taxonset;
    CodonSequenceAlignment* codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    // suffstats

    // suff stats for substitution paths
    // summed over all branches and over all sites
    PathSuffStat pathsuffstat;  // FIXME not component for now because undergoes copies

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function of omega
    OmegaPathSuffStat omegapathsuffstat;  // FIXME not component for now because undergoes copies

    // Poisson suffstats for substitution histories, as a function of branch lengths
    PoissonSuffStatBranchArray* lengthpathsuffstatarray;

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
    SingleOmegaModel(string datafile, string treefile) {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        cerr << "-- Number of sites: " << Nsite << endl;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();
    }

    void DeclareModel() {
        // branch lengths iid expo (gamma of shape 1 and scale lambda)
        // where lambda is a hyperparameter
        model.component<FWrapper<double>>("lambda", 10);
        model.component<BranchIIDGamma>("branchlength", tree, 1.0, 10);
        // model.component<PoissonSuffStatBranchArray>("lengthpathsuffstatarray", tree);

        // nucleotide exchange rates and equilibrium frequencies (stationary probabilities)
        model.component<Wrapper<vector<double>>>("nucrelrate", Nrr, 0);
        model.connect<DirichletSample>("nucrelrate", vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
        model.component<Wrapper<vector<double>>>("nucstat", Nnuc, 0);
        model.connect<DirichletSample>("nucstat", vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));

        // a nucleotide matrix (parameterized by nucrelrate and nucstat)
        model.component<GTRSubMatrix>("nucmatrix", Nnuc, true)
            .connect<Use<vector<double>>>("mRelativeRate", "nucrelrate")
            .connect<Use<vector<double>>>("CopyStationary", "nucstat");

        // omega has a Gamma prior
        // of mean omegahypermean and inverse shape parameter omegahyperinvshape
        model.component<FWrapper<double>>("omegahypermean", 1);
        model.component<FWrapper<double>>("omegahyperinvshape", 1);
        model.component<FWrapper<double>>("omega", 1);

        // a codon matrix (parameterized by nucmatrix and omega)
        // 1 = omega
        model.component<MGOmegaCodonSubMatrix>("codonmatrix", GetCodonStateSpace()->GetNstate(), 1)
            .connect<Set<CodonStateSpace*>>("statespace", GetCodonStateSpace())
            .connect<Use<SubMatrix>>("nucmatrix", "nucmatrix");

        model.component<PhyloProcess>("phyloprocess", tree, codondata)
            .set("siterate", static_cast<const Selector<double>*>(nullptr))
            .connect<Use<BranchSelector<double>>>("branchlength", "branchlength")
            .connect<Use<SubMatrix>>("submatrix", "nucmatrix");

        model.dot_to_file();

        assembly = unique_ptr<Assembly>(new Assembly(model));
        assembly->print_all();
        cout << (*assembly->at<Wrapper<vector<double>>>("nucrelrate"))[0] << endl;
        cout << (*assembly->at<Wrapper<vector<double>>>("nucstat"))[0] << endl;
        cout << *assembly->at<FWrapper<double>>("omega") << endl;
        cout << *assembly->at<FWrapper<double>>("lambda") << endl;
        cout << assembly->at<BranchIIDGamma>("branchlength").GetVal(2) << endl;
    }

    //! \brief unfold the model
    //!
    //! only at that step does the PhyloProcess create the whole structure of substitution mappings.
    void Unfold() {
        auto& phyloprocess = assembly->at<PhyloProcess>("phyloprocess");

        cerr << "-- unfold\n";
        phyloprocess.Unfold();
        cerr << phyloprocess.GetLogLikelihood() << '\n';
        cerr << "-- mapping substitutions\n";
        phyloprocess.ResampleSub();
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace* GetCodonStateSpace() const { return (CodonStateSpace*)codondata->GetStateSpace(); }

    //-------------------
    // Setting and updating
    // ------------------

    //! \brief set branch lengths to a new value
    //!
    //! Used in a multigene context.
    void SetBranchLengths(const BranchSelector<double>& inbranchlength) {
        assembly->at<BranchIIDGamma>("branchlength").Copy(inbranchlength);
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is notified,
    //! such that the matrix knows that it will have to recalculate whichever component is requested later on upon demand.
    void TouchNucMatrix() {
        auto& nucmatrix = assembly->at<GTRSubMatrix>("nucmatrix");
        nucmatrix.CopyStationary(&assembly->at<Wrapper<vector<double>>>("nucstat"));
        nucmatrix.CorruptMatrix();
    }

    //! \brief tell the codon matrix that its parameters have changed and that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is notified,
    //! such that the matrix knows that it will have to recalculate whichever component is requested later on upon demand.
    void TouchCodonMatrix() {
        auto& codonmatrix = assembly->at<MGOmegaCodonSubMatrix>("codonmatrix");
        auto& omega = *assembly->at<FWrapper<double>>("omega");
        codonmatrix.SetOmega(omega);
        codonmatrix.CorruptMatrix();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices() {
        TouchNucMatrix();
        TouchCodonMatrix();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove (defined in ProbModel),
    //! all of which require a void (*f)(void) function pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over all substitution histories)
    double GetLogLikelihood() const {
        auto& phyloprocess = assembly->at<PhyloProcess>("phyloprocess");
        return phyloprocess.GetLogLikelihood();
    }

    // //! return joint log prob (log prior + log likelihood)
    // double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    //! \brief log prior over hyperparameters of prior over branch lengths (here, lambda ~ exponential of rate 10)
    double BranchLengthsHyperLogPrior() const {
        // exponential of mean 10
        auto& lambda = *assembly->at<FWrapper<double>>("lambda");
        return -lambda / 10;
    }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const { return assembly->at<BranchIIDGamma>("branchlength").GetLogProb(); }

    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq. freqs. (nucstat) -- uniform Dirichlet in
    //! both cases
    double NucRatesLogPrior() const { return 0; }

    //! log prior over omega (gamma of mean omegahypermean and shape 1/omegahyperinvshape)
    double OmegaLogPrior() const {
        auto& omega = *assembly->at<FWrapper<double>>("omega");
        auto& omegahyperinvshape = *assembly->at<FWrapper<double>>("omegahyperinvshape");
        auto& omegahypermean = *assembly->at<FWrapper<double>>("omegahypermean");
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return Random::logGammaDensity(omega, alpha, beta);
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current substitution mapping
    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const { return lengthpathsuffstatarray; }

    //! \brief const acess to nuc-pathsuffstat
    //!
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate) and equilibrium frequencies (nucstat)
    //! conditional on the current substitution mapping.
    const NucPathSuffStat& GetNucPathSuffStat() const { return nucpathsuffstat; }

    //! \brief return log prob of the current substitution mapping, as a function of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the substitution mapping)
    //! and the codonmatrix.
    //! Both pathsuffstat and codonmatrix are assumed to be updated.
    double PathSuffStatLogProb() const {
        auto& codonmatrix = assembly->at<MGOmegaCodonSubMatrix>("codonmatrix");
        return pathsuffstat.GetLogProb(codonmatrix);
    }

    //! \brief return log prob of current substitution mapping, as a function of branch lengths
    //!
    //! Calculated using the lengthpathsuffstat
    //! (which summarizes all information about how the prob of the substitution mapping depends on branch lengths).
    //! lengthpathsuffstat is assumed to be updated.
    double BranchLengthsHyperSuffStatLogProb() const {
        auto& lambda = *assembly->at<FWrapper<double>>("lambda");
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    //! \brief return log prob of current substitution mapping, as a function of nucleotide parameters (nucrelrate and
    //! nucstat)
    //!
    //! Calculated using nucpathsuffstat
    //! (which summarizes all information about how the probability of the substitution mapping depends on nucleotide
    //! mutation rates)
    //! and the nucmatrix.
    //! Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesSuffStatLogProb() const {
        auto& nucmatrix = assembly->at<GTRSubMatrix>("nucmatrix");
        return nucpathsuffstat.GetLogProb(nucmatrix, *GetCodonStateSpace());
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! \brief log prob factor to be recomputed when moving branch lengths hyperparameters (here, lambda)
    //!
    //! simply: BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb()
    double BranchLengthsHyperLogProb() const { return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb(); }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation rate parameters (nucrelrate and nucstat)
    //!
    //! simply: NucRatesLogPrior() + NucRatesSuffStatLogProb();
    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() override {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter configuration
    void ResampleSub(double frac) {
        auto& phyloprocess = assembly->at<PhyloProcess>("phyloprocess");
        TouchMatrices();
        phyloprocess.Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();

            CollectPathSuffStat();

            MoveOmega();
            MoveNucRates();
        }
    }

    //! collect sufficient statistics for moving branch lengths (directly from the substitution mappings)
    void CollectLengthSuffStat() {
        auto& phyloprocess = assembly->at<PhyloProcess>("phyloprocess");
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(phyloprocess);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        assembly->at<BranchIIDGamma>("branchlength").GibbsResample(*lengthpathsuffstatarray);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda, based on suffstats for branch lengths)
    void MoveBranchLengthsHyperParameter() {
        auto& lambda = *assembly->at<FWrapper<double>>("lambda");
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(assembly->at<BranchIIDGamma>("branchlength"));
        ;
        ScalingMove(lambda, 1.0, 10, &SingleOmegaModel::BranchLengthsHyperLogProb, &SingleOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &SingleOmegaModel::BranchLengthsHyperLogProb, &SingleOmegaModel::NoUpdate, this);
        assembly->at<BranchIIDGamma>("branchlength").SetScale(lambda);
    }

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        auto& phyloprocess = assembly->at<PhyloProcess>("phyloprocess");
        pathsuffstat.Clear();
        pathsuffstat.AddSuffStat(phyloprocess);
    }

    //! Gibbs resample omega (based on sufficient statistics of current substitution mapping)
    void MoveOmega() {
        auto& omega = *assembly->at<FWrapper<double>>("omega");
        auto& omegahyperinvshape = *assembly->at<FWrapper<double>>("omegahyperinvshape");
        auto& omegahypermean = *assembly->at<FWrapper<double>>("omegahypermean");
        auto& codonmatrix = assembly->at<MGOmegaCodonSubMatrix>("codonmatrix");
        omegapathsuffstat.Clear();
        omegapathsuffstat.AddSuffStat(codonmatrix, pathsuffstat);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(), beta + omegapathsuffstat.GetBeta());
        TouchCodonMatrix();
    }

    //! collect sufficient statistics for moving nucleotide rates (based on generic sufficient statistics stored in
    //! pathsuffstat)
    void CollectNucPathSuffStat() {
        auto& codonmatrix = assembly->at<MGOmegaCodonSubMatrix>("codonmatrix");
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(codonmatrix, pathsuffstat);
    }

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        auto& nucrelrate = *assembly->at<Wrapper<vector<double>>>("nucrelrate");
        auto& nucstat = *assembly->at<Wrapper<vector<double>>>("nucstat");

        ProfileMove(nucrelrate, 0.1, 1, 3, &SingleOmegaModel::NucRatesLogProb, &SingleOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &SingleOmegaModel::NucRatesLogProb, &SingleOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &SingleOmegaModel::NucRatesLogProb, &SingleOmegaModel::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &SingleOmegaModel::NucRatesLogProb, &SingleOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &SingleOmegaModel::NucRatesLogProb, &SingleOmegaModel::TouchNucMatrix, this);

        TouchMatrices();
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream& os) const override {
        os << "#logprior\tlnL\tlength\tlambda\t";
        os << "omega\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream& os) const override {
        auto& lambda = *assembly->at<FWrapper<double>>("lambda");
        auto& omega = *assembly->at<FWrapper<double>>("omega");
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << assembly->at<BranchIIDGamma>("branchlength").GetTotalLength() << '\t';
        os << lambda << '\t';
        os << omega << '\t';
        os << Random::GetEntropy(*assembly->at<Wrapper<vector<double>>>("nucstat")) << '\t';
        os << Random::GetEntropy(*assembly->at<Wrapper<vector<double>>>("nucrelrate")) << '\n';
    }

    void Monitor(ostream& os) const override {}

    void ToStream(ostream& os) const override {
        auto& lambda = *assembly->at<FWrapper<double>>("lambda");
        auto& omega = *assembly->at<FWrapper<double>>("omega");
        os << lambda << '\n';
        os << assembly->at<BranchIIDGamma>("branchlength") << '\n';
        os << omega << '\n';
        os << *assembly->at<Wrapper<vector<double>>>("nucrelrate") << '\n';
        os << *assembly->at<Wrapper<vector<double>>>("nucstat") << '\n';
    }

    void FromStream(istream& is) override {
        auto& lambda = *assembly->at<FWrapper<double>>("lambda");
        auto& omega = *assembly->at<FWrapper<double>>("omega");
        is >> lambda;
        is >> assembly->at<BranchIIDGamma>("branchlength");
        is >> omega;
        is >> *assembly->at<Wrapper<vector<double>>>("nucrelrate");
        is >> *assembly->at<Wrapper<vector<double>>>("nucstat");
    }
};
