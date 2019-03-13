
/**
 * \brief The M2amodel of codeml (Muse and Gaut version)
 *
 * the omega_i's across sites (i=1..Nsite) are a mixture with 3 components
 * - omega0 < 1, with weight w0
 * - omega1 = 1, with weight w1
 * - omega2 > 1, with weight w2
 *
 * This model is used to test for the presence of positive selection (i.e to
 * test whether w2>0) in a gene and then to select those sites that have a dN/dS
 * > 1 (i.e. that are allocated to the third category of the mixture with high
 * post prob). Here, the model is parameterized as follows:
 * - omega0 = purom,
 * - omega1 = 1,
 * - omega2 = 1 + dposom,
 * where 0 < purom < 1 and dposom > 0;
 * purom has a beta prior (hyperparams: puromhypermean and puromhyperinvconc);
 * dposom has a gamma prior (hyperparams: dposomhypermean and
 * dposomhyperinvshape).
 *
 * The weights of the mixture are parameterized as follows:
 * - w0 = purw * (1 - posw)
 * - w1 = (1-purw) * (1-posw)
 * - w2 = posw
 * where 0<purw<1 and 0<=posw<1;
 * purw has a beta prior (hyperparams: purwhypermean and purwhyperinvconc);
 * the prior on posw is a mixture:
 * - with probability 1-pi, posw = 0
 * - with probability pi, 0 < posw < 1, in which case is it from a beta prior
 * (hyperparams: poswhypermean and poswhyperinvconc). Thus, setting pi = 0
 * imposes a model without positive selection.
 *
 * In total, the 9 hyperparameters of the mixture of omegas are:
 * puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape,
 * purwhypermean, purwhyperinvconc, pi, poswhypermean, poswhyperinvconc.
 * In a single-gene context, these hyperparameters are fixed;
 * in a multigene context, they can be either fixed or estimated across genes
 * (see MultiGeneCodonM2aModel).
 */

#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "M2aMix.hpp"
#include "Move.hpp"
#include "MultinomialAllocationVector.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "WhiteNoise.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"

// mode shared:      global
// mode shrunken:    gene specific, with hyperparameters estimated across genes
// mode independent: gene-specific, with fixed hyperparameters
enum param_mode_t { shared, shrunken, independent };

class CodonM2aModel : public ChainComponent {
  public:
    //-------------------
    // Constructors
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    CodonM2aModel(string datafile, string treefile, double pi)
        : datafile(datafile), treefile(treefile), pi(pi) {
        init();
        Update();
    }

    virtual ~CodonM2aModel() = default;

    void init();
    // CodonM2aModel(string datafile, string treefile, double inpi);

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        model_node(info, "nucstat", nucstat);
        model_node(info, "nucrelrate", nucrelrate);
        model_node(info, "branchlength", *branchlength);

        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        model_stat(info, "length", [this]() { return branchlength->GetTotalLength(); });
        // model_stat(info, "omega", omega);
        model_stat(info, "statent", [&]() { return Random::GetEntropy(nucstat); });
        model_stat(info, "rrent", [&]() { return Random::GetEntropy(nucrelrate); });
    }

    //! model allocation
    void Allocate();

    void ToStream(ostream &os) const {
        os << "CodonM2a" << '\t';
        os << datafile << '\t';
        os << treefile << '\t';
        os << pi << '\t';
        tracer->write_line(os);
    }

    CodonM2aModel(istream &is) {
        std::string model_name;
        is >> model_name;
        if (model_name != "CodonM2a") {
            std::cerr << "Expected CodonM2a for model name, got " << model_name << "\n";
            exit(1);
        }
        is >> datafile;
        is >> treefile;
        is >> pi;
        init();
        tracer->read_line(is);
        Update();
    }

    //-------------------
    // Accessors
    // ------------------

    //! number of aligned positions
    int GetNsite() const { return codondata->GetNsite(); }

    //! const access to codon state space
    const CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return value of omega_0 < 1
    double GetPurOm() const { return purom; }

    //! return value of omega_2 > 1
    double GetPosOm() const { return 1.0 + dposom; }

    //! return value of dposom  = omega_2 - 1 > 0
    double GetDPosOm() const { return dposom; }

    //! return proportion of sites under strictly purifying selection (under
    //! omega_0 < 1)
    double GetPurW() const { return purw; }

    //! return proportion of sites under positive selection
    double GetPosW() const { return posw; }

    //! whether branch lengths are fixed externally (e.g. when branch lengths are
    //! shared across genes in a multi-gene context)
    bool FixedBranchLengths() const { return blmode == 2; }

    //! whether nuc rates are fixed externally (e.g. when nuc rates are shared
    //! across genes in a multi-gene context)
    bool FixedNucRates() const { return nucmode == 2; }

    //-------------------
    // Setting and updating
    // ------------------

    // Setting model features and (hyper) parameters

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

    //! set branch lengths to a new value (multi-gene analyses)
    void SetBranchLengths(const BranchSelector<double> &inbranchlength);

    //! get a copy of branch lengths into array given as argument
    void GetBranchLengths(BranchArray<double> &inbranchlength) const;

    //! set branch lengths hyperparameters to a new value (multi-gene analyses)
    void SetBranchLengthsHyperParameters(
        const BranchSelector<double> &inblmean, double inblinvshape);

    //! resample all branches not conditioned by sequence data from prior (as indicated by
    //! lengthpathsuffstats)
    void ResampleEmptyBranches() { branchlength->ResampleEmptyBranches(*lengthpathsuffstatarray); }

    //! set nucleotide rates (relative exchangeabilities and eq. frequencies) to a
    //! new value (multi-gene analyses)
    void SetNucRates(const std::vector<double> &innucrelrate, const std::vector<double> &innucstat);

    //! get a copy of nucleotide rates into arrays given as arguments
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const;

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double> &innucrelratehypercenter,
        double innucrelratehyperinvconc, const std::vector<double> &innucstathypercenter,
        double innucstathyperinvconc);

    //! set omega mixture parameters to a new value
    void SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw);

    //! get omega mixture parameter values
    void GetMixtureParameters(
        double &inpurom, double &indposom, double &inpurw, double &inposw) const;

    //! set omega mixture hyperparameters to a new value
    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc,
        double indposomhypermean, double indposomhyperinvshape, double inpi, double inpurwhypermean,
        double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc);

    //-------------------
    // Matrices
    //-------------------

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update();

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name);

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void UpdateNucMatrix();

    //! \brief tell the codon matrices that their parameters have changed and that
    //! it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is
    //! notified, such that the matrices know that they will have to recalculate
    //! whichever component is requested later on upon demand.
    void UpdateCodonMatrices();

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to UpdateNucMatrix() and then UpdateCodonMatrices();
    void UpdateMatrices();

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //-------------------
    // Traces and Monitors
    // ------------------

    //! brief return current mean omega value
    double GetMeanOmega() const;

    /*
    void TraceHeader(ostream &os) const override;
    void Trace(ostream &os) const override;
    */

    //! \brief write current site post probs (of being under positive selection)
    //! on one line
    void TracePostProb(ostream &os) const;

    //! \brief get a copy of current site post probs (of being under positive
    //! selection) into array
    void GetSitesPostProb(double *array) const;

    // void Monitor(ostream& os) const override {}
    // void FromStream(istream &is) override;
    // void ToStream(ostream &os) const override;

    //-------------------
    // Likelihood
    //-------------------

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const {
        return GetLogPrior() + GetLogLikelihood();
        // return GetLogPrior() + GetIntegratedLogLikelihood();
    }

    //! return current value of likelihood, conditional on omega mixture
    //! allocations
    double GetLogLikelihood() const;

    //! return current value of likelihood, averaged over omega mixture
    //! allocations
    double GetIntegratedLogLikelihood() const;

    //-------------------
    // Priors
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const;

    //! log prior over branch lengths (iid exponential)
    double BranchLengthsLogPrior() const;

    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq.
    //! freqs. (nucstat) -- uniform Dirichlet in both cases
    double NucRatesLogPrior() const;

    //! log prior over omega mixture
    double OmegaLogPrior() const;

    //! beta prior for purom
    double PurOmegaLogPrior() const;

    //! gamma prior for dposom
    double PosOmegaLogPrior() const;

    //! beta prior for purw
    double PurWeightLogPrior() const;

    //! mixture of point mass at 0 (with prob pi) and Beta distribution (with prob
    //! 1 - pi) for posw
    double PosWeightLogPrior() const;

    //! Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogPrior() const;

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current
    //! substitution mapping
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const;

    //! \brief const acess to nuc-pathsuffstat
    //!
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate)
    //! and equilibrium frequencies (nucstat) conditional on the current
    //! substitution mapping.
    const NucPathSuffStat &GetNucPathSuffStat() const;

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrix. Both pathsuffstat and
    //! codonmatrix are assumed to be updated.
    double PathSuffStatLogProb() const;

    //! \brief return log prob of current substitution mapping, as a function of
    //! nucleotide parameters (nucrelrate and nucstat)
    //!
    //! Calculated using nucpathsuffstat
    //! (which summarizes all information about how the probability of the
    //! substitution mapping depends on nucleotide mutation rates) and the
    //! nucmatrix. Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesSuffStatLogProb() const;

    //! \brief return log prob of current substitution mapping, as a function of
    //! omega mixture configuration
    //!
    //! Calculated using siteomegapathsuffstatarray
    double OmegaPathSuffStatLogProb() const;

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    //! \brief log prob factor to be recomputed when moving parameters of omega
    //! mixture
    double OmegaLogProb() const { return OmegaLogPrior() + OmegaPathSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move();

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac);

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep);

    //
    // Branch Lengths
    //

    //! overall schedule branch length updatdes
    void MoveBranchLengths();

    //! Gibbs resample branch lengths (based on sufficient statistics)
    void ResampleBranchLengths();

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat();

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat();

    //
    // Omega mixture
    //

    //! collect sufficient statistics per component of omega mixtures
    void CollectComponentPathSuffStat();

    //! complete move schedule for omega mixture parameters
    void MoveOmega();

    //! collect sufficient statistics as a function of omega (per site)
    void CollectOmegaPathSuffStat();

    //! resample site allocations of omega mixture
    void ResampleAlloc();

    //! reversible jump move on posw
    double SwitchPosWeight(int nrep);

    //
    // nucleotide parameters
    //

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat();

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates();

    /*
    void SetLengthsFromTree() {
        cerr << "before set lengths: " << branchlength->GetTotalLength() << '\n';
        RecursiveSetLengthsFromTree(tree->GetRoot());
        cerr << "after set lengths: " << branchlength->GetTotalLength() << '\n';
    }

    void RecursiveSetLengthsFromTree(const Link *from) {
        if (!from->isRoot()) {
            double tmp = atof(from->GetBranch()->GetName().c_str());
            if (tmp <= 0) {
                cerr << "error: branch length is not positive: " << tmp << '\n';
                exit(1);
            }
            (*branchlength)[from->GetBranch()->GetIndex()] = tmp;
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveSetLengthsFromTree(link->Out());
        }
    }

    void FromStreamCodeML(istream &is);
    */

    //-------------------
    // Data structures
    // ------------------

  private:
    std::string datafile, treefile;
    std::unique_ptr<Tracer> tracer;
    std::unique_ptr<const Tree> tree;

    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    // Branch lengths
    double blhypermean;
    double blhyperinvshape;
    SimpleBranchArray<double> *blhypermeanarray;
    GammaWhiteNoise *branchlength;

    //
    // parameters of the distribution of omega across sites
    //

    // omega0 < 1: weight purw * (1 - posw)
    // omega1 = 1: weight (1-purw) * (1-posw)
    // omega2 > 1: weight posw

    // omega0 = purom
    double purom;

    // omega2 = 1 + dposom
    double dposom;

    double posw;
    double purw;

    M2aMix *componentomegaarray;
    MultinomialAllocationVector *sitealloc;
    mutable vector<vector<double>> sitepostprobarray;

    //
    // hyperparameters of the priors over the mixture parameters
    //

    // prior probability for the gene to be under positive selection (i.e. prior
    // prob that posw > 0)
    double pi;

    // Beta prior for purom (with hypermean and hyper inverse concentration)
    double puromhypermean;
    double puromhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape
    // parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    // Beta prior for purw
    double purwhypermean;
    double purwhyperinvconc;

    // Beta prior for posw (assuming posw>0)
    double poswhypermean;
    double poswhyperinvconc;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    // nucleotide rate parameters
    vector<double> nucrelrate;
    vector<double> nucstat;
    GTRSubMatrix *nucmatrix;

    // array of matrices across components of the mixture
    MGOmegaCodonSubMatrixArray *componentcodonmatrixarray;

    // arrays of matrices across sites (such as determined by the site allocations
    // to the mixture components) two versions differing only by their exact type

    // used for collecting omega suffstats: need to have access to the *codon*
    // matrix for each site
    MixtureSelector<MGOmegaCodonSubMatrix> *sitecodonmatrixarray;

    // used by PhyloProcess: has to be a Selector<SubMatrix>
    MixtureSelector<SubMatrix> *sitesubmatrixarray;

    PhyloProcess *phyloprocess;

    // suffstats

    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    OmegaPathSuffStatArray *siteomegapathsuffstatarray;
    PathSuffStatArray *sitepathsuffstatarray;
    PathSuffStatArray *componentpathsuffstatarray;

    NucPathSuffStat nucpathsuffstat;

    int blmode;
    int nucmode;
};
