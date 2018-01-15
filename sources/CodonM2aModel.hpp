
// this is the M2a model of codeml (Muse and Gaut version)
// this model assumes that the omega_i's across sites (i=1..Nsite)
// are a mixture with 3 components
//
// omega0 < 1, with weight w0
// omega1 = 1, with weight w1
// omega2 > 1, with weight w2
//
// this model is used to test for the presence of positive selection (i.e to test whether w2>0) in a gene
// and then to select those sites that have a dN/dS > 1 (i.e. that are allocated to the third category of the mixture with
// high post prob)
//
// here, the model is parameterized as follows:
// omega0 = purom,
// omega1 = 1,
// omega2 = 1 + dposom,
// where 0 < purom < 1 and dposom > 0
// purom has a beta prior (hyperparams: puromhypermean and puromhyperinvconc)
// dposom has a gamma prior (hyperparams: dposomhypermean and dposomhyperinvshape)
//
// the weights of the mixture are parameterized as follows:
// w0 = purw * (1 - posw)
// w1 = (1-purw) * (1-posw)
// w2 = posw
// where 0<purw<1 and 0<=posw<1
// purw has a beta prior (hyperparams: purwhypermean and purwhyperinvconc)
// the prior on posw is a mixture:
// - with probability 1-pi, posw = 0
// - with probability pi, 0 < posw < 1, in which case is it from a beta prior (hyperparams: poswhypermean and
// poswhyperinvconc)
// thus, setting pi = 0 imposes a model without positive selection
//
// in total, the 9 hyperparameters of the mixture of omegas are:
// puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape
// purwhypermean, purwhyperinvconc, pi, poswhypermean, poswhyperinvconc
// in a single-gene context, these hyperparameters are fixed
// in a multigene context, they can be either fixed or estimated across genes (see MultiGeneCodonM2aModel)

#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "M2aMix.hpp"
#include "MultinomialAllocationVector.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

class CodonM2aModel : public ProbModel {
  public:
    //-------------------
    // Constructors
    // ------------------

    CodonM2aModel(string datafile, string treefile, double inpi);
    void Allocate();
    void Unfold();

    //-------------------
    // Accessors
    // ------------------

    int GetNsite() const { return codondata->GetNsite(); }

    CodonStateSpace* GetCodonStateSpace() const { return (CodonStateSpace*)codondata->GetStateSpace(); }

    double GetPurOm() const { return purom; }

    double GetDPosOm() const { return dposom; }

    double GetPurW() const { return purw; }

    double GetPosW() const { return posw; }

    bool FixedBranchLengths() const { return blmode == 2; }

    bool FixedNucRates() const { return nucmode == 2; }

    //-------------------
    // Setting and updating
    // ------------------

    // Setting model features and (hyper) parameters

    // called upon constructing the model
    // mode == 2: global
    // mode == 1: gene specific, with hyperparameters estimated across genes
    // mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode) {
        blmode = inblmode;
        nucmode = innucmode;
    }

    void SetBranchLengths(const BranchSelector<double>& inbranchlength);
    void GetBranchLengths(BranchArray<double>& inbranchlength) const;

    void SetBranchLengthsHyperParameters(const BranchSelector<double>& inblmean, double inblinvshape);

    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat);
    void GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const;

    void SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc,
                                    const std::vector<double>& innucstathypercenter, double innucstathyperinvconc);

    void SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw);
    void GetMixtureParameters(double& inpurom, double& indposom, double& inpurw, double& inposw) const;

    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean,
                                   double indposomhyperinvshape, double inpi, double inpurwhypermean,
                                   double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc);

    // void Update() override {}
    void NoUpdate() {}

    //-------------------
    // Matrices
    //-------------------

    void Update() override;
    void UpdateNucMatrix();
    void UpdateCodonMatrices();
    void UpdateMatrices();

    //-------------------
    // Traces and Monitors
    // ------------------

    double GetMeanOmega() const;

    void TraceHeader(std::ostream& os) const override;
    void Trace(ostream& os) const override;
    void TracePostProb(ostream& os) const;
    void GetSitesPostProb(double* array) const;

    // void Monitor(ostream& os) const override {}
    void FromStream(istream& is) override;
    void ToStream(ostream& os) const override;

    //-------------------
    // Likelihood
    //-------------------

    double GetLogProb() const override { return GetLogPrior() + GetIntegratedLogLikelihood(); }

    double GetLogLikelihood() const;
    double GetIntegratedLogLikelihood() const;

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const;
    const NucPathSuffStat& GetNucPathSuffStat() const;

    double PathSuffStatLogProb() const;
    double LambdaHyperSuffStatLogProb() const;
    double NucRatesSuffStatLogProb() const;
    double OmegaPathSuffStatLogProb() const;

    //-------------------
    // Priors
    //-------------------

    double GetLogPrior() const;

    double LambdaHyperLogPrior() const;
    double BranchLengthsLogPrior() const;
    double NucRatesLogPrior() const;

    double OmegaLogPrior() const;

    // Beta prior for purifmean
    double PurOmegaLogProb() const;
    // Gamma prior for dposom
    double PosOmegaLogProb() const;
    // Beta prior for purw
    double PurWeightLogProb() const;
    // mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1 - pi) for posw
    double PosWeightLogProb() const;
    // Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogProb() const;

    //-------------------
    //  Log probs for MH moves
    //-------------------

    double LambdaHyperLogProb() const { return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb(); }

    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    double OmegaLogProb() const { return OmegaLogPrior() + OmegaPathSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    double Move();
    void ResampleSub(double frac);
    void MoveParameters(int nrep);

    //
    // Branch Lengths and hyperparam lambda
    //

    void MoveBranchLengths();
    void ResampleBranchLengths();
    void CollectLengthSuffStat();
    void MoveLambda();

    void CollectPathSuffStat();

    //
    // Omega mixture
    //

    void CollectComponentPathSuffStat();
    void MoveOmega();
    void CollectOmegaPathSuffStat();
    void ResampleAlloc();
    double DrawBetaPosWeight();
    double SwitchPosWeight(int nrep);

    //
    // nucleotide parameters
    //

    void CollectNucPathSuffStat();
    void MoveNucRates();

    //-------------------
    // Data structures
    // ------------------

  private:
    Tree* tree;
    FileSequenceAlignment* data;
    const TaxonSet* taxonset;
    CodonSequenceAlignment* codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    double lambda;
    BranchIIDGamma* blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise* branchlength;

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

    M2aMix* componentomegaarray;
    MultinomialAllocationVector* sitealloc;
    mutable vector<vector<double> > sitepostprobarray;

    //
    // hyperparameters of the priors over the mixture parameters
    //

    // prior probability for the gene to be under positive selection (i.e. prior prob that posw > 0)
    double pi;

    // Beta prior for purom (with hypermean and hyper inverse concentration)
    double puromhypermean;
    double puromhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape parameter)
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
    GTRSubMatrix* nucmatrix;

    // array of matrices across components of the mixture
    MGOmegaCodonSubMatrixArray* componentcodonmatrixarray;

    // arrays of matrices across sites (such as determined by the site allocations to the mixture components)
    // two versions differing only by their exact type

    // used for collecting omega suffstats: need to have access to the *codon* matrix for each site
    MixtureSelector<MGOmegaCodonSubMatrix>* sitecodonmatrixarray;

    // used by PhyloProcess: has to be a Selector<SubMatrix>
    MixtureSelector<SubMatrix>* sitesubmatrixarray;

    PhyloProcess* phyloprocess;

    // suffstats

    PoissonSuffStatBranchArray* lengthpathsuffstatarray;
    GammaSuffStat hyperlengthsuffstat;
    OmegaPathSuffStatArray* siteomegapathsuffstatarray;
    PathSuffStatArray* sitepathsuffstatarray;
    PathSuffStatArray* componentpathsuffstatarray;

    NucPathSuffStat nucpathsuffstat;

    int blmode;
    int nucmode;
};
