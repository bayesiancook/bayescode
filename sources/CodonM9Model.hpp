#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "IIDM9.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

/**
 * \brief A Muse and Gaut codon model with site-specific omega's
 * omega across sites from a mixture of 
 * - point mass at 0
 * - beta distribution between 0 and 1
 * - point mass at 1
 * - shifted gamma distribution for values > 1
 */

class CodonM9Model : public ProbModel {

    // tree and data
    const Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    const CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

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

    // Omega

    // weight of positive selection component
    double pi;
    double poswhypermean;
    double poswhyperinvconc;
	double posw;

    // all this with probability 1-posw
    // purifweight[0] : weight of point mass at omega = 0
    // purifweight[1] : weight for 0<omega<1 (beta distribution)
    // purifweight[2] : weight of point mass at omega = 1
    vector<double> purifweighthypercenter;
    double purifweighthyperinvconc;
    vector<double> purifweight;

    // mean and inverse concentration of the discretized Beta for 0 < omega < 1
    double purifmeanhypermean;
	double purifmeanhyperinvconc;
    double purifmean;

    double purifinvconchypermean;
    double purifinvconchyperinvshape;
    double purifinvconc;

    double posmeanhypermean;
    double posmeanhyperinvshape;
    double posmean;

    double posinvshapehypermean;
    double posinvshapehyperinvshape;
    double posinvshape;

    IIDM9 *omegaarray;

    // a codon matrix array (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixArray *codonmatrixarray;

    // PhyloProcess

    PhyloProcess *phyloprocess;

    // generic suff stats for substitution paths
    PathSuffStatArray *pathsuffstatarray;

    // or, alternatively, collected as a simple Poisson suff stat, as a function
    // of omega
    OmegaPathSuffStatArray *omegapathsuffstatarray;
    M9SuffStat omegahypersuffstat;

  public:

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    CodonM9Model(string datafile, string treefile, double inpi) {

        blmode = 0;
        nucmode = 0;

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        Tree* tmptree = new Tree(treefile);
        // check whether tree and data fits together
        tmptree->RegisterWith(taxonset);
        tmptree->SetIndices();
        tree = tmptree;

        Nbranch = tree->GetNbranch();

        pi = inpi;

        purifmeanhypermean = 0.5;
        purifmeanhyperinvconc = 0.5;
        purifinvconchypermean = 1.0;
        purifinvconchyperinvshape = 1.0;
        purifweighthypercenter.assign(3,1.0/3);
        purifweighthyperinvconc = 1.0;

        poswhypermean = 0.5;
        poswhyperinvconc = 0.1;
        posmeanhypermean = 1.0;
        posmeanhyperinvshape = 1.0;
        posinvshapehypermean = 1.0;
        posinvshapehyperinvshape = 1.0;
    }

    CodonM9Model(const CodonSequenceAlignment* incodondata, const Tree* intree, double inpi)   {

        blmode = 0;
        nucmode = 0;

        codondata = incodondata;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();

        pi = inpi;

        purifmeanhypermean = 0.5;
        purifmeanhyperinvconc = 0.5;
        purifinvconchypermean = 1.0;
        purifinvconchyperinvshape = 1.0;
        purifweighthypercenter.assign(3,1.0/3);
        purifweighthyperinvconc = 1.0;

        poswhypermean = 0.5;
        poswhyperinvconc = 0.1;
        posmeanhypermean = 1.0;
        posmeanhyperinvshape = 1.0;
        posinvshapehypermean = 1.0;
        posinvshapehyperinvshape = 1.0;
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

        // Omega

        purifweight.assign(3,0);
        copy(purifweighthypercenter.begin(), purifweighthypercenter.end(), purifweight.begin());
        posw = poswhypermean;

        purifmean = purifmeanhypermean;
        purifinvconc = purifinvconchypermean;

        posmean = posmeanhypermean;
        posinvshape = posinvshapehypermean;

        omegaarray = new IIDM9(Nsite, purifweight, posw, purifmean, purifinvconc, posmean, posinvshape);
        omegapathsuffstatarray = new OmegaPathSuffStatArray(Nsite);

        codonmatrixarray =
            new MGOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, omegaarray);

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, codonmatrixarray);
        pathsuffstatarray = new PathSuffStatArray(Nsite);

        phyloprocess->Unfold();
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    int GetNsite() const { return Nsite; }
    
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

    // Omega

    void GetMixtureParameters(double& inpurifmean, double& inpurifinvconc, double& inposw, double& inposmean, double& inposinvshape, vector<double>& inpurifweight)  const {
        inpurifmean = purifmean;
        inpurifinvconc = purifinvconc;
        inposw = posw;
        inposmean = posmean;
        inposinvshape = posinvshape;
        copy(purifweight.begin(), purifweight.end(), inpurifweight.begin());
    }

    void SetMixtureParameters(double inpurifmean, double inpurifinvconc, double inposw, double inposmean, double inposinvshape, const vector<double>& inpurifweight)    {
        purifmean = inpurifmean;
        purifinvconc = inpurifinvconc;
        posw = inposw;
        posmean = inposmean;
        posinvshape = inposinvshape;
        copy(inpurifweight.begin(), inpurifweight.end(), purifweight.begin());
        UpdateOmega();
    }

    void SetMixtureHyperParameters(
            double inpi, 
            double inpurifmeanhypermean, double inpurifmeanhyperinvconc,
            double inpurifinvconchypermean, double inpurifinvconchyperinvshape,
            const vector<double>& inpurifweighthypercenter, double inpurifweighthyperinvconc,
            double inposmeanhypermean, double inposmeanhyperinvshape,
            double inposinvshapehypermean, double inposinvshapehyperinvshape,
            double inposwhypermean, double inposwhyperinvconc)  {

        pi = inpi;

        purifmeanhypermean = inpurifmeanhypermean;
        purifmeanhyperinvconc = inpurifmeanhyperinvconc;
        purifinvconchypermean = inpurifinvconchypermean;
        purifinvconchyperinvshape = inpurifinvconchyperinvshape;
        purifweighthypercenter = inpurifweighthypercenter;
        purifweighthyperinvconc = inpurifweighthyperinvconc;

        posmeanhypermean = inposmeanhypermean;
        posmeanhyperinvshape = inposmeanhyperinvshape;
        posinvshapehypermean = inposinvshapehypermean;
        posinvshapehyperinvshape = inposinvshapehyperinvshape;

        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;

        if (!pi) {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }

        if (!purifmeanhyperinvconc) {
            purifmean = purifmeanhypermean;
        }
        
        if (! purifinvconchyperinvshape) {
            purifinvconc = purifinvconchypermean;
        }

        if (! purifweighthyperinvconc)  {
            copy(purifweighthypercenter.begin(), purifweighthypercenter.end(), purifweight.begin());
        }

        if (! posmeanhyperinvshape)  {
            posmean = posmeanhypermean;
        }

        if (! posinvshapehyperinvshape) {
            posinvshape = posinvshapehypermean;
        }

        if (! poswhyperinvconc) {
            posw = poswhypermean;
        }
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

    //! \brief tell the codon matrix that its parameters have changed and that it
    //! should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchCodonMatrices() {
        codonmatrixarray->UpdateCodonMatrices();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices() {
        TouchNucMatrix();
        TouchCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        UpdateOmega();
        TouchMatrices();
        ResampleSub(1.0);
    }

    void UpdateOmega()  {
        omegaarray->SetParameters(posw, purifmean, purifinvconc, posmean, posinvshape);
    }

    //-------------------
    // Posterior Predictive
    // ------------------

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name) override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        UpdateOmega();
        TouchMatrices();
        phyloprocess->PostPredSample(name);
    }

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
        total += OmegaHyperLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    // Branch lengths

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

    // Omega
    //! log prior over omega mixture
    double OmegaHyperLogPrior() const    {
        double total = 0;
        total += PurOmegaLogPrior();
        total += PosOmegaLogPrior();
        total += PurWeightLogPrior();
        total += PosWeightLogPrior();
        if (std::isinf(total))  {
            cerr << "omega hyper log prior: inf\n";
            cerr << PurOmegaLogPrior() << '\t' << PosOmegaLogPrior() << '\t' << PurWeightLogPrior() << '\t' << PosWeightLogPrior() << '\n';
            exit(1);
        }
        return total;
    }

    //! beta prior for purifmean
    //! gamma prior for purifinvconc
    double PurOmegaLogPrior() const {
        double total = 0;
        if (purifmeanhyperinvconc)  {
            double alpha = purifmeanhypermean / purifmeanhyperinvconc;
            double beta = (1 - purifmeanhypermean) / purifmeanhyperinvconc;
            total += Random::logBetaDensity(purifmean, alpha, beta);
        }
        if (purifinvconchyperinvshape)  {
            double alpha = 1.0 / purifinvconchyperinvshape;
            double beta = alpha / purifinvconchypermean;
            total += Random::logGammaDensity(purifinvconc, alpha, beta);
        }
        return total;
    }

    //! gamma prior for dposom
    double PosOmegaLogPrior() const {
        double total = 0;
        if (posmeanhyperinvshape)   {
            double shape = 1.0 / posmeanhyperinvshape;
            double scale = shape / posmeanhypermean;
            total += Random::logGammaDensity(posmean, shape, scale);
        }
        if (posinvshapehyperinvshape)   {
            double shape = 1.0 / posinvshapehyperinvshape;
            double scale = shape / posinvshapehypermean;
            total += Random::logGammaDensity(posinvshape, shape, scale);
        }
        return total;
    }

    //! dirichlet prior for purifweight
    double PurWeightLogPrior() const    {
        double total = 0;
        if (purifweighthyperinvconc)    {
            total += Random::logDirichletDensity(purifweight, purifweighthypercenter, 1.0/purifweighthyperinvconc);
        }
        return total;
    }

    //! mixture of point mass at 0 (with prob pi) and Beta distribution (with prob
    //! 1 - pi) for posw
    double PosWeightLogPrior() const    {
        if (posw) {
            if (!pi) {
                cerr << "in PosWeightLogProb: pi == 0 and posw > 0\n";
                exit(1);
            }

            if (! poswhyperinvconc)   {
                return 0;
            }
            double alpha = poswhypermean / poswhyperinvconc;
            double beta = (1 - poswhypermean) / poswhyperinvconc;
            return log(pi) + Random::logBetaDensity(posw, alpha, beta);
        } else {
            return log(1 - pi);
        }
    }

    double OmegaLogPrior() const { return omegaarray->GetLogProb(); }

    //! Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogPrior() const    {
        if (posw) {
            return log(pi);
        }
        return log(1 - pi);
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

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
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

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixarray, *pathsuffstatarray);
    }

    // Omega
    
    double OmegaHyperSuffStatLogProb() const {
        return omegahypersuffstat.GetLogProb(purifweight, posw, purifmean, purifinvconc, posmean, posinvshape);
    }


    // Paths 

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrix. Both pathsuffstat and
    //! codonmatrix are assumed to be updated.
    double PathSuffStatLogProb() const { return pathsuffstatarray->GetLogProb(*codonmatrixarray); }

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

    // for moving omegamean and omegainvshape
    double OmegaHyperLogProb() const { return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() override {
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

            MoveOmega();
            MoveOmegaHyperParameters();

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

    //! resample all branches not conditioned by sequence data from prior (as indicated by lengthpathsuffstats)
    void ResampleEmptyBranches()    {
        branchlength->ResampleEmptyBranches(*lengthpathsuffstatarray);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &CodonM9Model::LambdaHyperLogProb, &CodonM9Model::NoUpdate,
                    this);
        ScalingMove(lambda, 0.3, 10, &CodonM9Model::LambdaHyperLogProb, &CodonM9Model::NoUpdate,
                    this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    // Nucleotide rates

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &CodonM9Model::NucRatesLogProb,
                    &CodonM9Model::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &CodonM9Model::NucRatesLogProb,
                    &CodonM9Model::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &CodonM9Model::NucRatesLogProb,
                    &CodonM9Model::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &CodonM9Model::NucRatesLogProb,
                    &CodonM9Model::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &CodonM9Model::NucRatesLogProb,
                    &CodonM9Model::TouchNucMatrix, this);

        TouchMatrices();
    }

    // Omega

    void MoveOmega() {
        omegapathsuffstatarray->Clear();
        omegapathsuffstatarray->AddSuffStat(*codonmatrixarray, *pathsuffstatarray);
        omegaarray->MultipleTryMove(100, *omegapathsuffstatarray);
        TouchCodonMatrices();
    }

    void MoveOmegaHyperParameters()    {

        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegaarray);

        if (purifmeanhyperinvconc)  {
            SlidingMove(purifmean, 0.1, 10, 0, 1, &CodonM9Model::OmegaHyperLogProb, &CodonM9Model::NoUpdate, this);
        }
        if (purifinvconchyperinvshape)  {
            ScalingMove(purifinvconc, 0.1, 10, &CodonM9Model::OmegaHyperLogProb, &CodonM9Model::NoUpdate, this);
        }
        if (purifweighthyperinvconc)   {
            ProfileMove(purifweight, 0.1, 1, 10, &CodonM9Model::OmegaHyperLogProb, &CodonM9Model::NoUpdate, this);
        }
        if (pi != 0) {
            if (posmeanhyperinvshape)   {
                ScalingMove(posmean, 0.1, 10, &CodonM9Model::OmegaHyperLogProb, &CodonM9Model::NoUpdate, this);
            }
            if (posinvshapehyperinvshape)   {
                ScalingMove(posinvshape, 0.1, 10, &CodonM9Model::OmegaHyperLogProb, &CodonM9Model::NoUpdate, this);
            }
            if (poswhyperinvconc)   {
                SlidingMove(posw, 0.1, 10, 0, 1, &CodonM9Model::OmegaHyperLogProb, &CodonM9Model::NoUpdate, this);
            }
        }
        if ((pi != 0) && (pi != 1)) {
            SwitchPosWeight(10);
        }
        UpdateOmega();
    }

    //! reversible jump move on posw
    double SwitchPosWeight(int nrep)    {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double bkposw = posw;
            double deltalogprob = -PosSwitchLogPrior() - OmegaHyperLogProb();
            if (posw) {
                posw = 0;
            } else {
                double alpha = poswhypermean / poswhyperinvconc;
                double beta = (1 - poswhypermean) / poswhyperinvconc;
                posw = Random::BetaSample(alpha, beta);
            }
            deltalogprob += PosSwitchLogPrior() + OmegaHyperLogProb();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                posw = bkposw;
            }
            ntot++;
        }
        return nacc / ntot;
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    double GetEmpiricalPosFrac() const {
        double tot = 0;
        for (int i = 0; i < Nsite; i++) {
            if ((*omegaarray)[i] > 1.0) {
                tot++;
            }
        }
        return tot / Nsite;
    }

    void TraceOmega(ostream &os) const {
        for (int i = 0; i < GetNsite(); i++) {
            os << omegaarray->GetVal(i) << '\t';
        }
        os << '\n';
    }

    void GetSiteOmega(double *array) const {
        for (int i = 0; i < GetNsite(); i++) {
            array[i] = omegaarray->GetVal(i);
        }
    }

    double GetMeanOmega() const {
        return omegaarray->GetMean();
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "omegamean\t";
        os << "posfrac\t";
        os << "posw\tw0\tw1\tw2\t";
        os << "purifmean\tinvconc\tposmean\tinvshape\t";
        os << "statent\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << omegaarray->GetMean() << '\t';
        os << GetEmpiricalPosFrac() << '\t';
        os << posw << '\t' << purifweight[0] << '\t' << purifweight[1] << '\t' << purifweight[2] << '\t';
        os << purifmean << '\t' << purifinvconc << '\t' << posmean << '\t' << posinvshape << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << purifweight << '\t' << posw << '\t';
        os << purifmean << '\t' << purifinvconc << '\t' << posmean << '\t' << posinvshape << '\t';
        os << *omegaarray << '\t';
        os << nucstat << '\t';
        os << nucrelrate << '\t';
        os << lambda << '\t';
        os << *branchlength << '\n';
    }

    void FromStream(istream &is) override {
        is >> purifweight >> posw;
        is >> purifmean >> purifinvconc >> posmean >> posinvshape;
        is >> *omegaarray;
        is >> nucstat;
        is >> nucrelrate;
        is >> lambda;
        is >> *branchlength;
    }
};
