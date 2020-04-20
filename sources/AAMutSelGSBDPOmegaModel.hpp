#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "Chrono.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "MultinomialAllocationVector.hpp"
#include "Permutation.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "StickBreakingProcess.hpp"
#include "Tree.hpp"
#include "AADiffSelCodonMatrixBidimArray.hpp"
#include "EmpiricalSelACProfileBidimArray.hpp"


/**
 * \brief The mutation-selection model with constant fitness landscape over the
 * tree -- double Dirichlet process version.
 *
 * The model is parameterized by
 * - an fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for
 * j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process
 * (assumed homogeneous across sites and lineages)
 * - an array of site-specific amino-acid fitness profiles F_ia, for site i and
 * amino-acid a
 * - an omega multiplier, capturing deviations of the non-syn rate from the
 * model (see Rodrigue and Lartillot, 2107); this parameter is fixed to 1 by
 * default.
 *
 * Site-specific amino-acid fitness profiles are drawn from a Dirichlet process,
 * implemented using a truncated stick-breaking process, of concentration
 * parameter kappa, and truncated at Ncat.
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 * - omega: fixed to 1 or exponential of rate 1
 * - kappa: exponential of rate 0.1
 * - center of base distribution: uniform Dirichlet
 * - concentration of base distribution: exponential of mean 20.
 *
 * In a multi-gene context, shrinkage across genes can be applied to branch
 * lengths, omega, nucleotide rate parameters (rho and pi), and to the
 * parameters of the base distribution (center and concentration) -- see
 * MultiGeneAAMutSelGSBDPModel.
 *
 */

class AAMutSelGSBDPOmegaModel : public ProbModel {
    const Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    const CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    double lambda;
    BranchIIDGamma *blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise *branchlength;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStat hyperlengthsuffstat;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucstat;
    std::vector<double> nucrelrate;
    GTRSubMatrix *nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
    double omega;
    OmegaPathSuffStat omegapathsuffstat;

    double maxdposom;

    // base distribution G0 is itself a stick-breaking mixture of Dirichlet
    // distributions

    // aa fitness arrays across sites are a SBDP process of base G0 defined above
    int Ncat;
    double kappa;
    StickBreakingProcess *weight;
    OccupancySuffStat *occupancy;

    vector<double> basecenter;
    double baseconcentration;

    IIDDirichlet *componentaafitnessarray;
    DirichletSuffStat basesuffstat;

    int Gcat;
    double Ginvshape;
    vector<double> G;
    double psi;

    MultinomialAllocationVector *sitealloc;

    vector<double> Gweight;
    MultinomialAllocationVector* Galloc;

    UnconstrainedSelACProfileBidimArray* selacprofiles;
    AADiffSelCodonMatrixBidimArray* codonmatrices;

    DoubleMixtureSelector<SubMatrix> *sitesubmatrixarray;

    PhyloProcess *phyloprocess;

    PathSuffStatArray *sitepathsuffstatarray;
    PathSuffStatBidimArray *componentpathsuffstatbidimarray;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    // currently: shared across genes
    int blmode;
    // currently, free without shrinkage: shared across genes
    int nucmode;
    // currently: fixed or free with shrinkage
    int omegamode;

    double dposompi;
    double dposomhypermean;
    double dposomhyperinvshape;

    // 0: simple gamma prior
    // 1: mix of (1-pi) at 1 and pi at 1+d, with d ~
    // Gamma(dposomhypermean,dposomhyperinvshape) 2: mix of 2 (modal) gamma
    // distributions: one at 1 and another one with mean > 1
    int omegaprior;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions,
    //! such as specified by branch names)
    //! - inomegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    //! - Ncat: truncation of the first-level stick-breaking process (by default:
    //! 100)
    //! - baseNcat: truncation of the second-level stick-breaking process (by
    //! default: 1)

    AAMutSelGSBDPOmegaModel(string datafile, string treefile, int inomegamode, int inomegaprior, int inNcat, int inGcat) :
        basesuffstat(Naa) {
        blmode = 0;
        nucmode = 0;
        omegamode = inomegamode;
        omegaprior = inomegaprior;
        maxdposom = 0;

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        Ncat = inNcat;
        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100) {
                Ncat = 100;
            }
        }

        Gcat = inGcat;

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        cerr << "ncat : " << Ncat << '\n';
        cerr << "gcat : " << Gcat << '\n';

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        Tree* tmptree = new Tree(treefile);
        // check whether tree and data fits together
        tmptree->RegisterWith(taxonset);
        tmptree->SetIndices();
        tree = tmptree;

        Nbranch = tree->GetNbranch();
    }

    AAMutSelGSBDPOmegaModel(const CodonSequenceAlignment* incodondata, const Tree* intree, int inomegamode, int inomegaprior, int inNcat, int inGcat) :
        basesuffstat(Naa) {

        blmode = 0;
        nucmode = 0;
        omegamode = inomegamode;
        omegaprior = inomegaprior;
        maxdposom = 0;

        codondata = incodondata;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        Ncat = inNcat;
        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100) {
                Ncat = 100;
            }
        }

        Gcat = inGcat;

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        cerr << "ncat : " << Ncat << '\n';
        cerr << "gcat : " << Gcat << '\n';

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();
    }

    void SetMaxDPosOm(double inmax)  {
        maxdposom = inmax;
    }

    void SetBLMode(int mode) { blmode = mode; }

    void SetNucMode(int mode) { nucmode = mode; }

    void SetOmegaMode(int mode) { omegamode = mode; }

    void SetOmegaPrior(int inprior) { omegaprior = inprior; }

    //! allocate the model (data structures)
    void Allocate() {
        // branch lengths
        lambda = 10;
        blhypermean = new BranchIIDGamma(*tree, 1.0, lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree, *blhypermean, 1.0 / blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));
        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        psi = 1.0;
        Ginvshape = 1.0;
        G.assign(Gcat, 1.0);
        UpdateG();

        basecenter.assign(Naa, 1.0 / Naa);
        baseconcentration = double(Naa);

        // Ncat fitness profiles iid from the base distribution
        componentaafitnessarray = new IIDDirichlet(Ncat, basecenter, baseconcentration);

        // mixture weights (truncated stick breaking process)
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat, kappa);

        // site allocations to the mixture (multinomial allocation)
        sitealloc = new MultinomialAllocationVector(Nsite, weight->GetArray());

        // occupancy suff stats of site allocations (for resampling weights)
        occupancy = new OccupancySuffStat(Ncat);

        Gweight.assign(Gcat, 1.0/Gcat);
        Galloc = new MultinomialAllocationVector(Nsite,Gweight);

        selacprofiles = new UnconstrainedSelACProfileBidimArray(*componentaafitnessarray, G, psi);

        // global omega (fixed to 1 by default)
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        omega = 1.0;
        dposompi = 0.1;
        dposomhypermean = 1;
        dposomhyperinvshape = 0.5;

        codonmatrices = new AADiffSelCodonMatrixBidimArray(*selacprofiles, *GetCodonStateSpace(), *nucmatrix, omega);
        sitesubmatrixarray = new DoubleMixtureSelector<SubMatrix>(codonmatrices, sitealloc, Galloc);

        // create polyprocess

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, sitesubmatrixarray);
        phyloprocess->Unfold();

        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatbidimarray = new PathSuffStatBidimArray(Ncat,Gcat);
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //! return current omega value
    double GetOmega() const { return omega; }

    //! \brief const access to array of length-pathsuffstats across branches
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //-------------------
    // Setting and updating
    // ------------------

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

    //! set omega to new value (multi-gene analyses)
    void SetOmega(double inomega) {
        omega = inomega;
        UpdateCodonMatrices();
    }

    //! set omega hyperparams to new value (multi-gene analyses)
    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape) {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    void SetDPosOmHyperParameters(double inpi, double inmean, double ininvshape) {
        dposompi = inpi;
        dposomhypermean = inmean;
        dposomhyperinvshape = ininvshape;
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

    //! set nucleotide rates to a new value (multi-gene analyses)
    void SetNucRates(const std::vector<double> &innucrelrate,
                     const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

    //! copy nucleotide rates into vectors given as arguments (multi-gene
    //! analyses)
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    void SetBase(const vector<double>& center, double concentration)    {
        copy(center.begin(), center.end(), basecenter.begin());
        baseconcentration = concentration;
        componentaafitnessarray->SetConcentration(baseconcentration);
    }

    void SetG(double inGinvshape)   {
        Ginvshape = inGinvshape;
        UpdateSelAC();
    }

    void UpdateG()  {
        Random::DiscGamma(G,1.0/Ginvshape);
    }

    void UpdateSelAC()  {
        UpdateG();
        selacprofiles->SetPsi(psi);
        selacprofiles->Update();
        UpdateCodonMatrices();
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatrices() {
        codonmatrices->SetOmega(omega);
        codonmatrices->Corrupt();
    }

    void UpdateCodonMatrices(const Permutation& permut)   {
        for (int i=0; i<Ncat; i++)  {
            if (permut.GetVal(i) != i) {
                UpdateCodonMatrix(i);
            }
        }
    }

    //! \brief tell codon matrix k that its parameters have changed and that it
    //! should be updated
    void UpdateCodonMatrix(int k) {
        selacprofiles->UpdateRow(k);
        codonmatrices->CorruptRow(k); 
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void Update() override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        weight->SetKappa(kappa);
        componentaafitnessarray->SetConcentration(baseconcentration);
        UpdateOccupancies();
        UpdateSelAC();
        UpdateMatrices();
        ResampleSub(1.0);
    }

    void PostPred(string name) override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        weight->SetKappa(kappa);
        componentaafitnessarray->SetConcentration(baseconcentration);
        UpdateOccupancies();
        UpdateSelAC();
        UpdateMatrices();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior (up to some constant)
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) {
            total += BranchLengthsLogPrior();
        }
        if (nucmode < 2) {
            total += NucRatesLogPrior();
        }
        total += GLogPrior();
        total += BaseLogPrior();
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AALogPrior();
        if (omegamode < 2) {
            total += OmegaLogPrior();
        }
        return total;
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    double GLogPrior() const    {
        return -Ginvshape;
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here,
    //! lambda ~ exponential of rate 10)
    double LambdaHyperLogPrior() const { return -lambda / 10; }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const {
        double ret = branchlength->GetLogProb();
        if (blmode == 0) {
            ret += LambdaHyperLogPrior();
        }
        return ret;
    }

    //! log prior over omega (gamma of mean omegahypermean and inverse shape
    //! omegahyperinvshape)
    double OmegaLogPrior() const {
        double ret = 0;

        // gamma distribution over 0, +infty
        if (omegaprior == 0) {
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            ret = alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(omega) -
                  beta * omega;
        }
        
        // mixture of point mass at 1 and shifted gamma
        else if (omegaprior == 1) {
            if ((dposompi <= 0) || (dposompi >= 1)) {
                cerr << "error in omegalogprior: pi is not 0<pi<1\n";
                exit(1);
            }
            if (omega < 1.0) {
                cerr << "error in omegalogprior: omega < 1\n";
                exit(1);
            }
            double dposom = omega - 1.0;
            if (dposom == 0) {
                ret = log(1 - dposompi);
            } else {
                ret = log(dposompi);
                double alpha = 1.0 / dposomhyperinvshape;
                double beta = alpha / dposomhypermean;
                ret += alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(dposom) -
                       beta * dposom;
            }
        }
        
        // mixture of mass at 0 + gamma over 0,+infty in log
        else if (omegaprior == 2) {
            if ((dposompi <= 0) || (dposompi >= 1)) {
                cerr << "error in omegalogprior: pi is not 0<pi<1\n";
                exit(1);
            }
            if (omega < 1.0) {
                cerr << "error in omegalogprior: omega < 1\n";
                exit(1);
            }
            double dposom = log(omega);
            if (dposom == 0) {
                ret = log(1 - dposompi);
            } else {
                ret = log(dposompi);
                double val = log(1 + dposompi);
                double alpha = 1.0 / dposomhyperinvshape;
                double beta = alpha / dposomhypermean;
                ret += alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(val) -
                       beta * val;
            }
        }
        
        // mixture of mass at 1 + half cauchy (shifted)
        else if (omegaprior == 3) {
            if ((dposompi <= 0) || (dposompi >= 1)) {
                cerr << "error in omegalogprior: pi is not 0<pi<1\n";
                exit(1);
            }
            if (omega < 1.0) {
                cerr << "error in omegalogprior: omega < 1\n";
                exit(1);
            }
            double dposom = omega - 1.0;
            if (dposom == 0) {
                ret = log(1 - dposompi);
            } else {
                ret = log(dposompi);
                double gamma = 1.0 / dposomhyperinvshape;
                // truncated Cauchy, both sides (0,maxdposom)
                ret += log(Pi * gamma * (1 + (dposom/gamma)*(dposom/gamma)));
                if (maxdposom)  {
                    ret -= log(2.0 / Pi * atan(maxdposom/gamma));
                }
            }
        } else {
            cerr << "error in OmegaLogPrior: unrecognized prior mode\n";
            exit(1);
        }
        return ret;
    }

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate, nucrelratehypercenter,
                                             1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //! log prior over concentration parameters kappa of stick-breaking mixture of
    //! amino-acid profiles
    double StickBreakingHyperLogPrior() const { return -kappa / 10; }

    //! log prior over weights of stick breaking process of amino-acid profiles
    double StickBreakingLogPrior() const { return weight->GetLogProb(kappa); }

    //! log prior over base center and concentration parameters
    double BaseLogPrior() const {
        return -baseconcentration / Naa;
    }

    //! log prior of amino-acid fitness profiles
    double AALogPrior() const { return componentaafitnessarray->GetLogProb(); }

    //! log prior of amino-acid fitness profile k
    double AALogPrior(int k) const { return componentaafitnessarray->GetLogProb(k); }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return componentpathsuffstatbidimarray->GetLogProb(*codonmatrices);
    }

    //! return log prob of the substitution mappings over sites allocated to
    //! component k of the mixture
    double PathSuffStatLogProb(int k) const {
        return componentpathsuffstatbidimarray->GetRowLogProb(k,*codonmatrices);
    }

    //! return log prob of current branch lengths, as a function of branch lengths
    //! hyperparameter lambda
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    //! return log prob of first-level mixture components (i.e. all amino-acid
    //! profiles drawn from component k of the base distribution), as a function
    //! of the center and concentration parameters of this component
    double BaseSuffStatLogProb() const {
        return basesuffstat.GetLogProb(basecenter, baseconcentration);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! log prob factor to be recomputed when moving branch lengths hyperparameter
    //! lambda
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    //! log prob factor to be recomputed when moving nucleotide mutation rate
    //! parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + PathSuffStatLogProb(); }

    //! log prob factor to be recomputed when moving aa hyper params (center and
    //! concentration) for component k of base distribution
    double BaseLogProb() const { return BaseLogPrior() + BaseSuffStatLogProb(); }

    //! log prob factor to be recomputed when moving kappa, the concentration
    //! parameter of the first-level strick breaking process
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    double GLogProb() const { return GLogPrior() + PathSuffStatLogProb(); }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics if substitution mappings across sites
    void CollectSitePathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    //! gather site-specific sufficient statistics component-wise
    void CollectComponentPathSuffStat() {
        componentpathsuffstatbidimarray->Clear();
        componentpathsuffstatbidimarray->Add(*sitepathsuffstatarray, *sitealloc, *Galloc);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //-------------------
    //  Moves
    //-------------------

    //! complete MCMC move schedule
    double Move() override {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        // UpdateMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            if (blmode < 2) {
                MoveBranchLengths();
            }

            CollectSitePathSuffStat();
            CollectComponentPathSuffStat();

            if (nucmode < 2) {
                MoveNucRates();
            }

            if (omegamode < 2) {
                MoveOmega();
            }

            if (Gcat > 1)   {
                MoveGinvshape();
                ResampleGAlloc();
                CollectComponentPathSuffStat();
            }

            MoveAAMixture(3);

            MoveBase();

        }
    }

    double MoveGinvshape()  {
        ScalingMove(Ginvshape, 1.0, 3, &AAMutSelGSBDPOmegaModel::GLogProb, &AAMutSelGSBDPOmegaModel::UpdateSelAC, this);
        ScalingMove(Ginvshape, 0.1, 3, &AAMutSelGSBDPOmegaModel::GLogProb, &AAMutSelGSBDPOmegaModel::UpdateSelAC, this);
        return 1.0;
    }

    //! Gibbs resample mixture allocations
    void ResampleGAlloc() {
        vector<double> postprob(Gcat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetGAllocPostProb(i, postprob);
            Galloc->GibbsResample(i, postprob);
        }
    }

    //! get allocation posterior probabilities for a given site
    void GetGAllocPostProb(int site, vector<double> &postprob) {
        double max = 0;
        const vector<double> &w = Gweight;
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        int k = sitealloc->GetVal(site);
        for (int i = 0; i < Gcat; i++) {
            double tmp = suffstat.GetLogProb(codonmatrices->GetVal(k,i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) {
                max = tmp;
            }
        }

        double total = 0;
        for (int i = 0; i < Gcat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < Gcat; i++) {
            postprob[i] /= total;
        }
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MCMC move schedule on branch lengths
    void MoveBranchLengths() {
        ResampleBranchLengths();
        if (blmode == 0) {
            MoveLambda();
        }
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &AAMutSelGSBDPOmegaModel::LambdaHyperLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &AAMutSelGSBDPOmegaModel::LambdaHyperLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    //! MH move on omega
    void MoveOmega() {
        omegapathsuffstat.Clear();
        omegapathsuffstat.AddSuffStat(*codonmatrices, *componentpathsuffstatbidimarray);
        if (omegaprior == 0) {
            GibbsResampleOmega();
        } else {
            MultipleTryMoveOmega(100);
            if (omega != 1.0)   {
                MHMoveOmega(10, 1.0);
                MHMoveOmega(10, 0.1);
            }
        }
        UpdateCodonMatrices();
    }

    void GibbsResampleOmega() {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(),
                                    beta + omegapathsuffstat.GetBeta());
    }

    double OmegaLogProb(double dposom) const {
        return OmegaLogPrior() + omegapathsuffstat.GetLogProb(1.0 + dposom);
    }

    double MHMoveOmega(int nrep, double tuning)    {
        double dposom = omega - 1;
        double nacc = 0;
        double ntot = 0;
        for (int rep=0; rep<nrep; rep++)    {
            double deltalogprob = -OmegaLogProb(dposom);
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            dposom *= e;
            deltalogprob += OmegaLogProb(dposom);
            deltalogprob += m;
            int acc = (log(Random::Uniform()) < deltalogprob)
                && ((!maxdposom) || (dposom < maxdposom));
            if (acc)    {
                nacc++;
            }
            else    {
                dposom /= e;
            }
            ntot++;
        }
        omega = dposom + 1;
        return nacc/ntot;
    }

    double DrawPosOm() {
        double ret = 0;
        if (omegaprior == 0)    {
            cerr << "error: in draw pos om but not under mixture model\n";
            exit(1);
        }
        else if (omegaprior == 1)   {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            do  {
                ret = Random::Gamma(alpha, beta);
            } while ((maxdposom > 0) && (ret > maxdposom));
            if (!ret) {
                 ret = 1e-5;
            }
        }
        else if (omegaprior == 2)   {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            do  {
                ret = exp(Random::Gamma(alpha, beta)) - 1;
            } while ((maxdposom > 0) && (ret > maxdposom));
            if (!ret) {
                 ret = 1e-5;
            }
        }
        else if (omegaprior == 3)   {
            double gamma = 1.0 / dposomhyperinvshape;
            do  {
                ret = gamma * tan(Pi * Random::Uniform() / 2);
            } while ((maxdposom > 0) && (ret > maxdposom));
        }
        return ret;
    }

    int MultipleTryMoveOmega(int ntry) {

        double logp0 = omegapathsuffstat.GetLogProb(1.0);

        vector<double> logparray(ntry, 0);
        vector<double> dposomarray(ntry, 0);
        double max = 0;
        for (int i = 0; i < ntry; i++) {
            if ((!i) && (omega > 1.0)) {
                dposomarray[i] = omega - 1.0;
            } else {
                dposomarray[i] = DrawPosOm();
                /*
                */
            }
            logparray[i] = omegapathsuffstat.GetLogProb(1.0 + dposomarray[i]);
            if ((!i) || (max < logparray[i])) {
                max = logparray[i];
            }
        }

        double tot = 0;
        vector<double> cumulprob(ntry, 0);
        for (int i = 0; i < ntry; i++) {
            tot += exp(logparray[i] - max);
            cumulprob[i] = tot;
        }
        double logp1 = log(tot/ntry) + max;

        int accept = 0;
        int choice = 0;

        if (omega == 1.0) {
            double logratio = log(dposompi) + logp1 - log(1.0 - dposompi) - logp0;
            if (log(Random::Uniform()) < logratio) {
                choice = 1;
                accept = 1;
            } else {
                choice = 0;
            }
        } else {
            double logratio = log(1.0 - dposompi) + logp0 - log(dposompi) - logp1;
            if (log(Random::Uniform()) < logratio) {
                choice = 0;
                accept = 1;
            } else {
                choice = 1;
            }
        }

        if (choice == 1) {
            // randomly choose dposom among the ntry values in array
            double u = tot * Random::Uniform();
            int i = 0;
            while ((i < ntry) && (u > cumulprob[i])) {
                i++;
            }
            if (i == ntry) {
                cerr << "error in MultipleTryMoveOmega: overflow\n";
                exit(1);
            }
            omega = 1.0 + dposomarray[i];
        } else {
            omega = 1.0;
        }

        return accept;
    }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        ProfileMove(nucrelrate, 0.1, 1, 3, &AAMutSelGSBDPOmegaModel::NucRatesLogProb,
                    &AAMutSelGSBDPOmegaModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &AAMutSelGSBDPOmegaModel::NucRatesLogProb,
                    &AAMutSelGSBDPOmegaModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &AAMutSelGSBDPOmegaModel::NucRatesLogProb,
                    &AAMutSelGSBDPOmegaModel::UpdateMatrices, this);

        ProfileMove(nucstat, 0.1, 1, 3, &AAMutSelGSBDPOmegaModel::NucRatesLogProb,
                    &AAMutSelGSBDPOmegaModel::UpdateMatrices, this);
        ProfileMove(nucstat, 0.01, 1, 3, &AAMutSelGSBDPOmegaModel::NucRatesLogProb,
                    &AAMutSelGSBDPOmegaModel::UpdateMatrices, this);
    }

    //! MCMC module for the mixture amino-acid fitness profiles
    void MoveAAMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {

            MoveAAProfiles();
            ResampleEmptyComponents();

            ResampleAlloc();
            LabelSwitchingMove();
            CollectComponentPathSuffStat();

            ResampleWeights();
            MoveKappa();
        }
    }

    //! resample empty components of the mixture from prior
    void ResampleEmptyComponents() {
        componentaafitnessarray->PriorResample(*occupancy);
        for (int i=0; i<occupancy->GetSize(); i++)   {
            if (!occupancy->GetVal(i))  {
                UpdateCodonMatrix(i);
            }
        }
    }

    //! MH move on amino-acid fitness profiles (occupied components only)
    void MoveAAProfiles() {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double CompMoveAAProfiles(int nrep) {
        MoveAA(1.0, 1, nrep);
        MoveAA(1.0,3,nrep);
        MoveAA(0.3,3,nrep);
        MoveAA(0.1, 3, nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MulMoveAAProfiles(int nrep) {
        MoveAAGamma(3.0, nrep);
        MoveAAGamma(1.0, nrep);
        MoveAAGamma(0.3,nrep);
        MoveAAGamma(0.1,nrep);
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double MoveAA(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Naa];
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                vector<double> &aa = (*componentaafitnessarray)[i];
                for (int rep = 0; rep < nrep; rep++) {
                    for (int l = 0; l < Naa; l++) {
                        bk[l] = aa[l];
                    }
                    double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                    deltalogprob += loghastings;
                    UpdateCodonMatrix(i);
                    deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        for (int l = 0; l < Naa; l++) {
                            aa[l] = bk[l];
                        }
                        UpdateCodonMatrix(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! helper function: log density of 20 gammas
    double GammaAALogPrior(const vector<double> &x, const vector<double> &aacenter, double aaconc) {
        double total = 0;
        for (int l = 0; l < Naa; l++) {
            total += (aaconc * aacenter[l] - 1) * log(x[l]) - x[l] -
                     Random::logGamma(aaconc * aacenter[l]);
        }
        return total;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MoveAAGamma(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                double aaconc = baseconcentration;
                const vector<double> &aacenter = basecenter;

                vector<double> &aa = (*componentaafitnessarray)[i];
                vector<double> x(Naa, 0);
                double z = Random::sGamma(aaconc);
                for (int l = 0; l < Naa; l++) {
                    x[l] = z * aa[l];
                }

                double bkz = z;
                vector<double> bkx = x;
                vector<double> bkaa = aa;

                for (int rep = 0; rep < nrep; rep++) {
                    double deltalogprob =
                        -GammaAALogPrior(x, aacenter, aaconc) - PathSuffStatLogProb(i);

                    double loghastings = 0;
                    z = 0;
                    for (int l = 0; l < Naa; l++) {
                        double m = tuning * (Random::Uniform() - 0.5);
                        double e = exp(m);
                        x[l] *= e;
                        z += x[l];
                        loghastings += m;
                    }
                    for (int l = 0; l < Naa; l++) {
                        aa[l] = x[l] / z;
                        if (aa[l] < 1e-50) {
                            aa[l] = 1e-50;
                        }
                    }

                    deltalogprob += loghastings;

                    UpdateCodonMatrix(i);

                    deltalogprob += GammaAALogPrior(x, aacenter, aaconc) + PathSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                        bkaa = aa;
                        bkx = x;
                        bkz = z;
                    } else {
                        aa = bkaa;
                        x = bkx;
                        z = bkz;
                        UpdateCodonMatrix(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! Gibbs resample mixture allocations
    void ResampleAlloc() {
        vector<double> postprob(Ncat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetAllocPostProb(i, postprob);
            sitealloc->GibbsResample(i, postprob);
        }
        UpdateOccupancies();
    }

    //! update mixture occupancy suff stats (for resampling mixture weights)
    void UpdateOccupancies() {
        occupancy->Clear();
        occupancy->AddSuffStat(*sitealloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetAllocPostProb(int site, vector<double> &postprob) {
        double max = 0;
        const vector<double> &w = weight->GetArray();
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        int g = Galloc->GetVal(site);
        for (int i = 0; i < Ncat; i++) {
            double tmp = suffstat.GetLogProb(codonmatrices->GetVal(i,g));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) {
                max = tmp;
            }
        }

        double total = 0;
        for (int i = 0; i < Ncat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < Ncat; i++) {
            postprob[i] /= total;
        }
    }

    //! MCMC sequence for label switching moves
    void LabelSwitchingMove() {
        Permutation permut(Ncat);
        weight->LabelSwitchingMove(5, *occupancy, permut);
        sitealloc->Permute(permut);
        componentaafitnessarray->Permute(permut);
        UpdateCodonMatrices(permut);
    }

    //! Gibbs resample mixture weights (based on occupancy suff stats)
    void ResampleWeights() { weight->GibbsResample(*occupancy); }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa() {
        ScalingMove(kappa, 1.0, 10, &AAMutSelGSBDPOmegaModel::StickBreakingHyperLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);
        ScalingMove(kappa, 0.3, 10, &AAMutSelGSBDPOmegaModel::StickBreakingHyperLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);
        weight->SetKappa(kappa);
    }

    //! MCMC module for moving the center and concentration parameters of the
    //! components of the the base mixture
    void MoveBase() {

        basesuffstat.Clear();
        componentaafitnessarray->AddSuffStat(basesuffstat, *occupancy);

        ProfileMove(basecenter, 1.0, 1, 3, &AAMutSelGSBDPOmegaModel::BaseLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);
        ProfileMove(basecenter, 1.0, 3, 3, &AAMutSelGSBDPOmegaModel::BaseLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);
        ProfileMove(basecenter, 0.3, 1, 3, &AAMutSelGSBDPOmegaModel::BaseLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);

        ScalingMove(baseconcentration, 1.0, 10, &AAMutSelGSBDPOmegaModel::BaseLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);
        ScalingMove(baseconcentration, 0.3, 10, &AAMutSelGSBDPOmegaModel::BaseLogProb,
                    &AAMutSelGSBDPOmegaModel::NoUpdate, this);

        componentaafitnessarray->SetConcentration(baseconcentration);
        ResampleEmptyComponents();
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    //! return number of occupied components in first-level mixture (mixture of
    //! amino-acid fitness profiles)
    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                n++;
            }
        }
        return n;
    }

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const { return Random::GetEntropy(nucrelrate); }

    double GetPredictedDNDS() const  {
        double m0 = 0;
        double m1 = 0;
        // vector<vector<double>> dnds(Ncat, vector<double>(Gcat,0));
        for (int i=0; i<Ncat; i++) {
            for (int j=0; j<Gcat; j++)  {
                m0 += weight->GetVal(i) * Gweight[j];
                m1 += weight->GetVal(i) * Gweight[j] * codonmatrices->GetVal(i,j).GetPredictedDNDS();
            }
        }
        return m1 / m0;
    }

    //! return mean entropy of amino-acd fitness profiles
    double GetMeanRefAAEntropy() const { return componentaafitnessarray->GetMeanEntropy(); }

    double GetMeanAAEntropy() const {
        double m1 = 0;
        double m0 = 0;
        for (int i=0; i<Ncat; i++)   {
            for (int j=0; j<Gcat; j++)  {
                m0 += weight->GetVal(i) * Gweight[j];
                m1 += weight->GetVal(i) * Gweight[j] * Random::GetEntropy(selacprofiles->GetVal(i,j));
            }
        }
        return m1 / m0;
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "omega\t";
        os << "dnds\t";
        if (Gcat > 1)   {
            os << "ginvshape\t";
        }
        os << "ncluster\t";
        os << "kappa\t";
        os << "refaaent\t";
        os << "aaent\t";
        os << "aacenterent\t";
        os << "aaconc\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << omega << '\t';
        os << GetPredictedDNDS() << '\t';
        if (Gcat > 1)   {
            os << Ginvshape << '\t';
        }
        os << GetNcluster() << '\t';
        os << kappa << '\t';
        os << GetMeanRefAAEntropy() << '\t';
        os << GetMeanAAEntropy() << '\t';
        os << Random::GetEntropy(basecenter) << '\t';
        os << baseconcentration << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const override {
    }

    void FromStream(istream &is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2) {
            is >> nucrelrate;
            is >> nucstat;
        }
        is >> basecenter;
        is >> baseconcentration;
        is >> kappa;
        weight->FromStreamSB(is);
        // is >> *weight;
        is >> *componentaafitnessarray;
        is >> *sitealloc;
        if (omegamode < 2) {
            is >> omega;
        }
        if (Gcat > 1)   {
            is >> Ginvshape;
            is >> *Galloc;
        }
    }

    void ToStream(ostream &os) const override {
        if (blmode < 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        }
        if (nucmode < 2) {
            os << nucrelrate << '\t';
            os << nucstat << '\t';
        }
        os << basecenter << '\t';
        os << baseconcentration << '\t';
        os << kappa << '\t';
        weight->ToStreamSB(os);
        // os << *weight << '\t';
        os << *componentaafitnessarray << '\t';
        os << *sitealloc << '\t';
        if (omegamode < 2) {
            os << omega << '\t';
        }
        if (Gcat > 1)   {
            os << Ginvshape << '\t';
            os << *Galloc << '\t';
        }
    }

    //! return size of model, when put into an MPI buffer (in multigene context --
    //! only omegatree)
    unsigned int GetMPISize() const {
        int size = 0;
        if (blmode < 2) {
            size++;
            size += branchlength->GetMPISize();
        }
        if (nucmode < 2) {
            size += nucrelrate.size();
            size += nucstat.size();
        }
        size += Naa + 1;
        size++;
        size += weight->GetMPISizeSB();
        size += componentaafitnessarray->GetMPISize();
        size += sitealloc->GetMPISize();
        if (omegamode < 2) {
            size++;
        }
        if (Gcat > 1)   {
            size ++;
            size += Gcat;
        }
        return size;
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer &is) {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2) {
            is >> nucrelrate;
            is >> nucstat;
        }
        is >> basecenter;
        is >> baseconcentration;
        is >> kappa;
        weight->MPIGetSB(is);
        is >> *componentaafitnessarray;
        is >> *sitealloc;
        if (omegamode < 2) {
            is >> omega;
        }
        if (Gcat > 1)   {
            is >> Ginvshape;
            is >> *Galloc;
        }
    }

    //! write array into MPI buffer
    void MPIPut(MPIBuffer &os) const {
        if (blmode < 2) {
            os << lambda;
            os << *branchlength;
        }
        if (nucmode < 2) {
            os << nucrelrate;
            os << nucstat;
        }
        os << basecenter;
        os << baseconcentration;
        os << kappa;
        weight->MPIPutSB(os);
        os << *componentaafitnessarray;
        os << *sitealloc;
        if (omegamode < 2) {
            os << omega;
        }
        if (Gcat > 1)   {
            os << Ginvshape;
            os << *Galloc;
        }
    }
};
