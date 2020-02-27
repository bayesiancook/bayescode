
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

#include "grantham.hpp"
#include "AADiffSelCodonMatrixBidimArray.hpp"
#include "SelACProfileBidimArray.hpp"

/**
 * \brief The mutation-selection model with constant fitness landscape over the
 * tree -- SelAC version
 *
 * The model is parameterized by
 * - an fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for
 * j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process
 * (assumed homogeneous across sites and lineages)
 * - an array of site-specific amino-acid fitness profiles F_ia, for site i and
 * amino-acid a
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 *
 * In a multi-gene context, shrinkage across genes can be applied to branch
 * lengths, nucleotide rate parameters (rho and pi)
 *
 */

class DiscSelACModel : public ProbModel {

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

    double wcom;
    double wpol;
    double wvol;
    // distance matrix between amino acids 
    vector<double> aadist;

    // number of categories of the discretized gamma
    // a discretized array of regularly spaced X values
    // with associated weights

    // this array is constant
    // x_i = xmin + i/Gcat (xmax - mxin)
    int Gcat;
    double xmin;
    double xmax;
    vector<double> G;

    // stringency (level of expression)
    double logpsihypermean;
    double logpsihypervar;
    double logpsi;

    double Gvar;
    double Gvarhypermean;
    double Gvarhyperinvshape;

    // this array is variable, based on psi and Gvar
    vector<double> Gweight;

    SelACProfileBidimArray* selacprofiles;
    AADiffSelCodonMatrixBidimArray* codonmatrices;

    // site allocation to each of the 20 amino-acid categories
    vector<double> aaweights;
    MultinomialAllocationVector* alloc;
    MultinomialAllocationVector* Galloc;
    OccupancySuffStat* aaoccupancy;
    OccupancySuffStat* Goccupancy;

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
    DoubleMixtureSelector<SubMatrix> *sitesubmatrixarray;

    PhyloProcess *phyloprocess;

    PathSuffStatArray *sitepathsuffstatarray;
    PathSuffStatBidimArray *componentpathsuffstatbidimarray;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    int blmode;
    int nucmode;
    int aadistmode;
    int gvarmode;
    int psimode;

    // 0: grantham
    // 1: general distance
    int aadistmodel;

    Chrono aachrono;
    Chrono totchrono;

    vector<double> aadist_acc;
    vector<double> aadist_tot;
    vector<double> G_acc;
    vector<double> G_tot;
    vector<double> psi_acc;
    vector<double> psi_tot;

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
    //! - Gcat: number of bins of discretized gamma for G's across sites

    DiscSelACModel(string datafile, string treefile, int inaadistmodel, int inGcat, double inxmin, double inxmax) {
        blmode = 0;
        nucmode = 0;
        aadistmodel = inaadistmodel;

        aadistmode = 0;
        gvarmode = 0;
        psimode = 0;

        Gcat = inGcat;
        xmin = inxmin;
        xmax = inxmax;

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

        // Allocate();
    }

    DiscSelACModel(const CodonSequenceAlignment* incodondata, const Tree* intree, int inaadistmodel, int inGcat, double inxmin, double inxmax) {
        blmode = 0;
        nucmode = 0;
        aadistmodel = inaadistmodel;

        aadistmode = 0;
        gvarmode = 0;
        psimode = 0;

        codondata = incodondata;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();

        Gcat = inGcat;
        xmin = inxmin;
        xmax = inxmax;

        // Allocate();
    }

    //! \brief set estimation method for branch lengths
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //! \brief set estimation method for branch lengths
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetBLMode(int mode) { blmode = mode; }

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetNucMode(int mode) { nucmode = mode; }

    void SetAAMode(int inaadistmode, int ingvarmode, int inpsimode) {
        aadistmode = inaadistmode;
        gvarmode = ingvarmode;
        psimode = inpsimode;
    }

    //! allocate the model (data structures)
    void Allocate() {

        aadist_acc.assign(3,0);
        aadist_tot.assign(3,0.1);
        G_acc.assign(3,0);
        G_tot.assign(3,0.1);
        psi_acc.assign(3,0);
        psi_tot.assign(3,0.1);

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

        wcom = grantham_wcom;
        wpol = grantham_wpol;
        wvol = grantham_wvol;

        aadist.assign(Naarr, 1.0 / Naarr);
        Random::DirichletSample(aadist, vector<double>(Naarr, 1.0 / Naarr), ((double)Naarr));
        if (!aadistmodel)    {
            UpdateGrantham();
        }

        G.assign(Gcat, 1.0);
        for (int i=0; i<Gcat; i++)  {
            G[i] = exp(xmin + ((double) i)/Gcat*(xmax-xmin));
        }

        Gvar = 1.0;
        Gvarhypermean = 10.0;
        Gvarhyperinvshape = 1.0;

        logpsi = 0;
        logpsihypermean = 0;
        logpsihypervar = 10.0;

        Gweight.assign(Gcat,0);
        UpdateGweights();

        // with a Naarr multiplier <=> aadist is 1 on average
        selacprofiles = new SelACProfileBidimArray(aadist, G, ((double) Naarr));

        Galloc = new MultinomialAllocationVector(Nsite,Gweight);
        Goccupancy = new OccupancySuffStat(Gcat);

        aaweights.assign(Naa, 1.0/Naa);
        alloc = new MultinomialAllocationVector(Nsite,aaweights);
        aaoccupancy = new OccupancySuffStat(Naa);

        codonmatrices = new AADiffSelCodonMatrixBidimArray(*selacprofiles, *GetCodonStateSpace(), *nucmatrix, 1.0);
        sitesubmatrixarray = new DoubleMixtureSelector<SubMatrix>(codonmatrices, alloc, Galloc);

        // create polyprocess

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, sitesubmatrixarray);
        phyloprocess->Unfold();

        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatbidimarray = new PathSuffStatBidimArray(Naa, Gcat, GetCodonStateSpace()->GetNstate());
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

    void SetW(double inwcom, double inwpol) {
        wcom = inwcom;
        wpol = inwpol;
        UpdateAADist();
    }

    void SetAADist(const vector<double>& inaadist)  {
        aadist = inaadist;
        UpdateAADist();
    }

    void SetAAWeights(const vector<double>& inaaweights)  {
        aaweights = inaaweights;
        UpdateAADist();
    }

    void SetGvar(double inGvar)   {
        Gvar = inGvar;
        UpdateGweights();
    }

    void SetGvarHyperParameters(double inGvarhypermean, double inGvarhyperinvshape) {
        Gvarhypermean = inGvarhypermean;
        Gvarhyperinvshape = inGvarhyperinvshape;
    }

    void SetLogPsi(double inlogpsi)   {
        logpsi = inlogpsi;
        UpdateGweights();
    }

    void SetLogPsiHyperParameters(double inlogpsihypermean, double inlogpsihypervar) {
        logpsihypermean = inlogpsihypermean;
        logpsihypervar = inlogpsihypervar;
    }

    const OccupancySuffStat *GetAAOccupancies() const { return aaoccupancy; }

    const OccupancySuffStat *GetGOccupancies() const { return Goccupancy; }

    double GetLogPsi() const    {
        return logpsi;
    }

    double GetGvar() const  {
        return Gvar;
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

    const PathSuffStatBidimArray& GetComponentPathSuffStat() const { return *componentpathsuffstatbidimarray; }

    //! copy nucleotide rates into vectors given as arguments (multi-gene
    //! analyses)
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void UpdateGrantham()   {
        double tot = 0;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                double tcom = grantham_com[b] - grantham_com[a];
                double tpol = grantham_pol[b] - grantham_pol[a];
                double tvol = grantham_vol[b] - grantham_vol[a];
                double d = sqrt(wcom*tcom*tcom + wpol*tpol*tpol + wvol*tvol*tvol);
                aadist[rrindex(a,b)] = d;
                tot += d;
            }
        }
        tot /= Naarr;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                aadist[rrindex(a,b)] /= tot;
            }
        }
    }

    void UpdateGweights()  {
        double totw = 0;
        for (int i=0; i<Gcat; i++)  {
            double x = xmin + ((double) i)/Gcat * (xmax - xmin);
            double w = exp( - (x - logpsi)*(x-logpsi) / 2 / Gvar);
            totw += w;
            Gweight[i] = w;
        }
        for (int i=0; i<Gcat; i++)  {
            Gweight[i] /= totw;
        }
    }

    void UpdateAADist() {
        if (! aadistmodel)   {
            UpdateGrantham();
        }
        selacprofiles->Update();
        UpdateCodonMatrices();
    }

    void UpdateSelAC()  {

        UpdateGweights();
        UpdateGOccupancies();

        UpdateAADist();
        UpdateAAOccupancies();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatrices() {
        codonmatrices->Corrupt();
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
        UpdateSelAC();
        UpdateMatrices();
        ResampleSub(1.0);
    }

    void PostPred(string name) override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
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

        // branch lengths
        if (blmode < 2) {
            total += BranchLengthsLogPrior();
        }

        // nuc rates
        if (nucmode < 2) {
            total += NucRatesLogPrior();
        }

        if (aadistmode < 2) {
            total += AADistLogPrior();
        }

        if (gvarmode < 2)   {
            total += GvarHyperLogPrior();
        }

        if (psimode < 2)    {
            total += LogPsiHyperLogPrior();
        }

        return total;
    }

    double GvarHyperLogPrior() const   {
        double alpha = 1.0 / Gvarhyperinvshape;
        double beta = alpha / Gvarhypermean;
        return alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(Gvar) - beta * Gvar;
    }

    double LogPsiHyperLogPrior() const {
        return - (logpsi - logpsihypermean) * (logpsi - logpsihypermean) / 2 / logpsihypervar;
    }

    double AADistLogPrior() const {
        double total = 0;
        if (! aadistmodel)   {
            total -= log(wcom) + log(wpol);
        }
        return total;
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

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

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate, nucrelratehypercenter,
                                             1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return componentpathsuffstatbidimarray->GetLogProb(*codonmatrices);
    }

    double PathSuffStatLogProb(const vector<int>& mask) const  {
        double tot = 0;
        for (int a=0; a<Naa; a++)   {
            if (mask[a])   {
                tot += componentpathsuffstatbidimarray->GetRowLogProb(a,*codonmatrices);
            }
        }
        return tot;
    }

    double GAllocSuffStatLogProb() const {
        return Goccupancy->GetLogProb(Gweight);
    }

    //! return log prob of current branch lengths, as a function of branch lengths
    //! hyperparameter lambda
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
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

    double GvarLogProb() const { return GvarHyperLogPrior() + GAllocSuffStatLogProb(); }
    double PsiLogProb() const { return LogPsiHyperLogPrior() + GAllocSuffStatLogProb(); }
    double AADistLogProb() const { return AADistLogPrior() + PathSuffStatLogProb(); }

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
        componentpathsuffstatbidimarray->Add(*sitepathsuffstatarray, *alloc, *Galloc);
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
        UpdateMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            totchrono.Start();
            if (blmode < 2) {
                MoveBranchLengths();
            }

            CollectSitePathSuffStat();
            CollectComponentPathSuffStat();

            if (nucmode < 2) {
                MoveNucRates();
            }

            aachrono.Start();
            MoveSelAC();
            aachrono.Stop();

            totchrono.Stop();
        }
    }


    void MoveSelAC()    {

        if (aadistmode < 2) {
            if (! aadistmodel)   {
                MoveGranthamWeights();
            }
            else    {
                // MoveAADistOld();
                MoveAADist();
            }
        }

        if (gvarmode < 2)   {
            MoveGvar();
        }

        if (psimode < 2)    {
            MoveLogPsi();
        }

        ResampleAlloc();
        ResampleGAlloc();
        CollectComponentPathSuffStat();

        if (aadistmode < 2) {
            ResampleWeights();
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
        ScalingMove(lambda, 1.0, 10, &DiscSelACModel::LambdaHyperLogProb,
                    &DiscSelACModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &DiscSelACModel::LambdaHyperLogProb,
                    &DiscSelACModel::NoUpdate, this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        ProfileMove(nucrelrate, 0.1, 1, 3, &DiscSelACModel::NucRatesLogProb,
                    &DiscSelACModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &DiscSelACModel::NucRatesLogProb,
                    &DiscSelACModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &DiscSelACModel::NucRatesLogProb,
                    &DiscSelACModel::UpdateMatrices, this);

        ProfileMove(nucstat, 0.1, 1, 3, &DiscSelACModel::NucRatesLogProb,
                    &DiscSelACModel::UpdateMatrices, this);
        ProfileMove(nucstat, 0.01, 1, 3, &DiscSelACModel::NucRatesLogProb,
                    &DiscSelACModel::UpdateMatrices, this);
    }

    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Naa - i - 1) * i / 2 + j - i - 1
                       : (2 * Naa - j - 1) * j / 2 + i - j - 1;
    }

    void getaapair(int index, int& a, int& b)   {
        int range = Naa-1;
        int k = 0;
        a = 0;
        while (index >= k + range)   {
            k += range;
            range--;
            a++;
        }
        b = index - k + a + 1;
        if (rrindex(a,b) != index)  {
            cerr << "error in getaapair\n";
            exit(1);
        }
        if ((b < 1) || (b > 19) || (a > 18))  {
            cerr << "aa out of range : " << a << '\t' << b << '\n';
            exit(1);
        }
        /*
        cerr << index << '\t' << a << '\t' << b << '\n';
        cerr << rrindex(a,b) << '\n';
        cerr << '\n';
        */
    }

    double MoveAADistOld() {
        aadist_acc[0] += ProfileMove(aadist, 1.00, 1, 10, &DiscSelACModel::AADistLogProb, &DiscSelACModel::UpdateAADist, this);
        aadist_tot[0] ++;
        aadist_acc[1] += ProfileMove(aadist, 0.30, 3, 10, &DiscSelACModel::AADistLogProb, &DiscSelACModel::UpdateAADist, this);
        aadist_tot[1] ++;
        aadist_acc[2] += ProfileMove(aadist, 0.10, 3, 10, &DiscSelACModel::AADistLogProb, &DiscSelACModel::UpdateAADist, this);
        aadist_tot[2] ++;
        return 1.0;
    }

    double MoveAADist() {
        aadist_acc[0] += MoveAADist(10, 1, 1.0);
        aadist_tot[0] ++;
        aadist_acc[1] += MoveAADist(10, 3, 0.3);
        aadist_tot[1] ++;
        aadist_acc[2] += MoveAADist(10, 3, 0.1);
        aadist_tot[2] ++;
        return 1.0;
    }

    double MoveAADist(int nrep, int n, double tuning) {

        vector<double> bk(Naarr,0);
        int indices[2*n];
        double nacc = 0;

        for (int rep=0; rep<nrep; rep++)    {

            Random::DrawFromUrn(indices, 2 * n, Naarr);

            vector<int> mask(Naa,0);
            for (int i=0; i<2*n; i++)   {
                int a,b;
                getaapair(indices[i],a,b);
                mask[a] = 1;
                mask[b] = 1;
            }

            double deltalogprob = - PathSuffStatLogProb(mask);

            for (int i=0; i<2*n; i++) {
                bk[indices[i]] = aadist[indices[i]];
            }

            for (int i = 0; i < n; i++) {
                int i1 = indices[2 * i];
                int i2 = indices[2 * i + 1];
                double tot = aadist[i1] + aadist[i2];
                double x = aadist[i1];

                double h = tot * tuning * (Random::Uniform() - 0.5);
                x += h;
                while ((x < 0) || (x > tot)) {
                    if (x < 0) {
                        x = -x;
                    }
                    if (x > tot) {
                        x = 2 * tot - x;
                    }
                }
                aadist[i1] = x;
                aadist[i2] = tot - x;
            }

            UpdateAADist();
            deltalogprob += PathSuffStatLogProb(mask);

            int accept = (log(Random::Uniform()) < deltalogprob);
            if (accept) {
                nacc++;
            }
            else    {
                for (int i=0; i<2*n; i++) {
                    aadist[indices[i]] = bk[indices[i]];
                }
                UpdateAADist();
            }
        }
        return nacc / nrep;
    }

    double MoveGranthamWeights()    {
        aadist_acc[0] += ScalingMove(wcom, 1.0, 10, &DiscSelACModel::AADistLogProb, &DiscSelACModel::UpdateAADist, this);
        aadist_tot[0] ++;
        aadist_acc[0] += ScalingMove(wpol, 1.0, 10, &DiscSelACModel::AADistLogProb, &DiscSelACModel::UpdateAADist, this);
        aadist_tot[0] ++;
        return 1.0;
    }

    double MoveGvar()  {
        G_acc[0] += ScalingMove(Gvar, 1.0, 10, &DiscSelACModel::GvarLogProb, &DiscSelACModel::UpdateGweights, this);
        G_tot[0] ++;
        G_acc[1] += ScalingMove(Gvar, 0.1, 10, &DiscSelACModel::GvarLogProb, &DiscSelACModel::UpdateGweights, this);
        G_tot[1] ++;
        return 1.0;
    }

    double MoveLogPsi()  {
        psi_acc[0] += SlidingMove(logpsi, 1.0, 10, 0, 0, &DiscSelACModel::PsiLogProb, &DiscSelACModel::UpdateGweights, this);
        psi_tot[0] ++;
        psi_acc[1] += SlidingMove(logpsi, 0.1, 10, 0, 0, &DiscSelACModel::PsiLogProb, &DiscSelACModel::UpdateGweights, this);
        psi_tot[1] ++;
        return 1.0;
    }

    //! Gibbs resample mixture allocations
    void ResampleAlloc() {
        vector<double> postprob(Naa, 0);
        for (int i = 0; i < Nsite; i++) {
            GetAllocPostProb(i, postprob);
            alloc->GibbsResample(i, postprob);
        }
        UpdateAAOccupancies();
    }

    //! Gibbs resample mixture allocations
    void ResampleGAlloc() {
        vector<double> postprob(Gcat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetGAllocPostProb(i, postprob);
            Galloc->GibbsResample(i, postprob);
        }
        UpdateGOccupancies();
    }

    //! update mixture occupancy suff stats (for resampling mixture weights)
    void UpdateAAOccupancies() {
        aaoccupancy->Clear();
        aaoccupancy->AddSuffStat(*alloc);
    }

    void UpdateGOccupancies()   {
        Goccupancy->Clear();
        Goccupancy->AddSuffStat(*Galloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetAllocPostProb(int site, vector<double> &postprob) {
        double max = 0;
        const vector<double> &w = aaweights;
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        int g = Galloc->GetVal(site);
        for (int i = 0; i < Naa; i++) {
            double tmp = suffstat.GetLogProb(codonmatrices->GetVal(i,g));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) {
                max = tmp;
            }
        }

        double total = 0;
        for (int i = 0; i < Naa; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < Naa; i++) {
            postprob[i] /= total;
        }
    }

    //! get allocation posterior probabilities for a given site
    void GetGAllocPostProb(int site, vector<double> &postprob) {
        double max = 0;
        const vector<double> &w = Gweight;
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        int aa = alloc->GetVal(site);
        for (int i = 0; i < Gcat; i++) {
            double tmp = suffstat.GetLogProb(codonmatrices->GetVal(aa,i));
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

    //! Gibbs resample mixture weights (based on occupancy suff stats)
    void ResampleWeights() { 
        // UpdateAAOccupancies();
        double total = 0;
        for (int i=0; i<Naa; i++)   {
            aaweights[i] = Random::sGamma(1 + aaoccupancy->GetVal(i));
            total += aaweights[i];
        }
        for (int i=0; i<Naa; i++)   {
            aaweights[i] /= total;
        }
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const { return Random::GetEntropy(nucrelrate); }

    double GetPredictedDNDS() const  {
        double m0 = 0;
        double m1 = 0;
        vector<vector<double>> dnds(Naa, vector<double>(Gcat,0));
        for (int i=0; i<Naa; i++) {
            for (int j=0; j<Gcat; j++)  {
                m0 += aaweights[i] * Gweight[j];
                m1 += aaweights[i] * Gweight[j] * codonmatrices->GetVal(i,j).GetPredictedDNDS();
            }
        }
        return m1 / m0;
    }

    double GetMeanAAEntropy() const {
        double m1 = 0;
        double m0 = 0;
        for (int i=0; i<Naa; i++)   {
            for (int j=0; j<Gcat; j++)  {
                m0 += aaweights[i] * Gweight[j];
                m1 += aaweights[i] * Gweight[j] * Random::GetEntropy(selacprofiles->GetVal(i,j));
            }
        }
        return m1 / m0;
    }

    double GetVarAADist() const {
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<Naarr; i++) {
            m1 += aadist[i];
            m2 += aadist[i]*aadist[i];
        }
        m1 /= Naarr;
        m2 /= Naarr;
        m2 -= m1*m1;
        m2 *= Naarr*Naarr;
        return m2;
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "psi\t";
        os << "gvar\t";
        os << "dnds\t";
        if (! aadistmodel)   {
            os << "w_comp\t";
            os << "w_pol\t";
        }
        os << "distvar\t";
        os << "weightent\t";
        os << "aaent\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << 3 * branchlength->GetTotalLength() << '\t';
        os << exp(logpsi) << '\t';
        os << Gvar << '\t';
        os << GetPredictedDNDS() << '\t';
        if (! aadistmodel)   {
            os << wcom << '\t';
            os << wpol << '\t';
        }
        os << GetVarAADist() << '\t';
        os << Random::GetEntropy(aaweights) << '\t';
        os << GetMeanAAEntropy() << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void TraceGweights(ostream& os) const {
        os << Gweight << '\n';
    }

    void Monitor(ostream &os) const override {
        os << totchrono.GetTime() << '\t' << aachrono.GetTime() << '\n';
        os << "prop time in aa moves  : " << aachrono.GetTime() / totchrono.GetTime() << '\n';
        os << "aadist\n";
        for (size_t i=0; i<aadist_acc.size(); i++)  {
            os << aadist_acc[i] / aadist_tot[i] << '\n';
        }
        os << "G\n";
        for (size_t i=0; i<G_acc.size(); i++)  {
            os << G_acc[i] / G_tot[i] << '\n';
        }
        os << "psi\n";
        for (size_t i=0; i<psi_acc.size(); i++)  {
            os << psi_acc[i] / psi_tot[i] << '\n';
        }
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

        if (gvarmode < 2)   {
            is >> Gvar;
        }

        if (psimode < 2)    {
            is >> logpsi;
        }

        if (aadistmode < 2) {
            if (! aadistmodel)   {
                is >> wcom >> wpol;
            }
            else    {
                is >> aadist;
            }
            is >> aaweights;
        }

        is >> *alloc;
        is >> *Galloc;
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

        if (gvarmode < 2)   {
            os << Gvar << '\t';
        }

        if (psimode < 2)    {
            os << logpsi << '\t';
        }

        if (aadistmode < 2) {
            if (! aadistmodel)   {
                os << wcom << '\t' << wpol << '\t';
            }
            else    {
                os << aadist << '\t';
            }
            os << aaweights << '\t';
        }

        os << *alloc << '\t';
        os << *Galloc << '\t';
    }

    //! return size of model, when put into an MPI buffer
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

        if (gvarmode < 2)   {
            size++;
        }

        if (psimode < 2)    {
            size++;
        }

        if (aadistmode < 2) {
            if (! aadistmodel)   {
                size += 2;
            }
            else    {
                size += aadist.size();
            }
            size += aaweights.size();
        }

        size += alloc->GetMPISize();
        size += Galloc->GetMPISize();

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

        if (gvarmode < 2)   {
            is >> Gvar;
        }

        if (psimode < 2)    {
            is >> logpsi;
        }

        if (aadistmode < 2) {
            if (! aadistmodel)   {
                is >> wcom >> wpol;
            }
            else    {
                is >> aadist;
            }
            is >> aaweights;
        }

        is >> *alloc;
        is >> *Galloc;
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

        if (gvarmode < 2)   {
            os << Gvar;
        }

        if (psimode < 2)    {
            os << logpsi;
        }

        if (aadistmode < 2) {
            if (! aadistmodel)   {
                os << wcom << wpol;
            }
            else    {
                os << aadist;
            }
            os << aaweights;
        }

        os << *alloc;
        os << *Galloc;
    }
};

