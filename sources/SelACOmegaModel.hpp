
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
 * - an omega multiplier, capturing deviations of the non-syn rate from the
 * model (see Rodrigue and Lartillot, 2107); this parameter is fixed to 1 by
 * default.
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 * - omega: fixed to 1 or exponential of rate 1
 *
 * In a multi-gene context, shrinkage across genes can be applied to branch
 * lengths, omega, nucleotide rate parameters (rho and pi)
 *
 */

class SelACOmegaModel : public ProbModel {

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

    double wcom;
    double wpol;
    double wvol;
    // distance matrix between amino acids 
    vector<double> aadist;
    // number of categories of the discretized gamma
    int Gcat;
    double Ginvshape;
    vector<double> G;

    // stringency (level of expression)
    double psihypermean;
    double psihyperinvshape;
    double psi;

    SelACProfileBidimArray* selacprofiles;
    AADiffSelCodonMatrixBidimArray* codonmatrices;

    // site allocation to each of the 20 amino-acid categories
    vector<double> aaweight;
    vector<double> Gweight;
    MultinomialAllocationVector* alloc;
    MultinomialAllocationVector* Galloc;
    OccupancySuffStat* aaoccupancy;

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
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
    int aadistmode;

    // 0: grantham
    // 1: general distance
    int aadistmodel;

    double dposompi;
    double dposomhypermean;
    double dposomhyperinvshape;

    // 0: simple gamma prior
    // 1: mix of (1-pi) at 1 and pi at 1+d, with d ~
    // Gamma(dposomhypermean,dposomhyperinvshape) 2: mix of 2 (modal) gamma
    // distributions: one at 1 and another one with mean > 1
    int omegaprior;

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
    //! - inomegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    //! - Gcat: number of bins of discretized gamma for G's across sites

    SelACOmegaModel(string datafile, string treefile, int inaadistmodel, int inaadistmode, int inomegamode, int inomegaprior,
                            int inGcat) {
        blmode = 0;
        nucmode = 0;
        aadistmodel = inaadistmodel;
        aadistmode = inaadistmode;
        omegamode = inomegamode;
        omegaprior = inomegaprior;

        Gcat = inGcat;

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

    SelACOmegaModel(const CodonSequenceAlignment* incodondata, const Tree* intree, int inaadistmodel, int inaadistmode, int inomegamode, int inomegaprior,
                            int inGcat) {
        blmode = 0;
        nucmode = 0;
        aadistmodel = inaadistmodel;
        aadistmode = inaadistmode;
        omegamode = inomegamode;
        omegaprior = inomegaprior;

        Gcat = inGcat;

        codondata = incodondata;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();

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

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 3: fixed to 1
    //! - mode == 2: shared and estimated across genes: currently not implemented
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, either mode 3 and mode 0 can be used -- default
    //! mode is 3.
    void SetOmegaMode(int mode) { omegamode = mode; }

    //! \brief set prior for omega
    //
    //! 0: simple gamma prior
    //! 1: mix of (1-pi) at 1 and pi at 1+d, with d ~
    //! Gamma(dposomhypermean,dposomhyperinvshape)
    void SetOmegaPrior(int inprior) { omegaprior = inprior; }

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

        // global omega (fixed to 1 by default)
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        omega = 1.0;
        dposompi = 0.1;
        dposomhypermean = 1;
        dposomhyperinvshape = 0.5;

        wcom = grantham_wcom;
        wpol = grantham_wpol;
        wvol = grantham_wvol;

        aadist.assign(Naarr, 1.0);
        if (!aadistmodel)    {
            UpdateGrantham();
        }

        Ginvshape = 1.0;
        G.assign(Gcat, 1.0);
        UpdateG();
        psi = 0.5;
        psihypermean = 10.0;
        psihyperinvshape = 1;

        selacprofiles = new SelACProfileBidimArray(aadist,G,psi);

        aaweight.assign(Naa, 1.0/Naa);
        Gweight.assign(Gcat, 1.0/Gcat);
        alloc = new MultinomialAllocationVector(Nsite,aaweight);
        Galloc = new MultinomialAllocationVector(Nsite,Gweight);
        aaoccupancy = new OccupancySuffStat(Naa);

        codonmatrices = new AADiffSelCodonMatrixBidimArray(*selacprofiles, *GetCodonStateSpace(), *nucmatrix, omega);
        sitesubmatrixarray = new DoubleMixtureSelector<SubMatrix>(codonmatrices, alloc, Galloc);

        // create polyprocess

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, sitesubmatrixarray);
        phyloprocess->Unfold();

        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatbidimarray = new PathSuffStatBidimArray(Naa,Gcat);
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

    void SetW(double inwcom, double inwpol) {
        wcom = inwcom;
        wpol = inwpol;
        UpdateSelAC();
    }

    void SetAADist(const vector<double>& inaadist)  {
        aadist = inaadist;
        UpdateSelAC();
    }

    void SetAAWeight(const vector<double>& inaaweight)  {
        aaweight = inaaweight;
        UpdateSelAC();
    }

    void SetG(double inGinvshape)   {
        Ginvshape = inGinvshape;
        UpdateSelAC();
    }

    void SetPsi(double inpsi)   {
        psi = inpsi;
        UpdateSelAC();
    }

    void SetPsiHyperParameters(double inpsihypermean, double inpsihyperinvshape) {
        psihypermean = inpsihypermean;
        psihyperinvshape = inpsihyperinvshape;
    }

    const OccupancySuffStat *GetOccupancies() const { return aaoccupancy; }

    double GetPsi() const {
        return psi;
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

    void UpdateG()  {
        Random::DiscGamma(G,1.0/Ginvshape);
    }

    void UpdateSelAC()  {
        UpdateG();
        if (! aadistmodel)   {
            UpdateGrantham();
        }
        selacprofiles->SetPsi(psi);
        selacprofiles->Update();
        UpdateCodonMatrices();
    }

    void UpdateAA(int aa)   {
        selacprofiles->UpdateRow(aa);
        codonmatrices->CorruptRow(aa);
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatrices() {
        codonmatrices->SetOmega(omega);
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
        UpdateOccupancies();
        UpdateMatrices();
        ResampleSub(1.0);
    }

    void PostPred(string name) override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        UpdateSelAC();
        UpdateOccupancies();
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
            // psi, G, aadist
            total += PsiHyperLogPrior();
            total += AALogPrior();
            // allocation of sites to aa
            // total += alloc->GetLogProb();
        }

        // omega
        if (omegamode < 2) {
            total += OmegaLogPrior();
        }

        return total;
    }

    double PsiHyperLogPrior() const {
        if (psihypermean == 0)  {
            return -log(psi);
        }
        double alpha = 1.0 / psihyperinvshape;
        double beta = alpha / psihypermean;
        return alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(psi) - beta * psi;
    }

    double AALogPrior() const {
        double total = 0;
        total -= Ginvshape;
        total += PsiHyperLogPrior();
        if (! aadistmodel)   {
            total -= log(wcom) + log(wpol);
        }
        else    {
            for (int i=0; i<Naarr; i++) {
                total -= aadist[i];
            }
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

    //! log prior over omega (gamma of mean omegahypermean and inverse shape
    //! omegahyperinvshape)
    double OmegaLogPrior() const {
        double ret = 0;
        if (omegaprior == 0) {
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            ret = alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(omega) -
                  beta * omega;
        } else if (omegaprior == 1) {
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

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return componentpathsuffstatbidimarray->GetLogProb(*codonmatrices);
    }

    double AAPathSuffStatLogProb(int aa) const  {
        return componentpathsuffstatbidimarray->GetRowLogProb(aa,*codonmatrices);
    }

    void AAPathSuffStatLogProbs(vector<double>& logprobvector) const    {
        return componentpathsuffstatbidimarray->GetRowLogProbs(logprobvector,*codonmatrices);
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

    double AALogProb() const { return AALogPrior() + PathSuffStatLogProb(); }

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

            if (omegamode < 2) {
                MoveOmega();
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
                MoveAADist();
                AAPsiCompMove(1.0, 10);
            }
            if (Gcat > 1)   {
                MoveGinvshape();
            }
        }

        MovePsi();

        ResampleAlloc();
        ResampleGAlloc();

        if (aadistmode < 2) {
            CollectComponentPathSuffStat();
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
        ScalingMove(lambda, 1.0, 10, &SelACOmegaModel::LambdaHyperLogProb,
                    &SelACOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &SelACOmegaModel::LambdaHyperLogProb,
                    &SelACOmegaModel::NoUpdate, this);
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
        }
        UpdateCodonMatrices();
    }

    void GibbsResampleOmega() {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(),
                                    beta + omegapathsuffstat.GetBeta());
    }

    int MultipleTryMoveOmega(int ntry) {
        /*
        if (omega == 1.0)   {

            // initial prob: (1-pi)*p(D | omega == 1.0)
            // ntry values of dposom
            // final prob: pi* < p(D | omega == 1.0 + dposom) >
        }

        else    {

            // ntry values of dposom, with first value being equal to current omega
        - 1.0
            // initial prob: pi* < p(D | omega == 1.0 + dposom) >
            // final prob: (1-pi)*p(D | omega == 1.0)
        }
        // MH decision rule between omega == 1.0 or omega > 1.0
        // if decision -> omega > 1.0, then randomly sample dposom among ntry values
        and set omega = 1.0 + dposom
        */

        double logp0 = omegapathsuffstat.GetLogProb(1.0);

        double alpha = 1.0 / dposomhyperinvshape;
        double beta = alpha / dposomhypermean;

        vector<double> logparray(ntry, 0);
        vector<double> dposomarray(ntry, 0);
        double max = 0;
        for (int i = 0; i < ntry; i++) {
            if ((!i) && (omega > 1.0)) {
                dposomarray[i] = omega - 1.0;
            } else {
                dposomarray[i] = Random::Gamma(alpha, beta);
                if (!dposomarray[i]) {
                    dposomarray[i] = 1e-5;
                }
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

        // Multiple-Try Gibbs version
        /*
        double m = (logp0>logp1) ? logp0 : logp1;
        double q0 = (1-dposompi) * exp(logp0-m);
        double q1 = dposompi * exp(logp1-m);
        double p0 = q0 / (q0+q1);

        if (Random::Uniform() < p0) {
            omega = 1.0;
        }
        else    {
            // randomly choose dposom among the ntry values in array
        }
        */

        // Multiple-Try MH version
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
        ProfileMove(nucrelrate, 0.1, 1, 3, &SelACOmegaModel::NucRatesLogProb,
                    &SelACOmegaModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &SelACOmegaModel::NucRatesLogProb,
                    &SelACOmegaModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &SelACOmegaModel::NucRatesLogProb,
                    &SelACOmegaModel::UpdateMatrices, this);

        ProfileMove(nucstat, 0.1, 1, 3, &SelACOmegaModel::NucRatesLogProb,
                    &SelACOmegaModel::UpdateMatrices, this);
        ProfileMove(nucstat, 0.01, 1, 3, &SelACOmegaModel::NucRatesLogProb,
                    &SelACOmegaModel::UpdateMatrices, this);
    }

    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Naa - i - 1) * i / 2 + j - i - 1
                       : (2 * Naa - j - 1) * j / 2 + i - j - 1;
    }

    void MoveAADist() {
        aadist_acc[0] += MoveAAPair(1.0, 10);
        aadist_tot[0] ++;
        aadist_acc[1] += MoveAAPair(0.3, 10);
        aadist_tot[1] ++;
    }

    double AAPsiCompMove(double tuning, int nrep)    {

        double nacc = 0;
        double ntot = 0;

        for (int rep=0; rep<nrep; rep++)    {
            double delta = - AALogPrior();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            psi /= e;
            for (int i=0; i<Naarr; i++) {
                aadist[i] *= e;
            }
            delta += AALogPrior();
            delta += (Naarr-1)*m;

            int accepted = (log(Random::Uniform()) < delta);
            if (accepted) {
                nacc++;
            } else {
                psi *= e;
                for (int i=0; i<Naarr; i++) {
                    aadist[i] /= e;
                }
            }
            ntot++;
        }

        UpdateSelAC();

        return nacc / ntot;
    }

    double MoveAAPair(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        for (int rep=0; rep<nrep; rep++)    {
            int aa1 = (int) (Naa * Random::Uniform());
            int aa2 = (int) ((Naa-1) * Random::Uniform());
            if (aa2 >= aa1) aa2++;
            int i = rrindex(aa1,aa2);
            double bk = aadist[i];

            double delta = - AAPathSuffStatLogProb(aa1) - AAPathSuffStatLogProb(aa2);
            delta += aadist[i];
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            aadist[i] *= e;
            UpdateAA(aa1);
            UpdateAA(aa2);
            delta += AAPathSuffStatLogProb(aa1) + AAPathSuffStatLogProb(aa2);
            delta -= aadist[i];
            delta += m;

            int accepted = (log(Random::Uniform()) < delta);
            if (accepted) {
                nacc++;
            } else {
                aadist[i] = bk;
                UpdateAA(aa1);
                UpdateAA(aa2);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    /*
    double MoveAADist() {
        aadist_acc[0] += ProfileMove(aadist, 1.00, 1, 10, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[0] ++;
        aadist_acc[1] += ProfileMove(aadist, 0.30, 3, 10, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[1] ++;
        aadist_acc[2] += ProfileMove(aadist, 0.10, 3, 10, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[2] ++;
        return 1.0;
    }
    */

    double MoveGranthamWeights()    {
        aadist_acc[0] += ScalingMove(wcom, 1.0, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[0] ++;
        /*
        aadist_acc[1] += ScalingMove(wcom, 0.3, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[1] ++;
        aadist_acc[2] += ScalingMove(wcom, 0.1, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[2] ++;
        */
        aadist_acc[0] += ScalingMove(wpol, 1.0, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[0] ++;
        /*
        aadist_acc[1] += ScalingMove(wpol, 0.3, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[1] ++;
        aadist_acc[2] += ScalingMove(wpol, 0.1, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        aadist_tot[2] ++;
        */
        return 1.0;
    }

    double MoveGinvshape()  {
        G_acc[0] += ScalingMove(Ginvshape, 1.0, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        G_tot[0] ++;
        G_acc[1] += ScalingMove(Ginvshape, 0.1, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        G_tot[1] ++;
        return 1.0;
    }

    double MovePsi()  {
        psi_acc[0] += ScalingMove(psi, 1.0, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
        psi_tot[0] ++;
        psi_acc[1] += ScalingMove(psi, 0.1, 3, &SelACOmegaModel::AALogProb, &SelACOmegaModel::UpdateSelAC, this);
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
        // UpdateOccupancies();
    }

    //! Gibbs resample mixture allocations
    void ResampleGAlloc() {
        vector<double> postprob(Gcat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetGAllocPostProb(i, postprob);
            Galloc->GibbsResample(i, postprob);
        }
    }

    //! update mixture occupancy suff stats (for resampling mixture weights)
    void UpdateOccupancies() {
        aaoccupancy->Clear();
        aaoccupancy->AddSuffStat(*alloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetAllocPostProb(int site, vector<double> &postprob) {
        double max = 0;
        const vector<double> &w = aaweight;
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
        UpdateOccupancies();
        double total = 0;
        for (int i=0; i<Naa; i++)   {
            aaweight[i] = Random::sGamma(1 + aaoccupancy->GetVal(i));
            total += aaweight[i];
        }
        for (int i=0; i<Naa; i++)   {
            aaweight[i] /= total;
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
                m0 += aaweight[i] * Gweight[j];
                m1 += aaweight[i] * Gweight[j] * codonmatrices->GetVal(i,j).GetPredictedDNDS();
            }
        }
        return m1 / m0;
    }

    double GetMeanAAEntropy() const {
        double m1 = 0;
        double m0 = 0;
        for (int i=0; i<Naa; i++)   {
            for (int j=0; j<Gcat; j++)  {
                m0 += aaweight[i] * Gweight[j];
                m1 += aaweight[i] * Gweight[j] * Random::GetEntropy(selacprofiles->GetVal(i,j));
            }
        }
        return m1 / m0;
    }

    double GetMeanAADist() const    {
        double m1 = 0;
        for (int i=0; i<Naarr; i++) {
            m1 += aadist[i];
        }
        return m1 / Naarr;
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
        return m2;
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "psi\t";
        if (Gcat > 1)   {
            os << "ginvshape\t";
        }
        if (omegamode < 2)  {
            os << "omega\t";
        }
        os << "dnds\t";
        if (! aadistmodel)   {
            os << "w_comp\t";
            os << "w_pol\t";
        }
        os << "distmean\t";
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
        os << psi << '\t';
        if (Gcat > 1)   {
            os << Ginvshape << '\t';
        }
        if (omegamode < 2)  {
            os << omega << '\t';
        }
        os << GetPredictedDNDS() << '\t';
        if (! aadistmodel)   {
            os << wcom << '\t';
            os << wpol << '\t';
        }
        os << GetMeanAADist() << '\t';
        os << GetVarAADist() << '\t';
        os << Random::GetEntropy(aaweight) << '\t';
        os << GetMeanAAEntropy() << '\t';

        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
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

        if (aadistmode < 2) {
            is >> psi;
            is >> Ginvshape;
            if (! aadistmodel)   {
                is >> wcom >> wpol;
            }
            else    {
                is >> aadist;
            }
            is >> aaweight;
            is >> *alloc;
            is >> *Galloc;
        }

        if (omegamode < 2) {
            is >> omega;
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

        if (aadistmode < 2) {
            os << psi << '\t';
            os << Ginvshape << '\t';
            if (! aadistmodel)   {
                os << wcom << '\t' << wpol << '\t';
            }
            else    {
                os << aadist << '\t';
            }
            os << aaweight << '\t';
            os << *alloc << '\t';
            os << *Galloc << '\t';
        }

        if (omegamode < 2) {
            os << omega << '\t';
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
        if (aadistmode < 2) {
            size++;
            size++;
            if (! aadistmodel)   {
                size += 2;
            }
            else    {
                size += aadist.size();
            }
            size += aaweight.size();
            size += alloc->GetMPISize();
            size += Galloc->GetMPISize();
        }
        if (omegamode < 2) {
            size++;
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

        if (aadistmode < 2) {
            is >> psi;
            is >> Ginvshape;
            if (! aadistmodel)   {
                is >> wcom >> wpol;
            }
            else    {
                is >> aadist;
            }
            is >> aaweight;
            is >> *alloc;
            is >> *Galloc;
        }

        if (omegamode < 2) {
            is >> omega;
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

        if (aadistmode < 2) {
            os << psi;
            os << Ginvshape;
            if (! aadistmodel)   {
                os << wcom << wpol;
            }
            else    {
                os << aadist;
            }
            os << aaweight;
            os << *alloc;
            os << *Galloc;
        }

        if (omegamode < 2) {
            os << omega;
        }
    }
};
