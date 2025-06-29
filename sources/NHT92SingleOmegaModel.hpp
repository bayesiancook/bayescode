#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSuffStat.hpp"
#include "T92SubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "IIDBeta.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

#include "dSOmegaPathSuffStat.hpp"
#include "GCConsdSOmegaPathSuffStat.hpp"

#include "T92SubMatrixBranchArray.hpp"
#include "NucPathSuffStatBranchArray.hpp"
#include "NHCodonSubMatrixBranchArray.hpp"

class SingleOmegaModel : public ProbModel {
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
    int fixbl;

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

    double rootkappa;
    double kappahypermean;
    double kappahyperinvshape;
    BranchIIDGamma* branchkappa;

    GammaSuffStat kappahypersuffstat;

    double rootgamma;
    double gammahypermean;
    double gammahyperinvconc;
    BranchIIDBeta* branchgamma;

    BetaSuffStat gammahypersuffstat;

    T92SubMatrix* rootnucmatrix;
    T92SubMatrixBranchArray* nucmatrixarray;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    // includes root suffstat in addition to branch suff stats
    NucPathSuffStatBranchArray* nucpathsuffstat;

    // Omega

    double omegahypermean;
    double omegahyperinvshape;
    double omega;

    BranchHomogeneousSelector<double>* omegaselector;

    // codon matrices across branches (with same omega but different nuc matrices)
    NHMGOmegaCodonSubMatrixBranchArray *codonmatrixarray;
    // root matrix
    MGOmegaCodonSubMatrix* rootcodonmatrix;

    // PhyloProcess

    PhyloProcess *phyloprocess;

    // suff stats for substitution paths
    // summed over all branches and over all sites
    // PathSuffStat pathsuffstat;
    PathSuffStatNodeArray *pathsuffstatbrancharray;

    // or, alternatively, collected as a simple Poisson suff stat, as a function
    // of omega
    OmegaPathSuffStat omegapathsuffstat;


  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    SingleOmegaModel(string datafile, string treefile) {

        blmode = 0;
        fixbl = 0;
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
    }

    SingleOmegaModel(const CodonSequenceAlignment* incodondata, const Tree* intree) {

        blmode = 0;
        fixbl = 0;
        nucmode = 0;

        codondata = incodondata;
        Nsite = codondata->GetNsite();
        Ntaxa = codondata->GetNtaxa();
        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();
    }

    int GetNsite() const    {
        return Nsite;
    }

    int GetNbranch() const  {
        return Nbranch;
    }

    const Tree* GetTree() const {
        return tree;
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

        if (fixbl)  {
            for (int j=0; j<Nbranch; j++)  {
                (*branchlength)[j] = tree->GetBranchLength(j);
            }
        }

        // Nucleotide rates

        rootkappa = 2.0;
        kappahypermean = 2.0;
        kappahyperinvshape = 0.1;
        double kappaalpha = 1.0 / kappahyperinvshape;
        double kappabeta = kappaalpha / kappahypermean;
        branchkappa = new BranchIIDGamma(*tree, kappaalpha, kappabeta);
        branchkappa->SetAllBranches(2.0);

        rootgamma = 0.5;
        gammahypermean = 0.5;
        gammahyperinvconc = 0.1;
        double gammaalpha = gammahypermean / gammahyperinvconc;
        double gammabeta = (1-gammahypermean) / gammahyperinvconc;
        branchgamma = new BranchIIDBeta(*tree, gammaalpha, gammabeta);
        branchgamma->SetAllBranches(0.5);

        rootnucmatrix = new T92SubMatrix(rootkappa, rootgamma, true);
        nucmatrixarray = new T92SubMatrixBranchArray(branchkappa, branchgamma, true);

        nucpathsuffstat = new NucPathSuffStatBranchArray(*tree);

        // Omega

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        omega = 0.2;

        omegaselector = new BranchHomogeneousSelector<double>(tree, omega);

        rootcodonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), rootnucmatrix, omega);
        codonmatrixarray = new NHMGOmegaCodonSubMatrixBranchArray(GetCodonStateSpace(), nucmatrixarray, omegaselector);

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, codonmatrixarray, rootcodonmatrix);
        phyloprocess->Unfold();

        pathsuffstatbrancharray = new PathSuffStatNodeArray(*tree);
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return current value of omega
    double GetOmega() const { return omega; }

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

    void SetFixBL(int in)   {
        fixbl = in;
        if (fixbl)  {
            blmode = 2;
        }
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

    /*
    //! set nucleotide rates (relative exchangeabilities and eq. frequencies) to a
    //! new value (multi-gene analyses)
    void SetNucRates(double inrootkappa, double inrootgamma,
            const std::vector<double>& inkapparray, const std::vector<double>& inbranchgamma)    {
        rootkappa = inrootkappa;
        rootgamma = inrootgamma;
        branchkappa = inbranchkappa;
        branchgamma = inbranchgamma;
        TouchMatrices();
    }

    //! get a copy of nucleotide rates into arrays given as arguments
    void GetNucRates(double& inrootkappa, double& inrootgamma,
            std::vector<double>& inkapparray, std::vector<double>& inbranchgamma)    {

        inrootkappa = rootkappa;
        inrootgamma = rootgamma;
        inbranchkappa = branchkappa;
        inbranchgamma = branchgamma;
    }

    void SetNucRatesHyperParameters(double inkappamean, double inkappainvshape,
            double ingammamean, double ingammainvshape) {
        kappahypermean = inkappahypermean;
        kappahyperinvshape = inkappahyperinvshape;
        gammahypermean = ingammahypermean;
        gammahyperinvconc = ingammahyperinvconc;

        double kappaalpha = 1.0 / kappahyperinvshape;
        double kappabeta = kappaalpha / kappahypermean;
        branchkappa->SetShape(kappaalpha);
        branchkappa->SetScale(kappabeta);

        branchgamma->SetMeanInvConc(gammahypermean, gammahyperinvconc);
    }
    */

    // Omega

    //! \brief set omega to a new value
    //!
    //! Used in a multigene context.
    //! Notifies corruption to the codon matrix.
    void SetOmega(double inomega) {
        omega = inomega;
        TouchCodonMatrices();
    }

    //! \brief set the hyperparameters of the gamma prior over omega
    //!
    //! Used in a multigene context.
    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape) {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchNucMatrices() {
        TouchRootNucMatrix();
        nucmatrixarray->UpdateMatrices();
    }

    void TouchRootNucMatrix()	{
        rootnucmatrix->SetKappa(rootkappa);
        rootnucmatrix->SetGC(rootgamma);
        rootnucmatrix->CorruptMatrix();
    }

    void TouchBranchNucMatrix(int branchindex)  {
        (*nucmatrixarray)[branchindex].SetKappa(branchkappa->GetVal(branchindex));
        (*nucmatrixarray)[branchindex].SetGC(branchgamma->GetVal(branchindex));
        (*nucmatrixarray)[branchindex].CorruptMatrix();
    }

    //! \brief tell the codon matrix that its parameters have changed and that it
    //! should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchCodonMatrices() {
        rootcodonmatrix->SetOmega(omega);
        rootcodonmatrix->CorruptMatrix();
        codonmatrixarray->UpdateCodonMatrices();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices() {
        TouchNucMatrices();
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

        double kappaalpha = 1.0 / kappahyperinvshape;
        double kappabeta = kappaalpha / kappahypermean;
        branchkappa->SetShape(kappaalpha);
        branchkappa->SetScale(kappabeta);

        branchgamma->SetMeanInvConc(gammahypermean, gammahyperinvconc);

        TouchMatrices();
        ResampleSub(1.0);
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

        double kappaalpha = 1.0 / kappahyperinvshape;
        double kappabeta = kappaalpha / kappahypermean;
        branchkappa->SetShape(kappaalpha);
        branchkappa->SetScale(kappabeta);

        branchgamma->SetMeanInvConc(gammahypermean, gammahyperinvconc);

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
        total += OmegaLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    // Branch lengths

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

    double NucRatesHyperLogPrior() const    {
        double total = 0;
        total -= gammahyperinvconc;
        total -= kappahypermean / 5.0;
        total -= kappahyperinvshape;
        return total;
    }

    double KappaHyperLogPrior() const   {
        return -kappahypermean / 5.0 - kappahyperinvshape;
    }

    double GammaHyperLogPrior() const    {
        return -gammahyperinvconc;
    }

    double NucRatesLogPrior() const {
        double total = 0;

        // uniform prior over rootgamma
        // root kappa is fixed

        total += branchkappa->GetLogProb();
        total += branchgamma->GetLogProb();
        return total;
    }

    double RootNucRatesLogPrior() const {
        return 0;
    }

    double NucRatesLogPrior(int branchindex) const  {
        return branchkappa->GetLogProb(branchindex) + branchgamma->GetLogProb(branchindex);
    }

    // Omega

    double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return Random::logGammaDensity(omega, alpha, beta);
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    // Branch lengths

    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    // Nuc Rates

    double GammaHyperSuffStatLogProb() const {
        return gammahypersuffstat.GetMeanInvConcLogProb(gammahypermean, gammahyperinvconc);
    }

    double GammaHyperLogProb() const    {
        return GammaHyperLogPrior() + GammaHyperSuffStatLogProb();
    }

    double KappaHyperSuffStatLogProb() const    {
        double kappaalpha = 1.0 / kappahyperinvshape;
        double kappabeta = kappaalpha / kappahypermean;
        return kappahypersuffstat.GetLogProb(kappaalpha, kappabeta);
    }

    double KappaHyperLogProb() const    {
        return KappaHyperLogPrior() + KappaHyperSuffStatLogProb();
    }

    // const NucPathSuffStatBranchArray &GetNucPathSuffStatBranchArray() const { return *nucpathsuffstat; }

    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat->Clear();
        nucpathsuffstat->AddSuffStat(*codonmatrixarray, *rootcodonmatrix, *pathsuffstatbrancharray);
    }

    double NucRatesSuffStatLogProb(int branchindex) const   {
        return nucpathsuffstat->GetVal(branchindex).GetLogProb(nucmatrixarray->GetVal(branchindex), *GetCodonStateSpace());
    }

    double RootNucRatesSuffStatLogProb() const  {
        return nucpathsuffstat->GetRootVal().GetLogProb(*rootnucmatrix, *GetCodonStateSpace());
    }


    double NucRatesLogProb(int branchindex) const   {
        return NucRatesLogPrior(branchindex) + NucRatesSuffStatLogProb(branchindex);
    }

    double RootNucRatesLogProb() const  {
        return RootNucRatesLogPrior() + RootNucRatesSuffStatLogProb();
    }


    // Paths

    void CollectPathSuffStat() {
        pathsuffstatbrancharray->Clear();
        pathsuffstatbrancharray->AddSuffStat(*phyloprocess);
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrix. Both pathsuffstat and
    //! codonmatrix are assumed to be updated.
    double PathSuffStatLogProb() const { return pathsuffstatbrancharray->GetLogProb(*codonmatrixarray, *rootcodonmatrix); }

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

            if (!FixedNucRates()) {
                TouchMatrices();
                MoveNucRates();
                MoveNucRatesHyperParameters();
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
        ScalingMove(lambda, 1.0, 10, &SingleOmegaModel::LambdaHyperLogProb, &SingleOmegaModel::NoUpdate,
                    this);
        ScalingMove(lambda, 0.3, 10, &SingleOmegaModel::LambdaHyperLogProb, &SingleOmegaModel::NoUpdate,
                    this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    // Nucleotide rates

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        SlidingMove(rootgamma, 0.3, 3, 0, 1, &SingleOmegaModel::RootNucRatesLogProb,
                & SingleOmegaModel::TouchRootNucMatrix, this);
        SlidingMove(rootgamma, 0.1, 3, 0, 1, &SingleOmegaModel::RootNucRatesLogProb,
                & SingleOmegaModel::TouchRootNucMatrix, this);

        for (int j=0; j<GetNbranch(); j++)  {
            MoveBranchGamma(j, 0.3, 3);
            MoveBranchKappa(j, 0.3, 3);
            MoveBranchGamma(j, 0.1, 3);
            MoveBranchKappa(j, 0.1, 3);
        }

        // important for codon matrices (nucl. matrices are already updated)
        TouchMatrices();
    }

    double MoveBranchGamma(int branch, double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        double& x = (*branchgamma)[branch];
        double min = 0;
        double max = 1.0;

        for (int rep = 0; rep < nrep; rep++) {
            double bk = x;
            double deltalogprob = - NucRatesLogProb(branch);
            double m = tuning * (Random::Uniform() - 0.5);
            x += m;
            if (max > min) {
                while ((x < min) || (x > max)) {
                    if (x < min) {
                        x = 2 * min - x;
                    }
                    if (x > max) {
                        x = 2 * max - x;
                    }
                }
            }
            TouchBranchNucMatrix(branch);
            deltalogprob += NucRatesLogProb(branch);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x = bk;
                TouchBranchNucMatrix(branch);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveBranchKappa(int branch, double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        double& x = (*branchkappa)[branch];

        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -NucRatesLogProb(branch);
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            x *= e;
            TouchBranchNucMatrix(branch);
            deltalogprob += NucRatesLogProb(branch);
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x /= e;
                TouchBranchNucMatrix(branch);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    void MoveNucRatesHyperParameters()  {
        MoveKappaHyperParameters();
        MoveGammaHyperParameters();
    }

    void MoveKappaHyperParameters() {
        kappahypersuffstat.Clear();
        kappahypersuffstat.AddSuffStat(*branchkappa);
        ScalingMove(kappahypermean, 1.0, 10, &SingleOmegaModel::KappaHyperLogProb,
                    &SingleOmegaModel::NoUpdate, this);
        ScalingMove(kappahyperinvshape, 1.0, 10, &SingleOmegaModel::KappaHyperLogProb,
                    &SingleOmegaModel::NoUpdate, this);
        double alpha = 1.0 / kappahyperinvshape;
        double beta = alpha / kappahypermean;
        branchkappa->SetShape(alpha);
        branchkappa->SetScale(beta);
    }

    void MoveGammaHyperParameters() {
        gammahypersuffstat.Clear();
        branchgamma->AddSuffStat(gammahypersuffstat);
        ScalingMove(gammahypermean, 1.0, 10, &SingleOmegaModel::GammaHyperLogProb,
                    &SingleOmegaModel::NoUpdate, this);
        ScalingMove(gammahyperinvconc, 1.0, 10, &SingleOmegaModel::GammaHyperLogProb,
                    &SingleOmegaModel::NoUpdate, this);
        branchgamma->SetMeanInvConc(gammahypermean, gammahyperinvconc);
    }

    // Omega

    //! Gibbs resample omega (based on sufficient statistics of current
    //! substitution mapping)
    void MoveOmega() {
        omegapathsuffstat.Clear();
        omegapathsuffstat.AddSuffStat(*codonmatrixarray, *pathsuffstatbrancharray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(),
                                    beta + omegapathsuffstat.GetBeta());
        TouchCodonMatrices();
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    /*
    void AdddSOmegaPathSuffStat(dSOmegaPathSuffStatArray& into) const {
        PathSuffStatArray pathsuffstatarray(Nsite);
        pathsuffstatarray.Clear();
        pathsuffstatarray.AddSuffStat(*phyloprocess);
        double totlength = branchlength->GetTotalLength();
        into.AddSuffStat(
               [this](int i) -> const MGOmegaCodonSubMatrix& {return *codonmatrix;},
               [&path = pathsuffstatarray](int i) {return path.GetVal(i);},
               [&l = totlength](int i) {return l;},
               [&om = omega](int i) {return om;});
    }

    void AdddSOmegaPathSuffStat(dSOmegaPathSuffStatBranchArray& into) const {
        PathSuffStatNodeArray pathsuffstatarray(*tree);
        pathsuffstatarray.Clear();
        pathsuffstatarray.AddSuffStat(*phyloprocess);
        into.AddSuffStat(*codonmatrix, pathsuffstatarray, *branchlength, omega);
    }

    void AddGCConsdSOmegaPathSuffStat(GCConsdSOmegaPathSuffStatBranchArray& into) const {
        PathSuffStatNodeArray pathsuffstatarray(*tree);
        pathsuffstatarray.Clear();
        pathsuffstatarray.AddSuffStat(*phyloprocess);
        into.AddSuffStat(*codonmatrix, pathsuffstatarray, *branchlength, omega);
    }

    //! collect generic sufficient statistics from substitution mappings
    void AddNodePathSuffStat(RelativePathSuffStatNodeArray& into) const {
        PathSuffStatNodeArray pathsuffstatarray(*tree, codondata->GetNstate());
        pathsuffstatarray.AddSuffStat(*phyloprocess);

        RelativePathSuffStatNodeArray relpathsuffstatarray(*tree, codondata->GetNstate());
        relpathsuffstatarray.Clear();
        relpathsuffstatarray.AddSuffStat(pathsuffstatarray, *branchlength);

        into.Add(relpathsuffstatarray);
    }

    void AddDoubleCounts(vector<double>& meancounts) const  {
        vector<vector<int>> counts(GetNbranch(), vector<int>(GetNsite(),0));
        phyloprocess->AddSubstitutionCounts(counts);
        for (int j=0; j<GetNbranch(); j++)  {
            int tot = 0;
            int totdouble = 0;
            for (int i=0; i<GetNsite(); i++)    {
                if (counts[j][i])   {
                    tot++;
                    if (counts[j][i] > 1)   {
                        totdouble++;
                    }
                }
            }
            double f = 0;
            if (tot)    {
                f = double(totdouble) / tot;
            }
            meancounts[j] = f;
        }
    }
    */

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "omega\t";
        os << "meantstv\tinvshape\trootgc\tmeangc\tinvconc\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << omega << '\t';
        os << kappahypermean << '\t';
        os << kappahyperinvshape << '\t';
        os << rootgamma << '\t';
        os << gammahypermean << '\t' << gammahyperinvconc << '\n';
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << omega << '\t';
        os << kappahypermean << '\t';
        os << kappahyperinvshape << '\t';
        os << *branchkappa << '\t';
        os << rootgamma << '\t';
        os << gammahypermean << '\t';
        os << gammahyperinvconc << '\t';
        os << *branchgamma << '\t';
        os << lambda << '\t';
        os << *branchlength << '\n';
    }

    void FromStream(istream &is) override {
        is >> omega;
        is >> kappahypermean >> kappahyperinvshape >> *branchkappa;
        is >> rootgamma >> gammahypermean >> gammahyperinvconc >> *branchgamma;
        is >> lambda;
        is >> *branchlength;
    }
};
