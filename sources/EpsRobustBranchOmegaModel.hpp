#include "WhiteNoise.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Product.hpp"
#include "SimpleSubMatrixSelector.hpp"
#include "EpsRobustIIDGamma.hpp"
#include "Tree.hpp"

class EpsRobustBranchOmegaModel : public ProbModel {
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
    int omegamode;

    // Branch lengths

    double lambda;
    BranchIIDGamma *blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise *branchlength;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatbrancharray;

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

    double omegamean;
    double omegainvshape;
    double epsilon;
    EpsRobustBranchIIDGamma *omegabrancharray;

    // a codon matrix (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixBranchArray *codonmatrixbrancharray;
    MGOmegaCodonSubMatrix *rootcodonmatrix;

    PhyloProcess *phyloprocess;

    // suffstats
    PathSuffStatNodeArray *pathsuffstatbrancharray;
    OmegaPathSuffStatBranchArray *omegapathsuffstatbrancharray;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    EpsRobustBranchOmegaModel(string datafile, string treefile, double inepsilon)  {

        epsilon = inepsilon;

        blmode = 0;
        nucmode = 0;
        omegamode = 0;

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


    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    EpsRobustBranchOmegaModel(const CodonSequenceAlignment* incodondata, const Tree* intree, double inepsilon) {

        epsilon = inepsilon;

        blmode = 0;
        nucmode = 0;
        omegamode = 0;

	data = 0;
        codondata = incodondata;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();
    }

    //! model allocation
    void Allocate() {

        // Branch lengths

        lambda = 10.0;
        blhypermean = new BranchIIDGamma(*tree, 1.0, lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree, *blhypermean, 1.0 / blhyperinvshape);
        lengthpathsuffstatbrancharray = new PoissonSuffStatBranchArray(*tree);

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

        // omega across branches

        omegamean = 0.1;
        omegainvshape = 0.3;
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;

        omegabrancharray = new EpsRobustBranchIIDGamma(*tree, alpha, beta, epsilon);

        codonmatrixbrancharray = new MGOmegaCodonSubMatrixBranchArray(GetCodonStateSpace(), nucmatrix, omegabrancharray);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, 1.0);
        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, codonmatrixbrancharray, rootcodonmatrix);

        pathsuffstatbrancharray = new PathSuffStatNodeArray(*tree);
        omegapathsuffstatbrancharray = new OmegaPathSuffStatBranchArray(*tree);

        phyloprocess->Unfold();
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    const Tree& GetTree() const {return *tree;}

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

    void SetOmegaMode(int inomegamode)  {
        omegamode = inomegamode;
    }

    bool FixedOmega() const {
        return (omegamode == 2);
    }

    void SetBranchOmega(const BranchSelector<double>& inomegabrancharray)   {
        omegabrancharray->Copy(inomegabrancharray);
        FastUpdate();
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
        copy(innucrelrate.begin(), innucrelrate.end(), nucrelrate.begin());
        copy(innucstat.begin(), innucstat.end(), nucstat.begin());
        TouchMatrices();
    }

    //! get a copy of nucleotide rates into arrays given as arguments
    void GetNucRates(std::vector<double> &innucrelrate,
                                    std::vector<double> &innucstat) const {
        copy(nucrelrate.begin(), nucrelrate.end(), innucrelrate.begin());
        copy(nucstat.begin(), nucstat.end(), innucstat.begin());
    }

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double> &innucrelratehypercenter,
                                                   double innucrelratehyperinvconc,
                                                   const std::vector<double> &innucstathypercenter,
                                                   double innucstathyperinvconc) {
        copy(innucrelratehypercenter.begin(), innucrelratehypercenter.end(), nucrelratehypercenter.begin()); 
        nucrelratehyperinvconc = innucrelratehyperinvconc;
        copy(innucstathypercenter.begin(), innucstathypercenter.end(), nucstathypercenter.begin()); 
        nucstathyperinvconc = innucstathyperinvconc;
    }

    void SetOmegaHyperParam(double inomegamean, double inomegainvshape) {
        omegamean = inomegamean;
        omegainvshape = inomegainvshape;
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegabrancharray->SetShape(alpha);
        omegabrancharray->SetScale(beta);
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

    void TouchCodonMatrices()	{
        codonmatrixbrancharray->UpdateCodonMatrices();
        rootcodonmatrix->CorruptMatrix();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices() {
        TouchNucMatrix();
        TouchCodonMatrices();
    }

    void FastUpdate() {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegabrancharray->SetShape(alpha);
        omegabrancharray->SetScale(beta);
        TouchMatrices();
    }

    void Update() override {
        FastUpdate();
        ResampleSub(1.0);
    }

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name) override {
        FastUpdate();
        phyloprocess->PostPredSample(name);
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
        if (!FixedBranchLengths()) {
            total += BranchLengthsLogPrior();
        }
        if (!FixedNucRates()) {
            total += NucRatesLogPrior();
        }
        if (!FixedOmega())    {
            total += OmegaLogPrior();
            total += OmegaHyperLogPrior();
        }
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

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

    //! log prior over omega
    double OmegaLogPrior() const { return omegabrancharray->GetLogProb(); }

    double OmegaHyperLogPrior() const   {
        double total = 0;
        total -= omegamean;
        total -= omegainvshape;
        return total;
    }

    void UpdateOmega()  {
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegabrancharray->SetShape(alpha);
        omegabrancharray->SetScale(beta);
    }

    double OmegaHyperLogProb() const    {
        return OmegaHyperLogPrior() + OmegaLogPrior();
    }

    const OmegaPathSuffStatBranchArray& GetOmegaPathSuffStatBranchArray() const {
        return *omegapathsuffstatbrancharray;
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
        return lengthpathsuffstatbrancharray;
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatbrancharray->Clear();
        lengthpathsuffstatbrancharray->AddLengthPathSuffStat(*phyloprocess);
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

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrices. Both pathsuffstats and
    //! codonmatrices are assumed to be updated.
    double PathSuffStatLogProb() const {
        return pathsuffstatbrancharray->GetLogProb(*codonmatrixbrancharray, *rootcodonmatrix);
    }

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

            if (!FixedNucRates()) {
                MoveNucRates();
            }

            if (!FixedOmega())  {
                MoveOmega();
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
        branchlength->GibbsResample(*lengthpathsuffstatbrancharray);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &EpsRobustBranchOmegaModel::LambdaHyperLogProb, &EpsRobustBranchOmegaModel::NoUpdate,
                    this);
        ScalingMove(lambda, 0.3, 10, &EpsRobustBranchOmegaModel::LambdaHyperLogProb, &EpsRobustBranchOmegaModel::NoUpdate,
                    this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    // Nucleotide rates

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &EpsRobustBranchOmegaModel::NucRatesLogProb,
                    &EpsRobustBranchOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &EpsRobustBranchOmegaModel::NucRatesLogProb,
                    &EpsRobustBranchOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &EpsRobustBranchOmegaModel::NucRatesLogProb,
                    &EpsRobustBranchOmegaModel::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &EpsRobustBranchOmegaModel::NucRatesLogProb,
                    &EpsRobustBranchOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &EpsRobustBranchOmegaModel::NucRatesLogProb,
                    &EpsRobustBranchOmegaModel::TouchNucMatrix, this);

        TouchMatrices();
    }

    // Omega

    void MoveOmega()    {
        CollectOmegaSuffStat();
        MoveBranchOmega();
        MoveOmegaHyperParameters();
        TouchMatrices();
    }

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstatbrancharray->Clear();
        pathsuffstatbrancharray->AddSuffStat(*phyloprocess);
    }

    void CollectOmegaSuffStat() {
        omegapathsuffstatbrancharray->Clear();
        omegapathsuffstatbrancharray->AddSuffStat(*codonmatrixbrancharray, *rootcodonmatrix, *pathsuffstatbrancharray);
    }

    //! Gibbs resample omega (based on sufficient statistics of current
    //! substitution mapping)
    void MoveBranchOmega() {
        omegabrancharray->Move(1.0, 10, *omegapathsuffstatbrancharray);
        omegabrancharray->Move(0.1, 10, *omegapathsuffstatbrancharray);
        codonmatrixbrancharray->UpdateCodonMatrices();
    }

    void MoveOmegaHyperParameters() {
        ScalingMove(omegamean, 1.0, 10, &EpsRobustBranchOmegaModel::OmegaHyperLogProb,
                    &EpsRobustBranchOmegaModel::UpdateOmega, this);
        ScalingMove(omegamean, 0.3, 10, &EpsRobustBranchOmegaModel::OmegaHyperLogProb,
                    &EpsRobustBranchOmegaModel::UpdateOmega, this);
        ScalingMove(omegainvshape, 1.0, 10, &EpsRobustBranchOmegaModel::OmegaHyperLogProb,
                    &EpsRobustBranchOmegaModel::UpdateOmega, this);
        ScalingMove(omegainvshape, 0.3, 10, &EpsRobustBranchOmegaModel::OmegaHyperLogProb,
                    &EpsRobustBranchOmegaModel::UpdateOmega, this);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    double GetMeanOmega() const {
        return omegabrancharray->GetMean();
    }

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixbrancharray, *rootcodonmatrix, *pathsuffstatbrancharray);
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
        os << omegabrancharray->GetMean() << '\t' << omegabrancharray->GetVar() << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void TraceEpsRobustBranchOmega(ostream& os) const {
        for (int i=0; i<Nbranch; i++) {
            os << omegabrancharray->GetVal(i) << '\t';
        }
        os << '\n';
        os.flush();
    }

    const BranchSelector<double>& GetOmegaTree() const  {
        return *omegabrancharray;
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << lambda << '\t';
        os << *branchlength << '\t';
        os << omegamean << '\t' << omegainvshape << '\t';
        os << *omegabrancharray << '\t';
        os << nucrelrate << '\t';
        os << nucstat << '\n';
    }

    void FromStream(istream &is) override {
        is >> lambda;
        is >> *branchlength;
        is >> omegamean >> omegainvshape;
        is >> *omegabrancharray;
        is >> nucrelrate;
        is >> nucstat;
    }
};
