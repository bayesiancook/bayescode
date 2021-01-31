#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "BrownianTreeProcess.hpp"
#include "Chronogram.hpp"
#include "dSOmegaPathSuffStat.hpp"

class BrownianClockModel: public ProbModel {

    // (rooted) tree and data
    const Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    const CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    Chronogram* chronogram;

    double alphatau;
    double betatau;
    double tau;
    double rootrate;
    BrownianTreeProcess* lognoderate;
    BranchExpoLengthArray* branchlength;

    // Poisson suffstats for substitution histories, as a function of branch lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    PoissonSuffStat rootratesuffstat;

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
    double omegahypermean;
    double omegahyperinvshape;
    double omega;

    // a codon matrix (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrix *codonmatrix;

    // PhyloProcess

    PhyloProcess *phyloprocess;

    PathSuffStat pathsuffstat;
    OmegaPathSuffStat omegapathsuffstat;
    /*
    PathSuffStatNodeArray* pathsuffstatnodearray;
    dSOmegaPathSuffStatBranchArray* dsompathsuffstat;
    */

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    BrownianClockModel(string datafile, string treefile) {

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

    BrownianClockModel(const CodonSequenceAlignment* incodondata, const Tree* intree) {

        codondata = incodondata;
        Nsite = codondata->GetNsite();
        Ntaxa = codondata->GetNtaxa();
        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();
    }

    //! model allocation
    void Allocate() {

        chronogram = new Chronogram(*tree);
        alphatau = 1.0;
        betatau = 1.0;
        tau = 1.0;
        lognoderate = new BrownianTreeProcess(*chronogram, tau);
        rootrate = 0.1;
        branchlength = new BranchExpoLengthArray(*lognoderate, *chronogram, rootrate);
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

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        omega = 1.0;
        codonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, omega);

        /*
        pathsuffstatnodearray = new PathSuffStatNodeArray(*tree);
        dsompathsuffstat = new dSOmegaPathSuffStatBranchArray(*tree);
        */

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, codonmatrix);
        phyloprocess->Unfold();
    }

    //-------------------
    // Accessors
    // ------------------

    const Link* GetRoot() const {
        return tree->GetRoot();
    }

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return current value of omega
    double GetOmega() const { return omega; }

    //-------------------
    // Setting and updating
    // ------------------

    void TouchNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void TouchCodonMatrix() {
        codonmatrix->SetOmega(omega);
        codonmatrix->CorruptMatrix();
    }

    void TouchMatrices() {
        TouchNucMatrix();
        TouchCodonMatrix();
    }

    void NoUpdate() {}

    void Update() override {
        lognoderate->SetTau(tau);
        branchlength->SetRootRate(rootrate);
        branchlength->Update();
        TouchMatrices();
        ResampleSub(1.0);
    }

    //-------------------
    // Posterior Predictive
    // ------------------

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name) override {
        lognoderate->SetTau(tau);
        branchlength->Update();
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
        total += TauLogPrior();
        total += RootRateLogPrior();
        total += ChronoLogPrior();
        total += BrownianClockLogPrior();
        total += NucRatesLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    // Branch lengths

    double TauLogPrior() const {
        return -tau;
    }

    double RootRateLogPrior() const {
        return -rootrate;
    }

    double ChronoLogPrior() const   {
        return 0;
    }

    double BrownianClockLogPrior() const    {
        return lognoderate->GetLogProb();
    }

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

    //! log prior over omega (gamma of mean omegahypermean and shape
    //! 1/omegahyperinvshape)
    double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return Random::logGammaDensity(omega, alpha, beta);
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

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    void CollectRootRateSuffStat()  {
        rootratesuffstat.Clear();
        for (int i=0; i<Nbranch; i++)   {
            rootratesuffstat.AddSuffStat(lengthpathsuffstatarray->GetVal(i).GetCount(), lengthpathsuffstatarray->GetVal(i).GetBeta() * branchlength->GetVal(i));
        }
    }

    void GibbsResampleRootRate()    {
        double count = rootratesuffstat.GetCount();
        double beta = rootratesuffstat.GetBeta() / rootrate;
        rootrate = Random::GammaSample(1.0 + count, 1.0 + beta);
        branchlength->SetRootRate(rootrate);
    }

    // Nucleotide rates

    //! \brief const acess to nuc-pathsuffstat
    //!
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate)
    //! and equilibrium frequencies (nucstat) conditional on the current
    //! substitution mapping.
    const NucPathSuffStat &GetNucPathSuffStat() const { return nucpathsuffstat; }

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrix, pathsuffstat);
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

    // Paths

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstat.Clear();
        pathsuffstat.AddSuffStat(*phyloprocess);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // Nucleotide rates

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    double GetBranchSuffStatLogProb(const Link*from) const  {
        if (from->isRoot())   {
            return 0;
        }
        double suffstatlogprob = lengthpathsuffstatarray->GetVal(from->GetBranch()->GetIndex()).GetLogProb(branchlength->GetVal(from->GetBranch()->GetIndex()));
        return suffstatlogprob;
    }

    double GetNodeSuffStatLogProb(const Link* from) const {
        double total = GetBranchSuffStatLogProb(from);
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += GetBranchSuffStatLogProb(link->Out());
        }
        return total;
    }

    double GetNodeLogPrior(const Link* from) const  {
        return lognoderate->GetNodeLogProb(from);
    }

    double GetNodeLogProb(const Link* from) const   {
        return GetNodeLogPrior(from) + GetNodeSuffStatLogProb(from);
    }

    void NodeUpdate(const Link* from) {
        branchlength->LocalNodeUpdate(from);
    }

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

            CollectLengthSuffStat();
            MoveRootRate();
            MoveTimes();
            MoveRates();
            GibbsResampleTau();

            CollectPathSuffStat();
            MoveOmega();

            TouchMatrices();
            MoveNucRates();
        }
    }

    // Times and Rates

    void MoveRootRate() {
        CollectRootRateSuffStat();
        GibbsResampleRootRate();
    }

    void MoveTimes()    {
        RecursiveMoveTimes(1.0, GetRoot());
    }

    void RecursiveMoveTimes(double tuning, const Link* from)    {
        if (! from->isRoot())   {
            LocalMoveTime(tuning, from);
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveMoveTimes(tuning, link->Out());
        }
        if (! from->isRoot())   {
            LocalMoveTime(tuning, from);
        }
    }

    // only changes the dS
    // so, just get the 
    double LocalMoveTime(double tuning, const Link* from) {
        double logprob1 = GetNodeLogProb(from);
        double bk = chronogram->GetVal(from->GetNode()->GetIndex());
        double loghastings = chronogram->LocalProposeMove(from, tuning);
        NodeUpdate(from);
        double logprob2 = GetNodeLogProb(from);

        double deltalogprob = logprob2 - logprob1 + loghastings;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted)   {
            (*chronogram)[from->GetNode()->GetIndex()] = bk;
            NodeUpdate(from);
        }
        return ((double) accepted);
    }

    void MoveRates()    {
        RecursiveMoveRates(1.0, GetRoot());
    }

    void RecursiveMoveRates(double tuning, const Link* from)    {
        LocalMoveRate(tuning, from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveMoveRates(tuning, link->Out());
        }
        LocalMoveRate(tuning, from);
    }

    double LocalMoveRate(double tuning, const Link* from)   {
        double logprob1 = GetNodeLogProb(from);
        double bk = chronogram->GetVal(from->GetNode()->GetIndex());
        double loghastings = lognoderate->LocalProposeMove(from->GetNode()->GetIndex(), tuning);
        NodeUpdate(from);
        double logprob2 = GetNodeLogProb(from);

        double deltalogprob = logprob2 - logprob1 + loghastings;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted)   {
            (*lognoderate)[from->GetNode()->GetIndex()] = bk;
            NodeUpdate(from);
        }
        return ((double) accepted);
    }

    void GibbsResampleTau()  {
        double s2 = lognoderate->GetSumOfContrasts();
        tau = Random::Gamma(alphatau + 0.5*Nbranch, betatau + 0.5*s2);
        lognoderate->SetTau(tau);
    }

    // Nucleotide rates

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &BrownianClockModel::NucRatesLogProb,
                    &BrownianClockModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &BrownianClockModel::NucRatesLogProb,
                    &BrownianClockModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &BrownianClockModel::NucRatesLogProb,
                    &BrownianClockModel::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &BrownianClockModel::NucRatesLogProb,
                    &BrownianClockModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &BrownianClockModel::NucRatesLogProb,
                    &BrownianClockModel::TouchNucMatrix, this);

        TouchMatrices();
    }

    // Omega

    //! Gibbs resample omega (based on sufficient statistics of current
    //! substitution mapping)
    void MoveOmega() {
        omegapathsuffstat.Clear();
        omegapathsuffstat.AddSuffStat(*codonmatrix, pathsuffstat);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(),
                                    beta + omegapathsuffstat.GetBeta());
        TouchCodonMatrix();
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\t";
        os << "length\t";
        os << "tau\t";
        os << "rootrate\t";
        os << "omega\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << tau << '\t';
        os << rootrate << '\t';
        os << omega << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << omega << '\t';
        os << nucstat << '\t';
        os << nucrelrate << '\t';
        os << *chronogram << '\t';
        os << tau << '\t';
        os << rootrate << '\t';
        os << *lognoderate << '\t';
    }

    void FromStream(istream &is) override {
        is >> omega;
        is >> nucstat;
        is >> nucrelrate;
        is >> *chronogram;
        is >> tau;
        is >> rootrate;
        is >> *lognoderate;
    }
};
