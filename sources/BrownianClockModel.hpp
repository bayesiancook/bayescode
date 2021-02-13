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
#include "PoissonSuffStat.hpp"
#include "CodonSuffStat.hpp"
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

    double tauds;
    double rootds;

    double tauom;
    double rootomega;

    BrownianTreeProcess* logds;
    BrownianTreeProcess* logomega;
    BranchExpoLengthArray* branchlength;
    BranchExpoMeanArray* branchomega;

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

    // a codon matrix (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixBranchArray* codonmatrixarray;
    MGOmegaCodonSubMatrix* rootcodonmatrix;

    // PhyloProcess

    PhyloProcess *phyloprocess;

    PathSuffStatNodeArray* pathsuffstatarray;
    dSOmegaPathSuffStatBranchArray* dsompathsuffstatarray;
    PoissonSuffStat rootdssuffstat;
    OmegaPathSuffStat rootomegasuffstat;


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

        tauds = 1.0;
        tauom = 1.0;

        logds = new BrownianTreeProcess(*chronogram, tauds);
        rootds = 0.1;
        branchlength = new BranchExpoLengthArray(*logds, *chronogram, rootds);

        logomega = new BrownianTreeProcess(*chronogram, tauom);
        rootomega = 0.1;
        branchomega = new BranchExpoMeanArray(*logomega, rootomega);

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

        codonmatrixarray = new MGOmegaCodonSubMatrixBranchArray(GetCodonStateSpace(), nucmatrix, branchomega);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(GetCodonStateSpace(), nucmatrix, 1.0);

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, codonmatrixarray, rootcodonmatrix);
        phyloprocess->Unfold();

        pathsuffstatarray = new PathSuffStatNodeArray(*tree);
        dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);
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

    //-------------------
    // Setting and updating
    // ------------------

    void TouchNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void TouchCodonMatrices() {
        codonmatrixarray->UpdateCodonMatrices();
        rootcodonmatrix->CorruptMatrix();
    }

    void TouchMatrices() {
        TouchNucMatrix();
        TouchCodonMatrices();
    }

    void NoUpdate() {}

    void Update() override {
        logds->SetTau(tauds);
        branchlength->SetRootShift(rootds);
        branchlength->Update();

        logomega->SetTau(tauom);
        branchomega->SetRootShift(rootomega);
        branchomega->Update();

        TouchMatrices();
        ResampleSub(1.0);
    }

    //-------------------
    // Posterior Predictive
    // ------------------

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name) override {
        logds->SetTau(tauds);
        branchlength->SetRootShift(rootds);
        branchlength->Update();

        logomega->SetTau(tauom);
        branchomega->SetRootShift(rootomega);
        branchomega->Update();

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
        total += ChronoLogPrior();

        total += TaudSLogPrior();
        total += BrowniandSLogPrior();
        total += RootdSLogPrior();

        total += TauOmegaLogPrior();
        total += BrownianOmegaLogPrior();
        total += RootOmegaLogPrior();

        total += NucRatesLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    // Branch lengths

    double TaudSLogPrior() const {
        return -tauds;
    }

    double TauOmegaLogPrior() const {
        return -tauom;
    }

    double RootdSLogPrior() const {
        return -rootds;
    }

    double RootOmegaLogPrior() const    {
        return -rootomega;
    }

    double ChronoLogPrior() const   {
        return 0;
    }

    double BrowniandSLogPrior() const    {
        return logds->GetLogProb();
    }

    double BrownianOmegaLogPrior() const    {
        return logomega->GetLogProb();
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

    void CollectPathSuffStat() {
        pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    void CollectdSOmegaPathSuffStat() {
        dsompathsuffstatarray->Clear();
        dsompathsuffstatarray->AddSuffStat(*codonmatrixarray, *pathsuffstatarray);
    }

    void CollectRootdSSuffStat()  {
        rootdssuffstat.Clear();
        for (int i=0; i<Nbranch; i++)   {
            const dSOmegaPathSuffStat& ss = dsompathsuffstatarray->GetVal(i);
            double omega = branchomega->GetVal(i);
            double length = branchlength->GetVal(i);
            if (std::isnan(length)) {
                cerr << "in collect root dsss : length is nan\n";
                exit(1);
            }
            if (std::isnan(omega)) {
                cerr << "in collect root dsss : omega is nan\n";
                exit(1);
            }
            rootdssuffstat.AddSuffStat(ss.GetCount(), ss.GetBeta(omega) * length);
        }
    }

    void GibbsResampleRootdS()    {
        double count = rootdssuffstat.GetCount();
        double beta = rootdssuffstat.GetBeta() / rootds;
        cerr << "gibbs resample before \n";
        cerr << rootdssuffstat.GetBeta() << '\t' << rootds << '\n';
        rootds = Random::GammaSample(1.0 + count, 1.0 + beta);
        if (std::isnan(rootds)) {
            cerr << "rootds is nan\n";
            cerr << count << '\t' << beta << '\n';
            exit(1);
        }
        branchlength->SetRootShift(rootds);
        branchlength->Update();
    }

    void CollectRootOmegaSuffStat()  {
        rootomegasuffstat.Clear();
        for (int i=0; i<Nbranch; i++)   {
            const dSOmegaPathSuffStat& ss = dsompathsuffstatarray->GetVal(i);
            double omega = branchomega->GetVal(i);
            double length = branchlength->GetVal(i);
            if (std::isnan(length)) {
                cerr << "in collect root domss : length is nan\n";
                exit(1);
            }
            if (std::isnan(omega)) {
                cerr << "in collect root domss : omega is nan\n";
                exit(1);
            }
            rootomegasuffstat.PoissonSuffStat::AddSuffStat(ss.GetNonSynCount(), ss.GetNonSynBeta()*length*omega);
        }
    }

    void GibbsResampleRootOmega()    {
        double count = rootomegasuffstat.GetCount();
        double beta = rootomegasuffstat.GetBeta() / rootomega;
        rootomega = Random::GammaSample(1.0 + count, 1.0 + beta);
        if (std::isnan(rootomega)) {
            cerr << "rootomega is nan\n";
            cerr << count << '\t' << beta << '\n';
            exit(1);
        }
        branchomega->SetRootShift(rootomega);
        branchomega->Update();
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
        nucpathsuffstat.AddSuffStat(*codonmatrixarray, *rootcodonmatrix, *pathsuffstatarray);
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
        int index = from->GetBranch()->GetIndex();
        double suffstatlogprob = dsompathsuffstatarray->GetVal(index).GetLogProb(branchlength->GetVal(index), branchomega->GetVal(index));
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
        return logds->GetNodeLogProb(from) + logomega->GetNodeLogProb(from);
    }

    double GetNodeLogProb(const Link* from) const   {
        return GetNodeLogPrior(from) + GetNodeSuffStatLogProb(from);
    }

    void NodeUpdate(const Link* from) {
        branchlength->LocalNodeUpdate(from);
        branchomega->LocalNodeUpdate(from);
    }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() override {
        cerr << "resample sub\n";
        ResampleSub(1.0);
        cerr << "move param\n";
        MoveParameters(30);
        cerr << "move ok\n";
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

            cerr << "collect\n";
            CollectPathSuffStat();
            CollectdSOmegaPathSuffStat();

            cerr << "times\n";
            MoveTimes();

            cerr << "collect\n";
            CollectdSOmegaPathSuffStat();
            MoveRootdS();
            cerr << "collect\n";
            CollectdSOmegaPathSuffStat();
            MovedS();
            GibbsResampleTaudS();

            cerr << "collect\n";
            CollectdSOmegaPathSuffStat();
            MoveRootOmega();
            cerr << "collect\n";
            CollectdSOmegaPathSuffStat();
            MoveOmega();
            GibbsResampleTauOmega();

            cerr << "touch\n";
            TouchMatrices();
            cerr << "nuc\n";
            MoveNucRates();
            cerr << "nuc ok\n";

            cerr << rootds << '\t' << rootomega << '\n';
        }
    }

    // Times and Rates

    void MoveRootdS() {
        CollectRootdSSuffStat();
        GibbsResampleRootdS();
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

    void MovedS()    {
        RecursiveMovedS(1.0, GetRoot());
    }

    void RecursiveMovedS(double tuning, const Link* from)    {
        LocalMovedS(tuning, from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveMovedS(tuning, link->Out());
        }
        LocalMovedS(tuning, from);
    }

    double LocalMovedS(double tuning, const Link* from)   {
        double logprob1 = GetNodeLogProb(from);
        double bk = chronogram->GetVal(from->GetNode()->GetIndex());
        double loghastings = logds->LocalProposeMove(from->GetNode()->GetIndex(), tuning);
        NodeUpdate(from);
        double logprob2 = GetNodeLogProb(from);

        double deltalogprob = logprob2 - logprob1 + loghastings;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted)   {
            (*logds)[from->GetNode()->GetIndex()] = bk;
            NodeUpdate(from);
        }
        return ((double) accepted);
    }

    void GibbsResampleTaudS()  {
        double s2 = logds->GetSumOfContrasts();
        tauds = Random::Gamma(1.0 + 0.5*Nbranch, 1.0 + 0.5*s2);
        logds->SetTau(tauds);
    }

    void MoveOmega()    {
        RecursiveMoveOmega(1.0, GetRoot());
    }

    void MoveRootOmega() {
        CollectRootOmegaSuffStat();
        GibbsResampleRootOmega();
    }

    void RecursiveMoveOmega(double tuning, const Link* from)    {
        LocalMoveOmega(tuning, from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveMoveOmega(tuning, link->Out());
        }
        LocalMoveOmega(tuning, from);
    }

    double LocalMoveOmega(double tuning, const Link* from)   {
        double logprob1 = GetNodeLogProb(from);
        double bk = chronogram->GetVal(from->GetNode()->GetIndex());
        double loghastings = logomega->LocalProposeMove(from->GetNode()->GetIndex(), tuning);
        NodeUpdate(from);
        double logprob2 = GetNodeLogProb(from);

        double deltalogprob = logprob2 - logprob1 + loghastings;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted)   {
            (*logomega)[from->GetNode()->GetIndex()] = bk;
            NodeUpdate(from);
        }
        return ((double) accepted);
    }

    void GibbsResampleTauOmega()  {
        double s2 = logomega->GetSumOfContrasts();
        tauom = Random::Gamma(1.0 + 0.5*Nbranch, 1.0 + 0.5*s2);
        logomega->SetTau(tauom);
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

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\t";
        os << "length\t";
        os << "meanomega\t";
        os << "tauds\t";
        os << "rootds\t";
        os << "tauom\t";
        os << "rootom\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << branchomega->GetMean() << '\t';
        os << tauds << '\t';
        os << rootds << '\t';
        os << tauom << '\t';
        os << rootomega << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << nucstat << '\t';
        os << nucrelrate << '\t';
        os << *chronogram << '\t';
        os << tauds << '\t';
        os << rootds << '\t';
        os << *logds << '\t';
        os << tauom << '\t';
        os << rootomega << '\t';
        os << *logomega << '\n';
    }

    void FromStream(istream &is) override {
        is >> nucstat;
        is >> nucrelrate;
        is >> *chronogram;
        is >> tauds;
        is >> rootds;
        is >> *logds;
        is >> tauom;
        is >> rootomega;
        is >> *logomega;
    }
};
