#include "CodonSequenceAlignment.hpp"
#include "ContinuousData.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "MultivariateBrownianTreeProcess.hpp"
#include "Chronogram.hpp"
#include "PoissonSuffStat.hpp"
#include "CodonSuffStat.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include "InverseWishart.hpp"

class CoevolModel: public ProbModel {

    const Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    const CodonSequenceAlignment *codondata;
    const ContinuousData* contdata;

    int Nsite;
    int Ntaxa;
    int Nbranch;
    int Ncont;
    int L;
    int dSindex;
    int omindex;


    Chronogram* chronogram;

    int df;
    vector<double> kappa;
    InverseWishart* sigma;

    vector<double> rootmean;
    vector<double> rootvar;

    MultivariateBrownianTreeProcess* process;
    MVBranchExpoLengthArray* branchlength;
    MVBranchExpoMeanArray* branchomega;

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

    CoevolModel(string datafile, string contdatafile, string treefile, string rootfile) {

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}
        L = 2;
        dSindex = 0;
        omindex = 0;

        // get tree from file (newick format)
        Tree* tmptree = new Tree(treefile);
        // check whether tree and data fits together
        tmptree->RegisterWith(taxonset);
        tmptree->SetIndices();
        tree = tmptree;

        Nbranch = tree->GetNbranch();

        ifstream is(rootfile.c_str());
        int dim;
        is >> dim;
        if (dim != Ncont + L)   {
            cerr << "error in root file: non matching dimension\n";
            exit(1);
        }
        rootmean.assign(dim,0);
        rootvar.assign(dim,0);
        for (int i=0; i<dim; i++)   {
            is >> rootmean[i] >> rootvar[i];
        }
    }

    CoevolModel(const CodonSequenceAlignment* incodondata, const ContinuousData* incontdata, const Tree* intree) {

        codondata = incodondata;
        Nsite = codondata->GetNsite();
        Ntaxa = codondata->GetNtaxa();
        taxonset = codondata->GetTaxonSet();

        contdata = incontdata;
        Ncont = contdata->GetNsite();

        tree = intree;
        Nbranch = tree->GetNbranch();
    }

    //! model allocation
    void Allocate() {

        cerr << "allocate\n";

        chronogram = new Chronogram(*tree);

        kappa.assign(Ncont+L, 1.0);
        df = 0;
        sigma = new InverseWishart(kappa, df);

        process = new MultivariateBrownianTreeProcess(*chronogram, *sigma, rootmean, rootvar);
        for (int i=0; i<Ncont; i++)	{
            process->SetAndClamp(*contdata, L+i, i);
        }

        branchlength = new MVBranchExpoLengthArray(*process, *chronogram, dSindex);
        branchomega = new MVBranchExpoMeanArray(*process, omindex);

        cerr << "total length : " << branchlength->GetTotalLength() << '\n';
        cerr << "mean omega   : " << branchomega->GetMean() << '\n';

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

        cerr << "allocate ok\n";
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
        branchlength->Update();
        branchomega->Update();
        TouchMatrices();
        ResampleSub(1.0);
    }

    void PostPred(string name) override {
        branchlength->Update();
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
        total += KappaLogPrior();
        total += SigmaLogPrior();
        total += BrownianProcessLogPrior();
        total += NucRatesLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    // Branch lengths

    double ChronoLogPrior() const   {
        return 0;
    }

    double KappaLogPrior() const    {
        return 0;
    }

    double SigmaLogPrior() const    {
        return sigma->GetLogProb();
    }

    double BrownianProcessLogPrior() const    {
        return process->GetLogProb();
    }

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

    const NucPathSuffStat &GetNucPathSuffStat() const { return nucpathsuffstat; }

    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixarray, *rootcodonmatrix, *pathsuffstatarray);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    double BranchSuffStatLogProb(const Link*from) const  {
        if (from->isRoot())   {
            return 0;
        }
        int index = from->GetBranch()->GetIndex();
        double suffstatlogprob = dsompathsuffstatarray->GetVal(index).GetLogProb(branchlength->GetVal(index), branchomega->GetVal(index));
        return suffstatlogprob;
    }

    double NodeSuffStatLogProb(const Link* from) const {
        double total = BranchSuffStatLogProb(from);
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += BranchSuffStatLogProb(link->Out());
        }
        return total;
    }

    /*
    double dSOmegaSuffStatLogProb() const    {
        return RecursivedSOmegaSuffStatLogProb(GetRoot());
    }

    double RecursivedSOmegaSuffStatLogProb(const Link* from) const {
        double total = 0;
        if (! from->isRoot())   {
            total += BranchSuffStatLogProb(from);
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += RecursivedSOmegaSuffStatLogProb(link->Out());
        }
        return total;
    }

    double dSOmegaLogProb() const    {
        return ChronoLogPrior() + BrownianProcessLogPrior() + RootValLogPrior() + dSOmegaSuffStatLogProb();
    }
    */

    double NodeLogPrior(const Link* from) const  {
        return process->GetNodeLogProb(from);
    }

    double NodeLogProb(const Link* from) const   {
        return NodeLogPrior(from) + NodeSuffStatLogProb(from);
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

            CollectPathSuffStat();

            CollectdSOmegaPathSuffStat();
            MoveTimes();
            MoveBrownianProcess();
            // MoveSigma();

            CollectNucPathSuffStat();
            TouchMatrices();
            MoveNucRates();
            TouchMatrices();
        }
    }

    // Times and Rates

    void MoveTimes()    {
        chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
    }

    void MoveBrownianProcess()  {
        for (int i=0; i<L+Ncont; i++)   {
            process->SingleNodeMove(i, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
        }
    }

    void MoveNucRates() {

        ProfileMove(nucrelrate, 0.1, 1, 3, &CoevolModel::NucRatesLogProb,
                    &CoevolModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &CoevolModel::NucRatesLogProb,
                    &CoevolModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &CoevolModel::NucRatesLogProb,
                    &CoevolModel::TouchNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &CoevolModel::NucRatesLogProb,
                    &CoevolModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &CoevolModel::NucRatesLogProb,
                    &CoevolModel::TouchNucMatrix, this);

        // TouchMatrices();
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\t";
        os << "length\t";
        os << "meanomega\t";
        for (int i=0; i<process->GetDim(); i++) {
            for (int j=i+1; j<process->GetDim(); j++)   {
                os << "s_" << i << "_" << j << '\t';
            }
        }
        for (int i=0; i<process->GetDim(); i++) {
            os << "s_" << i << "_" << i << '\t';
        }
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << branchomega->GetMean() << '\t';
        for (int i=0; i<process->GetDim(); i++) {
            for (int j=i+1; j<process->GetDim(); j++)   {
                os << (*sigma)(i,j) << '\t';
            }
        }
        for (int i=0; i<process->GetDim(); i++) {
            os << (*sigma)(i,i) << '\t';
        }
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << nucstat << '\t';
        os << nucrelrate << '\t';
        os << *chronogram << '\t';
        os << *sigma << '\t';
        os << *process << '\t';
    }

    void FromStream(istream &is) override {
        is >> nucstat;
        is >> nucrelrate;
        is >> *chronogram;
        is >> *sigma;
        is >> *process;
    }
};
