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
#include "RelativePathSuffStat.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include "InverseWishart.hpp"

class CoevolModel: public ProbModel {

    const Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    const CodonSequenceAlignment *codondata;
    const ContinuousData* contdata;

    string suffstatfile;
    string dsomsuffstatfile;
    int mappingapprox;

    int nucmode;
    int coevolmode;

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
    RelativePathSuffStatNodeArray* relpathsuffstatarray;
    dSOmegaPathSuffStatBranchArray* dsompathsuffstatarray;
    MultivariateNormalSuffStat* browniansuffstat;

    int relative;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    CoevolModel(string datafile, string contdatafile, string treefile, string rootfile, string insuffstatfile, string indsomsuffstatfile, GeneticCodeType codetype) {

        relative = 1;
        coevolmode = 0;
        suffstatfile = insuffstatfile;
        dsomsuffstatfile = indsomsuffstatfile;
        if (dsomsuffstatfile != "None")  {
            mappingapprox = 2;
        } else if (suffstatfile != "None")  {
            mappingapprox = 1;
        }
        else    {
            mappingapprox = 0;
        }

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true, codetype);

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
        omindex = 1;

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
            cerr << dim << '\t' << Ncont << '\t' << L << '\n';
            exit(1);
        }
        rootmean.assign(dim,0);
        rootvar.assign(dim,0);
        for (int i=0; i<dim; i++)   {
            is >> rootmean[i] >> rootvar[i];
        }
    }

    CoevolModel(const CodonSequenceAlignment* incodondata, const ContinuousData* incontdata, const Tree* intree, const vector<double>& inrootmean, const vector<double>& inrootvar) {

        coevolmode = 2;

        codondata = incodondata;
        Nsite = codondata->GetNsite();
        Ntaxa = codondata->GetNtaxa();
        taxonset = codondata->GetTaxonSet();

        contdata = incontdata;

		if (contdata)   {
			Ncont = contdata->GetNsite();
		}
		else	{
			Ncont = 0;
		}

        L = 2;
        dSindex = 0;
        omindex = 1;

        tree = intree;
        Nbranch = tree->GetNbranch();

        rootmean.assign(Ncont+L,0);
        rootvar.assign(Ncont+L,0);
        for (int i=0; i<Ncont+L; i++)   {
            rootmean[i] = inrootmean[i];
            rootvar[i] = inrootvar[i];
        }
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

        if (mappingapprox == 2) {
            dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);
            browniansuffstat = new MultivariateNormalSuffStat(process->GetDim());
            ifstream is(dsomsuffstatfile.c_str());
            is >> *dsompathsuffstatarray;
        }

        else    {

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
            relpathsuffstatarray = 0;
            if (relative)   {
                relpathsuffstatarray = new RelativePathSuffStatNodeArray(*tree, codondata->GetNstate());
            }
            dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);
            browniansuffstat = new MultivariateNormalSuffStat(process->GetDim());
            cerr << "allocate ok\n";
        }
    }

    void SetCoevolMode(int inmode)  {
        coevolmode = inmode;
    }

    void SetNucMode(int innucmode)  {
        nucmode = innucmode;
    }

    void SetCoevol(const NodeSelector<double>& inchronogram, const NodeSelector<vector<double> >& inprocess)   {
        chronogram->Copy(inchronogram);
        process->Copy(inprocess);
        FastUpdate();
    }

    //-------------------
    // Accessors
    // ------------------

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

    const dSOmegaPathSuffStatBranchArray& GetdSOmegaPathSuffStatBranchArray() const {
        return *dsompathsuffstatarray;
    }

    const Tree& GetTree() const {
        return *tree;
    }

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

    void FastUpdate() {
        branchlength->Update();
        branchomega->Update();
        if (mappingapprox < 2)  {
            TouchMatrices();
        }
    }

    void Update() override {
        FastUpdate();
        if (!mappingapprox) {
            ResampleSub(1.0);
        }
    }

    void PostPred(string name) override {
        FastUpdate();
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
        if (coevolmode < 2) {
            total += ChronoLogPrior();
            total += KappaLogPrior();
            total += SigmaLogPrior();
            total += BrownianProcessLogPrior();
        }
        if (mappingapprox < 2)  {
            if (!FixedNucRates()) {
                total += NucRatesLogPrior();
            }
        }
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { 
        if (mappingapprox == 2) {
            return dSOmPathSuffStatLogProb();
        }
        return phyloprocess->GetLogLikelihood(); 
    }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    double ChronoLogPrior() const   {
        return 0;
    }

    double KappaLogPrior() const    {
        double total = 0;
        for (int i=0; i<sigma->GetDim(); i++)   {
            total -= kappa[i] / 10;
        }
        return total;
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
        if (relative)   {
            relpathsuffstatarray->Clear();
            relpathsuffstatarray->AddSuffStat(*pathsuffstatarray, *branchlength);
        }
    }

    void CollectdSOmegaPathSuffStat() {
        dsompathsuffstatarray->Clear();
        if (relative)   {
            dsompathsuffstatarray->AddSuffStat(*codonmatrixarray, *relpathsuffstatarray, *branchomega);
        }
        else    {
            dsompathsuffstatarray->AddSuffStat(*codonmatrixarray, *pathsuffstatarray, *branchlength, *branchomega);
        }
    }

    const NucPathSuffStat &GetNucPathSuffStat() const { return nucpathsuffstat; }

    void CollectNucPathSuffStat() {
        TouchMatrices();
        nucpathsuffstat.Clear();
        if (relative)   {
            nucpathsuffstat.AddSuffStat(*codonmatrixarray, *rootcodonmatrix, *relpathsuffstatarray, *branchlength);
        }
        else    {
            nucpathsuffstat.AddSuffStat(*codonmatrixarray, *rootcodonmatrix, *pathsuffstatarray);
        }
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    double dSOmPathSuffStatLogProb() const {
        double total = 0;
        for (int index=0; index<tree->GetNbranch(); index++)    {
            total += dsompathsuffstatarray->GetVal(index).GetLogProb(branchlength->GetVal(index), branchomega->GetVal(index));
        }
        return total;
    }

    double BranchSuffStatLogProb(const Link*from) const  {
        if (from->isRoot())   {
            return 0;
        }
        int index = from->GetBranch()->GetIndex();
        return dsompathsuffstatarray->GetVal(index).GetLogProb(branchlength->GetVal(index), branchomega->GetVal(index));
    }

    double NodeSuffStatLogProb(const Link* from) const {
        double total = BranchSuffStatLogProb(from);
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += BranchSuffStatLogProb(link->Out());
        }
        return total;
    }

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

    double KappaSuffStatLogProb() const {
        return sigma->GetLogProb();
    }

    double KappaLogProb() const {
        return KappaLogPrior() + KappaSuffStatLogProb();
    }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() override {
        if (! mappingapprox)    {
            ResampleSub(1.0);
        }
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        TouchMatrices();
        phyloprocess->Move(frac);
        if (relative)   {
            CollectPathSuffStat();
        }
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            if (! mappingapprox)    {
                if (! relative) {
                    CollectPathSuffStat();
                }
            }
            if (mappingapprox < 2)  {
                if (!FixedNucRates()) {
                    MoveNuc();
                }
            }
            if (coevolmode < 2) {
                MoveCoevol();
            }
        }
    }

    void MoveCoevol()   {
        if (mappingapprox < 2)  {
            CollectdSOmegaPathSuffStat();
        }
        for (int rep=0; rep<5; rep++)   {
            MoveTimes();
            MoveBrownianProcess();
            MoveSigma();
            MoveKappa();
        }
        if (mappingapprox < 2)  {
            TouchMatrices();
        }
    }

    // Times and Rates

    void MoveTimes()    {
        chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
    }

    void MoveBrownianProcess()  {
        for (int i=0; i<L+Ncont; i++)   {
            process->SingleNodeMove(i, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
            process->SingleNodeMove(i, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
        }
    }

    void MoveSigma()    {
        browniansuffstat->Clear();
        process->AddSuffStat(*browniansuffstat);
        sigma->GibbsResample(*browniansuffstat);
    }

    void MoveKappa()    {
        int nrep = 10;
        double tuning = 1.0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<sigma->GetDim(); i++)   {
                double logprob1 = KappaLogProb();
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                kappa[i] *= e;
                double logprob2 = KappaLogProb();
                double loghastings = m;
                double deltalogprob = logprob2 - logprob1 + loghastings;
                int accept = (log(Random::Uniform()) < deltalogprob);
                if (! accept)   {
                    kappa[i] /= e;
                }
            }
        }
    }

    void MoveNuc() {
        CollectNucPathSuffStat();
        MoveNucRates();
        TouchMatrices();
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
        for (int i=0; i<process->GetDim(); i++) {
            os << "k_" << i << '\t';
        }
        if (mappingapprox < 2)  {
            os << "statent\t";
            os << "rrent\t";
        }
        os << '\n';
    }

    double GetMeanOmega() const	{
        return branchomega->GetMean();
    }

    const BranchSelector<double>& GetOmegaTree() const  {
        return *branchomega;
    }

    const MultivariateBrownianTreeProcess* GetProcess() const {
        return process;
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
        for (int i=0; i<process->GetDim(); i++) {
            os << kappa[i] << '\t';
        }
        if (mappingapprox < 2)  {
            os << Random::GetEntropy(nucstat) << '\t';
            os << Random::GetEntropy(nucrelrate) << '\t';
        }
        os << '\n';
    }

    void TracedSOmegaPathSuffStat(ostream& os) const    {
        if (mappingapprox < 2)  {
            os << *dsompathsuffstatarray;
            os << '\n';
            os.flush();
        }
    }

    void TraceRelativePathSuffStat(ostream& os) const   {
        if (! mappingapprox)  {
            os << *relpathsuffstatarray;
            os << '\n';
            os.flush();
        }
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        if (mappingapprox < 2)  {
            os << nucstat << '\t';
            os << nucrelrate << '\t';
        }
        os << *chronogram << '\t';
        os << kappa << '\t';
        os << *sigma << '\t';
        os << *process << '\n';
    }

    void FromStream(istream &is) override {
        if (mappingapprox < 2)  {
            is >> nucstat;
            is >> nucrelrate;
        }
        is >> *chronogram;
        is >> kappa;
        is >> *sigma;
        is >> *process;
    }
};
