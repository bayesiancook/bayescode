#include "ContinuousData.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "MultivariateBrownianTreeProcess.hpp"
#include "Chronogram.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include "InverseWishart.hpp"

class FastCoevolModel: public ProbModel {

    const Tree *tree;
    const TaxonSet *taxonset;
    const ContinuousData* contdata;
    string dsomsuffstatfile;

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

    dSOmegaPathSuffStatBranchArray* dsompathsuffstatarray;
    MultivariateNormalSuffStat* browniansuffstat;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    FastCoevolModel(string contdatafile, string treefile, string rootfile, string indsomsuffstatfile)   {

        dsomsuffstatfile = indsomsuffstatfile;

	contdata = new FileContinuousData(contdatafile);
	Ncont = contdata->GetNsite();
        Ntaxa = contdata->GetNtaxa();
        taxonset = contdata->GetTaxonSet();

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

        dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);
        browniansuffstat = new MultivariateNormalSuffStat(process->GetDim());
        ifstream is(dsomsuffstatfile.c_str());
        is >> *dsompathsuffstatarray;

        cerr << "allocate ok\n";
    }

    void SetCoevol(const NodeSelector<double>& inchronogram, const NodeSelector<vector<double> >& inprocess)   {
        chronogram->Copy(inchronogram);
        process->Copy(inprocess);
        Update();
    }

    //-------------------
    // Accessors
    // ------------------

    const dSOmegaPathSuffStatBranchArray& GetdSOmegaPathSuffStatBranchArray() const {
        return *dsompathsuffstatarray;
    }

    const Tree& GetTree() const {
        return *tree;
    }

    const Link* GetRoot() const {
        return tree->GetRoot();
    }

    void NoUpdate() {}

    void Update() override {
        branchlength->Update();
        branchomega->Update();
    }

    void PostPred(string name) override {
        Update();
        cerr << "in post pred\n";
        exit(1);
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
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { 
        return dSOmPathSuffStatLogProb();
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

    //-------------------
    //  Log probs for MH moves
    //-------------------

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
        MoveParameters(30);
        return 1.0;
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveCoevol();
        }
    }

    void MoveCoevol()   {
        for (int rep=0; rep<5; rep++)   {
            MoveTimes();
            MoveBrownianProcess();
            MoveSigma();
            MoveKappa();
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
        os << '\n';
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << *chronogram << '\t';
        os << kappa << '\t';
        os << *sigma << '\t';
        os << *process << '\n';
    }

    void FromStream(istream &is) override {
        is >> *chronogram;
        is >> kappa;
        is >> *sigma;
        is >> *process;
    }
};
