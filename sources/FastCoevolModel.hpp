#include "ContinuousData.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "MultivariateBrownianTreeProcess.hpp"
#include "Chronogram.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include "InverseWishart.hpp"
#include "ChronoWhiteNoise.hpp"

class FastCoevolModel: public ProbModel {

    int wndsmode;
    int wnommode;

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

    double nuds;
    ChronoGammaWhiteNoise* wnds;
    PoissonSuffStatBranchArray* wndssuffstatbrancharray;

    double nuom;
    ChronoGammaWhiteNoise* wnom;
    PoissonSuffStatBranchArray* wnomsuffstatbrancharray;

    dSOmegaPathSuffStatBranchArray* dsompathsuffstatarray;
    MultivariateNormalSuffStat* browniansuffstat;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    FastCoevolModel(string contdatafile, string treefile, string rootfile, string indsomsuffstatfile, int inwndsmode, int inwnommode)   {

        wndsmode = inwndsmode;
        wnommode = inwnommode;

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
        // chronogram->SetAgesFromTree();

        kappa.assign(Ncont+L, 1.0);
        df = 0;
        sigma = new InverseWishart(kappa, df);
        /*
        for (int i=0; i<Ncont+L; i++)   {
            for (int j=0; j<Ncont+L; j++)   {
                sigma->setval(i,j,0);
            }
            sigma->setval(i,i,0.1);
        }
        */

        process = new MultivariateBrownianTreeProcess(*chronogram, *sigma, rootmean, rootvar);
        for (int i=0; i<Ncont; i++)	{
            process->SetAndClamp(*contdata, L+i, i);
        }

        branchlength = new MVBranchExpoLengthArray(*process, *chronogram, dSindex);
        branchomega = new MVBranchExpoMeanArray(*process, omindex);

        cerr << "total length : " << branchlength->GetTotalLength() << '\n';
        cerr << "mean omega   : " << branchomega->GetMean() << '\n';

        nuds = nuom = 1.0;
        wnds = wnom = 0;
        if (wndsmode)   {
            wnds = new ChronoGammaWhiteNoise(*tree, *chronogram, nuds, 3);
            wndssuffstatbrancharray = new PoissonSuffStatBranchArray(*tree);
        }
        if (wnommode)   {
            wnom = new ChronoGammaWhiteNoise(*tree, *chronogram, nuom, 3);
            wnomsuffstatbrancharray = new PoissonSuffStatBranchArray(*tree);
        }

        browniansuffstat = new MultivariateNormalSuffStat(process->GetDim());

        dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);
        ifstream is(dsomsuffstatfile.c_str());
        is >> *dsompathsuffstatarray;

        cerr << "allocate ok\n";
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

    const MultivariateBrownianTreeProcess& GetProcess() const {
        return *process;
    }

    int GetNcont() const    {
        return Ncont;
    }

    int GetDim() const  {
        return sigma->GetDim();
    }

    const CovMatrix& GetCovMatrix() const   {
        return *sigma;
    }

    const Chronogram& GetChronogram() const   {
        return *chronogram;
    }

    const Link* GetRoot() const {
        return tree->GetRoot();
    }

    void NoUpdate() {}

    void Update() override {
        branchlength->Update();
        branchomega->Update();
        if (wndsmode)   {
            wnds->SetShape(1.0 / nuds);
        }
        if (wnommode)  {
            wnom->SetShape(1.0 / nuom);
        }
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
        if (wndsmode)   {
            total += WNdSHyperLogPrior();
            total += WNdSLogPrior();
        }
        if (wnommode)   {
            total += WNOmHyperLogPrior();
            total += WNOmLogPrior();
        }
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

    double WNdSHyperLogPrior() const    {
        return -nuds/100;
    }

    double WNdSLogPrior() const {
        return wnds->GetLogProb();
    }

    double WNOmHyperLogPrior() const    {
        return -nuom/100;
    }

    double WNOmLogPrior() const {
        return wnom->GetLogProb();
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    double dSOmPathSuffStatLogProb() const {
        double total = 0;
        for (int index=0; index<tree->GetNbranch(); index++)    {
            double bl = branchlength->GetVal(index);
            if (wndsmode)   {
                bl *= wnds->GetVal(index);
            }
            double om = branchomega->GetVal(index);
            if (wnommode)   {
                om *= wnom->GetVal(index);
            }
            total += dsompathsuffstatarray->GetVal(index).GetLogProb(bl, om);
        }
        return total;
    }

    double BranchSuffStatLogProb(const Link*from) const  {
        if (from->isRoot())   {
            return 0;
        }
        int index = from->GetBranch()->GetIndex();
        double bl = branchlength->GetVal(index);
        if (wndsmode)   {
            bl*= wnds->GetVal(index);
        }
        double om = branchomega->GetVal(index);
        if (wnommode)   {
            om *= wnom->GetVal(index);
        }
        return dsompathsuffstatarray->GetVal(index).GetLogProb(bl, om);
    }

    double NodeSuffStatLogProb(const Link* from) const {
        double total = BranchSuffStatLogProb(from);
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += BranchSuffStatLogProb(link->Out());
        }
        return total;
    }

    double BranchSuffStatLogProbdSIntegrated(const Link*from) const  {
        if (from->isRoot())   {
            return 0;
        }
        int index = from->GetBranch()->GetIndex();
        double om = branchomega->GetVal(index);
        if (wnommode)   {
            om *= wnom->GetVal(index);
        }
        return dsompathsuffstatarray->GetVal(index).GetLogProbdSIntegrated(branchlength->GetVal(index), om, chronogram->GetDeltaTime(from), nuds);
    }

    double NodeSuffStatLogProbdSIntegrated(const Link* from) const {
        double total = BranchSuffStatLogProbdSIntegrated(from);
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += BranchSuffStatLogProbdSIntegrated(link->Out());
        }
        return total;
    }

    double BranchSuffStatLogProbOmIntegrated(const Link*from) const  {
        if (from->isRoot())   {
            return 0;
        }
        int index = from->GetBranch()->GetIndex();
        double bl = branchlength->GetVal(index);
        if (wndsmode)   {
            bl*= wnds->GetVal(index);
        }
        return dsompathsuffstatarray->GetVal(index).GetLogProbOmIntegrated(bl, branchomega->GetVal(index), chronogram->GetDeltaTime(from), nuom);
    }

    double NodeSuffStatLogProbOmIntegrated(const Link* from) const {
        double total = BranchSuffStatLogProbOmIntegrated(from);
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += BranchSuffStatLogProbOmIntegrated(link->Out());
        }
        return total;
    }

    double NodeLogPrior(const Link* from) const  {
        return process->GetNodeLogProb(from);
    }

    double NodeLogProb(const Link* from) const   {
        return NodeLogPrior(from) + NodeSuffStatLogProb(from);
    }

    double NodeLogProbdSIntegrated(const Link* from) const   {
        return NodeLogPrior(from) + NodeSuffStatLogProbdSIntegrated(from);
    }

    double NodeLogProbOmIntegrated(const Link* from) const   {
        return NodeLogPrior(from) + NodeSuffStatLogProbOmIntegrated(from);
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

    void NudSUpdate()   {
        wnds->SetShape(1.0 / nuds);
    }

    void NuOmUpdate()   {
        wnom->SetShape(1.0 / nuom);
    }

    double NudSHyperLogProb() const {
        return WNdSHyperLogPrior() + WNdSLogPrior();
    }

    double NuOmHyperLogProb() const {
        return WNOmHyperLogPrior() + WNOmLogPrior();
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
            if (wndsmode)   {
                MoveNudS();
            }
            if (wnommode)   {
                MoveNuOm();
            }
            MoveSigma();
            MoveKappa();
        }
    }

    // Times and Rates

    void MoveTimes()    {
        if (wndsmode)   {
            chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );
            ResampleWNdS();
        }
        else    {
            chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
        }
    }

    void MoveBrownianProcess()  {

        if (wndsmode)   {
            process->SingleNodeMove(0, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );
            process->SingleNodeMove(0, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );
            ResampleWNdS();
        }
        else    {
            process->SingleNodeMove(0, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
            process->SingleNodeMove(0, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
        }

        if (wnommode)   {
            process->SingleNodeMove(1, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbOmIntegrated(from);} );
            process->SingleNodeMove(1, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbOmIntegrated(from);} );
            ResampleWNOm();
        }
        else    {
            process->SingleNodeMove(1, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
            process->SingleNodeMove(1, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
        }

        for (int i=L; i<L+Ncont; i++)   {
            process->SingleNodeMove(i, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogPrior(from);} );
            process->SingleNodeMove(i, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogPrior(from);} );
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

    void MoveNudS()  {
        ScalingMove(nuds, 1.0, 10, &FastCoevolModel::NudSHyperLogProb, &FastCoevolModel::NudSUpdate, this);
        ScalingMove(nuds, 0.3, 10, &FastCoevolModel::NudSHyperLogProb, &FastCoevolModel::NudSUpdate, this);
    }

    void MoveNuOm()  {
        ScalingMove(nuom, 1.0, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
        ScalingMove(nuom, 0.3, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
    }

    void ResampleWNdS() {
        wndssuffstatbrancharray->Clear();
        dsompathsuffstatarray->AddWNdSSuffStat(*wndssuffstatbrancharray, *branchlength, *branchomega, *wnom);
        wnds->GibbsResample(*wndssuffstatbrancharray);
    }

    void ResampleWNOm() {
        wnomsuffstatbrancharray->Clear();
        dsompathsuffstatarray->AddWNOmSuffStat(*wnomsuffstatbrancharray, *branchlength, *branchomega, *wnds);
        wnom->GibbsResample(*wnomsuffstatbrancharray);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    void PrintEntries(ostream& os) const   {
        os << "dS\n";
        os << "dN/dS\n";
        for (int i=0; i<GetNcont(); i++)    {
            os << contdata->GetCharacterName(i) << '\n';
        }
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL";
        os << "\tlength";
        os << "\tmeanomega";
        if (wndsmode)   {
            os << "\tnuds";
        }
        if (wnommode)   {
            os << "\tnuom";
        }
        for (int i=0; i<process->GetDim(); i++) {
            for (int j=i+1; j<process->GetDim(); j++)   {
                os << "\ts_" << i << "_" << j;
            }
        }
        for (int i=0; i<process->GetDim(); i++) {
            os << "\ts_" << i << "_" << i;
        }
        for (int i=0; i<process->GetDim(); i++) {
            os << "\tk_" << i;
        }
        os << '\n';
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength();
        os << '\t' << branchomega->GetMean();
        if (wndsmode)   {
            os << '\t' << nuds;
        }
        if (wnommode)   {
            os << '\t' << nuom;
        }
        for (int i=0; i<process->GetDim(); i++) {
            for (int j=i+1; j<process->GetDim(); j++)   {
                os << '\t' << (*sigma)(i,j);
            }
        }
        for (int i=0; i<process->GetDim(); i++) {
            os << '\t' << (*sigma)(i,i);
        }
        for (int i=0; i<process->GetDim(); i++) {
            os << '\t' << kappa[i];
        }
        os << '\n';
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << *chronogram;
        os << '\t' << kappa;
        os << '\t' << *sigma;
        os << '\t' << *process;
        if (wndsmode)   {
            os << '\t' << nuds ;
            os << '\t' << *wnds;
        }
        if (wnommode)   {
            os << '\t' << nuom;
            os << '\t' << *wnom;
        }
        os << '\n';
    }

    void FromStream(istream &is) override {
        is >> *chronogram;
        is >> kappa;
        is >> *sigma;
        is >> *process;
        if (wndsmode)   {
            is >> nuds;
            is >> *wnds;
        }
        if (wnommode)   {
            is >> nuom;
            is >> *wnom;
        }
    }
};
