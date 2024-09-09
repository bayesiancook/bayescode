#include "ContinuousData.hpp"
#include "IntegerData.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "MultivariateBrownianTreeProcess.hpp"
#include "Chronogram.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include "InverseWishart.hpp"
#include "ChronoWhiteNoise.hpp"

// 
// X_1: log Ne
// X_2: log u
// X_3: log dN/dS
// X_4: log tau

// changing Brownian process at a tip
// changing a global parameter: requires to scan over all tips
// currently: internal indexing scheme is based on node indices

class BranchdSArray : public SimpleBranchArray<double>    {

    public:

    BranchdSArray(const NodeSelector<vector<double>>& innodetree, const NodeSelector<double>& inchrono, int inu_idx, int ingentime_idx, double inscale) :
        SimpleBranchArray<double>(innodetree.GetTree()),
        nodetree(innodetree),
        chrono(inchrono),
        u_idx(inu_idx),
        gentime_idx(ingentime_idx),
        scale(inscale) {
            Update();
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    double GetTotalLength() const {
        return RecursiveGetTotalLength(GetRoot());
    }

    double RecursiveGetTotalLength(const Link* from) const {
        double tot = 0;
        if (! from->isRoot())   {
            tot += GetVal(from->GetBranch()->GetIndex());
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            tot += RecursiveGetTotalLength(link->Out());
        }
        return tot;
    }

    void Update()   {
        RecursiveUpdate(GetRoot());
    }

    void RecursiveUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveUpdate(link->Out());
        }
    }

    void LocalUpdate(const Link* from)  {
        if (!from->isRoot()) {

            // dt is in Mya
            // rate should be in # subs per Mya = # subs per gen / # Mya per gen
            // # Mya per gen = # years per gen / 1e6
            // assuming gen time tau is measured in days:
            // dS = (u / tau * 1e6 * 365) * (dt * root_age)

            const vector<double>&  xup = nodetree.GetVal(from->GetNode()->GetIndex());
            const vector<double>& xdown = nodetree.GetVal(from->Out()->GetNode()->GetIndex());

            double synup = scale * exp(xup[u_idx] - xup[gentime_idx]);
            double syndown = scale * exp(xdown[u_idx] - xdown[gentime_idx]);
            double syn = 0.5 * (synup + syndown);
            double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            if (dt <= 0)    {
                cerr << "error: negative time on chronogram\n";
                exit(1);
            }
            (*this)[from->GetBranch()->GetIndex()] = syn * dt;
        }
    }

    void LocalNodeUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            LocalUpdate(link->Out());
        }
    }

    private:
    const NodeSelector<vector<double>>& nodetree;
    const NodeSelector<double>& chrono;
    int u_idx;
    int gentime_idx;
    double scale;
};

class BranchdNdSArray : public SimpleBranchArray<double>    {

    public:

    BranchdNdSArray(const NodeSelector<vector<double>>& innodetree, int inNe_idx, double inA, double inalpha) :
        SimpleBranchArray<double>(innodetree.GetTree()),
        nodetree(innodetree),
        Ne_idx(inNe_idx), A(inA), alpha(inalpha), om_min(1e-4)  {
            Update();
    }

    void SetParameters(double inA, double inalpha)  {
        A = inA;
        alpha = inalpha;
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    double GetMean() const  {
        return GetTotal() / GetTree().GetNbranch();
    }

    double GetTotal() const {
        return RecursiveGetTotal(GetRoot());
    }

    double RecursiveGetTotal(const Link* from) const {
        double tot = 0;
        if (! from->isRoot())   {
            tot += GetVal(from->GetBranch()->GetIndex());
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            tot += RecursiveGetTotal(link->Out());
        }
        return tot;
    }

    void Update()   {
        RecursiveUpdate(GetRoot());
    }

    void RecursiveUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveUpdate(link->Out());
        }
    }

    void LocalUpdate(const Link* from)  {
        if (!from->isRoot()) {
            double om_up = A * exp(-alpha * nodetree.GetVal(from->GetNode()->GetIndex())[Ne_idx]);
            double om_down = A * exp(-alpha * nodetree.GetVal(from->Out()->GetNode()->GetIndex())[Ne_idx]);
            double mean = 0.5 * (om_up + om_down);
            (*this)[from->GetBranch()->GetIndex()] = mean + om_min;
        }
    }

    void LocalNodeUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            LocalUpdate(link->Out());
        }
    }

    private:
    const NodeSelector<vector<double>>& nodetree;
    int Ne_idx;
    double A, alpha;
    double om_min;
};

class pNpS {

    public:

    pNpS(int inNtaxa, const NodeSelector<vector<double>>& innodetree, int inNe_idx, int inu_idx, double ininvshape, const IntegerData& from) :
        nodetree(innodetree),
        Ne_idx(inNe_idx), u_idx(inu_idx),
        Ntaxa(inNtaxa),
        Ks(Ntaxa, 0),
        Kn(Ntaxa, 0),
        Ls(Ntaxa, 0),
        Ln(Ntaxa, 0),
        taxon(Ntaxa, ""),
        invshape(ininvshape)    {
            RegisterWithData(from);
    }

    ~pNpS() {
    }

    int GetNtaxa() const    {
        return Ntaxa;
    }

    void SetInvShape(double ininvshape)    {
        invshape = ininvshape;
    }

    void RegisterWithData(const IntegerData& data)   {
        int k = 0;
        int n = 0;
        RecursiveRegisterWithData(nodetree.GetTree().GetRoot(), data, k, n);
		cerr << " polymorphism data: " << n-k << " out of " << n << " missing\n";
    }

    void RecursiveRegisterWithData(const Link* from, const IntegerData& data, int& k, int& n)   {
		if(from->isLeaf()){
			n++;
			int tax = data.GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			if (tax != -1)	{
                int idx = from->GetNode()->GetIndex();
                if (Ks[idx] != -1)  {
                    Ks[idx] = data.GetState(tax, 0);
                    Ls[idx] = data.GetState(tax, 1);
                    Kn[idx] = data.GetState(tax, 2);
                    Ln[idx] = data.GetState(tax, 3);
                }
                taxon[idx] = from->GetNode()->GetName();
                k++;
			}
			else	{
				cerr << "set and clamp : " << from->GetNode()->GetName() << " not found\n";
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegisterWithData(link->Out(), data, k, n);
		}
	}

    double GetLogProb() {
        double tot = 0;
        for (int i=0; i<GetNtaxa(); i++)    {
            tot += GetLogProb(i);
        }
        return tot;
    }

    double GetLogProb(int tax)  {
        double ret = 0;
        double Nl = exp(nodetree.GetVal(tax)[Ne_idx]);
        double u = exp(nodetree.GetVal(tax)[u_idx]);
        double alpha = 1.0 / invshape;
        double beta = alpha / Nl;
        ret += alpha * log(beta) - Random::logGamma(alpha);
        ret -= (alpha + Ks[tax]) * log(beta + 4*u*Ls[tax]) - Random::logGamma(alpha + Ks[tax]);
        if (Ks[tax])    {
            ret += Ks[tax] * log(4*u*Ls[tax]);
        }
        ret += Random::logGamma(Ks[tax] + 1);
        if (std::isinf(ret))    {
            cerr << "in pnps get log prob: inf\n";
            exit(1);
        }
        return ret;
    }

    double GetMeanLogNs() const {
        vector<double> Ns(Ntaxa,0);
        ResampleNe(Ns);
        double tot = 0;
        for (int i=0; i<Ntaxa; i++) {
            tot += log(Ns[i]);
        }
        tot /= Ntaxa;
        return tot;
    }

    void ResampleNe(vector<double>& Ns) const {
        for (int tax=0; tax<Ntaxa; tax++) {
            double shape = 1.0 / invshape + Ks[tax];
            double N = exp(nodetree.GetVal(tax)[Ne_idx]);
            double u = exp(nodetree.GetVal(tax)[u_idx]);
            double scale = shape/N + 4*u*Ls[tax];
            Ns[tax] = Random::GammaSample(shape,scale);
        }
    }

    void GetTaxonList(vector<string>& taxlist) const    {
        for (int tax=0; tax<Ntaxa; tax++) {
            taxlist[tax] = taxon[tax];
        }
    }

    void GetStats(vector<vector<double>>& stats, bool withlog, double macro_A, double macro_alpha) const   {
        for (int tax=0; tax<Ntaxa; tax++) {
            double Nl = exp(nodetree.GetVal(tax)[Ne_idx]);
            double u = exp(nodetree.GetVal(tax)[u_idx]);

            double shape0 = 1.0 / invshape;
            double scale0 = shape0/Nl;

            double shapeS = shape0 + Ks[tax];
            double scaleS = scale0 + 4*u*Ls[tax];
            double Ns = Random::GammaSample(shapeS,scaleS);

            double shapeN = shape0 + Kn[tax];
            double scaleN = scale0 + 4*u*Ln[tax];
            double ps = shapeS/scaleS;
            double pnps = (shapeN/scaleN) / (shapeS/scaleS);

            double dnds = macro_A * exp(-macro_alpha*nodetree.GetVal(tax)[Ne_idx]);

            if (withlog)    {
                stats[0][tax] = log(Nl) / log(10.0);
                stats[1][tax] = log(Ns) / log(10.0);
                stats[2][tax] = log(u) / log(10.0);
                stats[3][tax] = log(4*Ns*u) / log(10.0);
                stats[4][tax] = log(ps) / log(10.0);
                stats[5][tax] = log(pnps) / log(10.0);
                stats[6][tax] = log(dnds) / log(10.0);
            }
            else    {
                stats[0][tax] = Nl;
                stats[1][tax] = Ns;
                stats[2][tax] = u;
                stats[3][tax] = 4*Ns*u;
                stats[4][tax] = ps;
                stats[5][tax] = pnps;
                stats[6][tax] = dnds;
            }
        }
    }


    const NodeSelector<vector<double>>& nodetree;
    int Ne_idx;
    int u_idx;
    int Ntaxa;
    vector<int> Ks;
    vector<int> Kn;
    vector<int> Ls;
    vector<int> Ln;
    vector<string> taxon;
    double invshape;
};

class FastCoevolModel: public ProbModel {

    int wndsmode;
    int wnommode;

    const Tree *tree;
    const TaxonSet *taxonset;
    const ContinuousData* contdata;
    // data about number of syn and non syn targets and observed polymorphisms
    // Ks Ls Kn Ln
    const IntegerData* polydata;

    pNpS* pnps;

    string dsomsuffstatfile;

    int Ntaxa;
    int Nbranch;
    
    // includes generation time (first trait in data matrix of continuous characters)
    int Ncont;

    int u_idx, Ne_idx, gentime_idx;

    // age of the root, in Mya
    double root_age;
    // u is measured in mutation rate per gen and per base pair 
    // generations are measured in years: rate_scale is 10^6 * root_age
    double rate_scale;

    Chronogram* chronogram;

    int df;
    vector<double> kappa;
    InverseWishart* sigma;

    vector<double> rootmean;
    vector<double> rootvar;

    MultivariateBrownianTreeProcess* process;
    BranchdSArray* branchlength;
    BranchdNdSArray* branchomega;

    double nuds;
    double nuds2;
    double nuds3;
    ChronoGammaWhiteNoise* wnds;
    PoissonSuffStatBranchArray* wndssuffstatbrancharray;

    double nuom;
    double nuom2;
    double nuom3;
    ChronoGammaWhiteNoise* wnom;
    PoissonSuffStatBranchArray* wnomsuffstatbrancharray;

    dSOmegaPathSuffStatBranchArray* dsompathsuffstatarray;
    MultivariateNormalSuffStat* browniansuffstat;

    double macro_A;
    double macro_alpha;

    double tipNe_invshape;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    FastCoevolModel(string contdatafile, string polydatafile, string treefile, string rootfile, string indsomsuffstatfile, double inroot_age, int inwndsmode, int inwnommode)   {

        wndsmode = inwndsmode;
        wnommode = inwnommode;
        if (wndsmode != wnommode)   {
            cerr << "error: currently, same overdispersion mode for dS and dN/dS\n";
            exit(1);
        }

        dsomsuffstatfile = indsomsuffstatfile;

        contdata = new FileContinuousData(contdatafile);
        Ncont = contdata->GetNsite();
        Ntaxa = contdata->GetNtaxa();
        taxonset = contdata->GetTaxonSet();
        
        polydata = new FileIntegerData(polydatafile);
        
        Ne_idx = 0;
        u_idx = 1;
        gentime_idx = 2;

        root_age = inroot_age;
        rate_scale = 1e6 * 365.0 * root_age;

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
        if (dim != Ncont + 2)   {
            cerr << "error in root file: non matching dimension\n";
            cerr << dim << '\t' << Ncont << '\n';
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

        kappa.assign(Ncont+2, 1.0);
        df = 0;
        sigma = new InverseWishart(kappa, df);

        process = new MultivariateBrownianTreeProcess(*chronogram, *sigma, rootmean, rootvar);
        for (int i=0; i<Ncont; i++)	{
            process->SetAndClamp(*contdata, 2+i, i);
        }

        branchlength = new BranchdSArray(*process, *chronogram, u_idx, gentime_idx, rate_scale);

        macro_A = 1.0;
        macro_alpha = 0.08;
        branchomega = new BranchdNdSArray(*process, Ne_idx, macro_A, macro_alpha);

        cerr << "total length : " << branchlength->GetTotalLength() << '\n';
        cerr << "mean omega   : " << branchomega->GetMean() << '\n';

        nuds = nuds2 = nuds3 = nuom = nuom2 = nuom3 = 1.0;
        if (wndsmode == 2)  {
            nuds = 0.01;
        }
        if (wnommode == 2)  {
            nuom = 0.01;
        }

        wnds = wnom = 0;
        if (wndsmode)   {
            wnds = new ChronoGammaWhiteNoise(*tree, *chronogram, wndsmode);
            wndssuffstatbrancharray = new PoissonSuffStatBranchArray(*tree);
        }
        if (wnommode)   {
            wnom = new ChronoGammaWhiteNoise(*tree, *chronogram, wnommode);
            wnomsuffstatbrancharray = new PoissonSuffStatBranchArray(*tree);
        }

        tipNe_invshape = 1.0;
        pnps = new pNpS(Ntaxa, *process, Ne_idx, u_idx, tipNe_invshape, *polydata);

        browniansuffstat = new MultivariateNormalSuffStat(process->GetDim());

        dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);

        ifstream is(dsomsuffstatfile.c_str());
        string tmp;
        is >> tmp;
        if (tmp != "counts_dS") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treedscount(is);
        treedscount.SetIndices();

        is >> tmp;
        if (tmp != "counts_dS_norm") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treedsbeta(is);
        treedsbeta.SetIndices();

        is >> tmp;
        if (tmp != "counts_dN") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treedncount(is);
        treedncount.SetIndices();

        is >> tmp;
        if (tmp != "counts_dN_norm") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treednbeta(is);
        treednbeta.SetIndices();

        dsompathsuffstatarray->Add(treedscount, treedsbeta, treedncount, treednbeta);
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

    int GetNcont() const    {
        return Ncont;
    }

    int GetDim() const  {
        return sigma->GetDim();
    }

    int GetNtaxa() const    {
        return Ntaxa;
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
        UpdateMacro();
        UpdateMicro();
        if (wndsmode == 4)  {
            wnds->SetVar(nuds, nuds2, nuds3);
        }
        else if (wndsmode)   {
            wnds->SetVar(nuds);
        }
        if (wnommode == 4)  {
            wnom->SetVar(nuom, nuom2, nuom3);
        }
        else if (wnommode)  {
            wnom->SetVar(nuom);
        }
    }

    void UpdateMacro()  {
        branchlength->Update();
        branchomega->SetParameters(macro_A, macro_alpha);
        branchomega->Update();
    }

    void UpdateMicro()  {
        pnps->SetInvShape(tipNe_invshape);
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
        if (std::isinf(total))  {
            cerr << "log prior is inf before macro\n";
            exit(1);
        }
        if (std::isinf(total))  {
            cerr << "log prior is inf before tip\n";
            exit(1);
        }
        total += MacroHyperLogPrior();
        total += TipNeHyperLogPrior();
        if (wndsmode)   {
            total += WNdSHyperLogPrior();
            total += WNdSLogPrior();
        }
        if (wnommode)   {
            total += WNOmHyperLogPrior();
            total += WNOmLogPrior();
        }
        if (std::isinf(total))  {
            cerr << "log prior is inf\n";
            cerr << WNdSHyperLogPrior() << '\n';
            cerr << WNdSLogPrior() << '\n';
            cerr << WNOmHyperLogPrior() << '\n';
            cerr << WNOmLogPrior() << '\n';
            exit(1);
        }
        if (std::isnan(total))  {
            cerr << "log prior is nan\n";
            exit(1);
        }
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { 
        return GetMacroLogLikelihood() + GetMicroLogLikelihood();
    }

    double GetMacroLogLikelihood() const    {
        double ret = dSOmPathSuffStatLogProb();
        if (std::isinf(ret))    {
            cerr << "macro log likelihood is inf\n";
            exit(1);
        }
        if (std::isnan(ret))    {
            cerr << "macro log likelihood is nan\n";
            exit(1);
        }
        return ret;
    }
    
    double GetMicroLogLikelihood() const    {
        double ret = pnps->GetLogProb();
        if (std::isinf(ret))    {
            cerr << "micro log likelihood is inf\n";
            exit(1);
        }
        if (std::isnan(ret))    {
            cerr << "micro log likelihood is nan\n";
            exit(1);
        }
        return ret;
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

    double MacroHyperLogPrior() const   {
        return - macro_A/10 - macro_alpha/10;
    }

    double TipNeHyperLogPrior() const   {
        return - tipNe_invshape;
    }

    double WNdSHyperLogPrior() const    {
        if (wndsmode == 4)  {
            return - (nuds + nuds2 + nuds3)/100;
        }
        return -nuds/100;
    }

    double WNdSLogPrior() const {
        return wnds->GetLogProb();
    }

    double WNOmHyperLogPrior() const    {
        if (wnommode == 4)  {
            return - (nuom + nuom2 + nuom3)/100;
        }
        return -nuom/100;
    }

    double WNOmLogPrior() const {
        return wnom->GetLogProb();
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // when moving tipNe_invshape
    double TipNeHyperLogProb() const    {
        return TipNeHyperLogPrior() + pnps->GetLogProb();
    }

    // when moving long term Ne or u 
    double PnPsLogProb(const Link* from) const  {
        if (! from->isLeaf())   {
            return 0;
        }
        return pnps->GetLogProb(from->GetNode()->GetIndex());
    }

    // when moving macro_A and macro_alpha
    double MacroHyperLogProb() const    {
        return MacroHyperLogPrior() + SuffStatLogProbOmIntegrated(); 
    }

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
            total += dsompathsuffstatarray->GetVal(index).GetLogProb(bl, om, 1.0);
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
        return dsompathsuffstatarray->GetVal(index).GetLogProb(bl, om, 1.0);
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
        double nu = 0;
        double dt = chronogram->GetDeltaTime(from);
        switch(wndsmode)    {
            case 1:
                nu = nuds;
                break;
            case 2:
                nu = nuds/dt;
                break;
            case 3:
                nu = nuds*dt;
                break;
            case 4:
                nu = nuds/dt + nuds2 + nuds3*dt;
                break;
            default:
                cerr << "error: unrecocgnized variance mode for overdispersion\n";
                cerr << wndsmode << '\n';
                exit(1);
        }
        return dsompathsuffstatarray->GetVal(index).GetLogProbdSIntegrated(branchlength->GetVal(index), om, chronogram->GetDeltaTime(from), nu, 1.0);
    }

    double NodeSuffStatLogProbdSIntegrated(const Link* from) const {
        double total = BranchSuffStatLogProbdSIntegrated(from);
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += BranchSuffStatLogProbdSIntegrated(link->Out());
        }
        return total;
    }

    double SuffStatLogProbOmIntegrated() const  {
        return RecursiveSuffStatLogProbOmIntegrated(GetRoot());
    }

    double RecursiveSuffStatLogProbOmIntegrated(const Link* from) const {
        double total = 0;
        if (! from->isRoot())   {
            total += BranchSuffStatLogProbOmIntegrated(from);
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            total += RecursiveSuffStatLogProbOmIntegrated(link->Out());
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
        double nu = 0;
        double dt = chronogram->GetDeltaTime(from);
        switch(wnommode)    {
            case 1:
                nu = nuom;
                break;
            case 2:
                nu = nuom/dt;
                break;
            case 3:
                nu = nuom*dt;
                break;
            case 4:
                nu = nuom/dt + nuom2 + nuom3*dt;
                break;
            default:
                cerr << "error: unrecocgnized variance mode for overdispersion\n";
                cerr << wnommode << '\n';
                exit(1);
        }
        return dsompathsuffstatarray->GetVal(index).GetLogProbOmIntegrated(bl, branchomega->GetVal(index), chronogram->GetDeltaTime(from), nu, 1.0);
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

    double BranchLogPrior(const Link* from) const   {
        return process->GetLocalLogProb(from);
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

    void BranchUpdate(const Link* from) {
        branchlength->LocalUpdate(from);
        branchomega->LocalUpdate(from);
    }

    double KappaSuffStatLogProb() const {
        return sigma->GetLogProb();
    }

    double KappaLogProb() const {
        return KappaLogPrior() + KappaSuffStatLogProb();
    }

    void NudSUpdate()   {
        if (wndsmode == 4)  {
            wnds->SetVar(nuds, nuds2, nuds3);
        }
        else    {
            wnds->SetVar(nuds);
        }
    }

    void NuOmUpdate()   {
        if (wnommode == 4)  {
            wnom->SetVar(nuom, nuom2, nuom3);
        }
        else    {
            wnom->SetVar(nuom);
        }
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
        MoveParameters(300);
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
            MoveMacroHyper();
            MoveTipNeHyper();
        }
    }

    void MoveMacroHyper()    {
        // move macro_A and alpha
        ScalingMove(macro_A, 1.0, 10, &FastCoevolModel::MacroHyperLogProb, &FastCoevolModel::UpdateMacro, this);
        ScalingMove(macro_A, 0.3, 10, &FastCoevolModel::MacroHyperLogProb, &FastCoevolModel::UpdateMacro, this);
        ScalingMove(macro_alpha, 1.0, 10, &FastCoevolModel::MacroHyperLogProb, &FastCoevolModel::UpdateMacro, this);
        ScalingMove(macro_alpha, 0.3, 10, &FastCoevolModel::MacroHyperLogProb, &FastCoevolModel::UpdateMacro, this);
        ResampleWNOm();
    }

    void MoveTipNeHyper()   {
        ScalingMove(tipNe_invshape, 1.0, 10, &FastCoevolModel::TipNeHyperLogProb,
                    &FastCoevolModel::UpdateMicro, this);
        ScalingMove(tipNe_invshape, 0.1, 10, &FastCoevolModel::TipNeHyperLogProb,
                    &FastCoevolModel::UpdateMicro, this);
    }

    // Times and Rates

    void MoveTimes()    {
        if (wndsmode)   {
            chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from) + wnom->GetBranchLogProb(from);} );
           ResampleWNdS();
        }
        else    {
            cerr << "moves without white noise not yet implemented\n";
            exit(1);
        }
    }

    void MoveBrownianProcess()  {

        // 0 : long-term Ne
        // 1 : u
        // 2 : gen time
        //
        // long-term Ne impacts: log prob of pS and value of omega
        // u impacts: value of dS, log prob of pS
        // gen time impacts: value of dS

        // moving long-term Ne
        if (wnommode)   {
            process->SingleNodeMove(0, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbOmIntegrated(from) + PnPsLogProb(from);} );
            process->SingleNodeMove(0, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbOmIntegrated(from) + PnPsLogProb(from);} );
            ResampleWNOm();
        }
        else    {
            cerr << "moves without white noise not yet implemented\n";
            exit(1);
        }

        // moving u and gen time
        if (wndsmode)   {
            process->SingleNodeMove(1, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from) + PnPsLogProb(from);} );
            process->SingleNodeMove(1, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from) + PnPsLogProb(from);} );

            process->SingleNodeMove(2, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );
            process->SingleNodeMove(2, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );

            ResampleWNdS();
        }
        else    {
            cerr << "moves without white noise not yet implemented\n";
            exit(1);
        }

        // moving other quantitative traits
        for (int i=3; i<2+Ncont; i++)   {
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
        if (wndsmode == 4)  {
            ScalingMove(nuds2, 1.0, 10, &FastCoevolModel::NudSHyperLogProb, &FastCoevolModel::NudSUpdate, this);
            ScalingMove(nuds2, 0.3, 10, &FastCoevolModel::NudSHyperLogProb, &FastCoevolModel::NudSUpdate, this);
            ScalingMove(nuds3, 1.0, 10, &FastCoevolModel::NudSHyperLogProb, &FastCoevolModel::NudSUpdate, this);
            ScalingMove(nuds3, 0.3, 10, &FastCoevolModel::NudSHyperLogProb, &FastCoevolModel::NudSUpdate, this);
        }
    }

    void MoveNuOm()  {
        ScalingMove(nuom, 1.0, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
        ScalingMove(nuom, 0.3, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
        if (wnommode == 4)  {
            ScalingMove(nuom2, 1.0, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
            ScalingMove(nuom2, 0.3, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
            ScalingMove(nuom3, 1.0, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
            ScalingMove(nuom3, 0.3, 10, &FastCoevolModel::NuOmHyperLogProb, &FastCoevolModel::NuOmUpdate, this);
        }
    }

    void ResampleWNdS() {
        wndssuffstatbrancharray->Clear();
        dsompathsuffstatarray->AddWNdSSuffStat(*wndssuffstatbrancharray, *branchlength, *branchomega, *wnom);
        wnds->GibbsResample(*wndssuffstatbrancharray, 1.0);
    }

    void ResampleWNOm() {
        wnomsuffstatbrancharray->Clear();
        dsompathsuffstatarray->AddWNOmSuffStat(*wnomsuffstatbrancharray, *branchlength, *branchomega, *wnds);
        wnom->GibbsResample(*wnomsuffstatbrancharray, 1.0);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    void GetStats(vector<vector<double>>& stats, bool withlog) const  {
        pnps->GetStats(stats, withlog, macro_A, macro_alpha);
    }

    void GetTaxonList(vector<string>& taxlist) const    {
        pnps->GetTaxonList(taxlist);
    }

    void PrintEntries(ostream& os) const   {
        os << "longtermNe\n";
        os << "u\n";
        for (int i=0; i<GetNcont(); i++)    {
            os << contdata->GetCharacterName(i) << '\n';
        }
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnLmacro\tlnLmicro";
        os << "\tlength";
        os << "\tmeanomega";
        os << "\tmeanlog10Nelong";
        os << "\tmeanlog10Neshort";
        os << "\tmeanlog10u";
        os << "\tmacroA\tmacroalpha";
        os << "\ttipNeinvshape";
        if (wndsmode == 4)  {
            os << "\tds_va\tds_vb\tds_vc";
        }
        else if (wndsmode)   {
            os << "\tnuds";
        }
        if (wnommode == 4)  {
            os << "\tom_va\tom_vb\tom_vc";
        }
        else if (wnommode)   {
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

    const MultivariateBrownianTreeProcess& GetProcess() const {
        return *process;
    }

    const ChronoGammaWhiteNoise& GetSynDev() const  {
        return *wnds;
    }

    const ChronoGammaWhiteNoise& GetOmDev() const  {
        return *wnom;
    }

    void Trace(ostream &os) const override {
        os.precision(12);
        os << GetLogPrior() << '\t';
        os << GetMacroLogLikelihood() << '\t';
        os << GetMicroLogLikelihood() << '\t';
        os << branchlength->GetTotalLength();
        os << '\t' << branchomega->GetMean();
        os << '\t' << process->GetMean(0) / log(10.0);
        os << '\t' << pnps->GetMeanLogNs() / log(10.0);
        os << '\t' << process->GetMean(1) / log(10.0);
        os << '\t' << macro_A << '\t' << macro_alpha;
        os << '\t' << tipNe_invshape;
        if (wndsmode == 4)  {
            double v1,v2,v3;
            v1 = v2 = v3 = 0;
            wnds->VariancePartition(v1,v2,v3);
            double tot = v1 + v2 + v3;
            os << '\t' << nuds << '\t' << nuds2 << '\t' << nuds3;
            os << '\t' << v1/tot << '\t' << v2/tot << '\t' << v3/tot;
        }
        else if (wndsmode)   {
            os << '\t' << nuds;
        }
        if (wnommode == 4)  {
            double v1,v2,v3;
            v1 = v2 = v3 = 0;
            wnom->VariancePartition(v1,v2,v3);
            double tot = v1 + v2 + v3;
            os << '\t' << nuom << '\t' << nuom2 << '\t' << nuom3;
            os << '\t' << v1/tot << '\t' << v2/tot << '\t' << v3/tot;
        }
        else if (wnommode)   {
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

    void Monitor(ostream &os) const override {
    }

    void ToStream(ostream &os) const override {
        os << *chronogram;
        os << '\t' << kappa;
        os << '\t' << *sigma;
        os << '\t' << *process;
        if (wndsmode)   {
            os << '\t' << nuds;
            if (wndsmode == 4)  {
                os << '\t' << nuds2 << '\t' << nuds3;
            }
            os << '\t' << *wnds;
        }
        if (wnommode)   {
            os << '\t' << nuom;
            if (wnommode == 4)  {
                os << '\t' << nuom2 << '\t' << nuom3;
            }
            os << '\t' << *wnom;
        }
        os << '\t' << macro_A << '\t' << macro_alpha;
        os << '\t' << tipNe_invshape;
        os << '\n';
    }

    void FromStream(istream &is) override {
        is >> *chronogram;
        is >> kappa;
        is >> *sigma;
        is >> *process;
        if (wndsmode)   {
            is >> nuds;
            if (wndsmode == 4)  {
                is >> nuds2 >> nuds3;
            }
            is >> *wnds;
        }
        if (wnommode)   {
            is >> nuom;
            if (wnommode == 4)  {
                is >> nuom2 >> nuom3;
            }
            is >> *wnom;
        }
        is >> macro_A >> macro_alpha;
        is >> tipNe_invshape;
    }
};
