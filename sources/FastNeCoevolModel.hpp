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
// X_3: log tau
// theta = 4 Ne u = 4 exp(X_1 + X_2)
// pnps = A Ne^-alpha = A * exp(-alpha * X_1)
// dnds = 

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

class TipNe : public SimpleArray<double> {

    public:

    TipNe(const NodeSelector<vector<double>>& innodetree, int inNe_idx, double ininvshape) :
        SimpleArray<double>((innodetree.GetNnode()+1)/2,0),
        nodetree(innodetree),
        Ne_idx(inNe_idx),
        invshape(ininvshape) {
            Sample();
    }

    ~TipNe() {}

    void SetInvShape(double in) {
        invshape = in;
    }

    double GetMeanLog() const   {
        double tot = 0;
        for (int i=0; i<GetSize(); i++) {
            tot += log(GetVal(i));
        }
        return tot / GetSize();
    }

    void Sample()   {
        for (int i=0; i<GetSize(); i++) {
            Sample(i);
        }
    }

    void Sample(int i)  {
        double longtermNe = exp(nodetree.GetVal(i)[Ne_idx]);
        double shape = 1.0 / invshape;
        double scale = shape / longtermNe;
        (*this)[i] = Random::Gamma(shape, scale);
    }

    double GetLogProb() {
        double tot = 0;
        for (int i=0; i<GetSize(); i++) {
            tot += GetLogProb(i);
        }
        return tot;
    }

    double GetLogProb(int i)    {
        double longtermNe = exp(nodetree.GetVal(i)[Ne_idx]);
        double shape = 1.0 / invshape;
        double scale = shape / longtermNe;
        return Random::logGammaDensity(GetVal(i), shape, scale);
    }

    // Ns ~ Gamma(mean = Nl, shape = a)
    // Ns/Nl ~ Gamma(mean = 1, shape = a)
    void AddSuffStat(GammaSuffStat& suffstat) const {
        for (int i=0; i<GetSize(); i++) {
            double longtermNe = exp(nodetree.GetVal(i)[Ne_idx]);
            double x = GetVal(i) / longtermNe;
            suffstat.AddSuffStat(x, log(x));
        }
    }

    template<class Update, class LogProb> 
    double ScalingMove(double tuning, int nrep, Update update, LogProb logprob) {
        double acc = 0;
        for (int i=0; i<GetSize(); i++) {
            acc += ScalingMove(i, tuning, nrep, update, logprob);
        }
        return acc/GetSize()/nrep;
    }

    template<class Update, class LogProb> 
    double ScalingMove(int i, double tuning, int nrep, Update update, LogProb logprob)  {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = - GetLogProb(i) - logprob(i);
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            (*this)[i] *= e;
            update(i);
            deltalogprob += GetLogProb(i) + logprob(i);
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                (*this)[i] /= e;
                update(i);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    private:

    const NodeSelector<vector<double>>& nodetree;
    int Ne_idx;
    double invshape;
};

class TipTheta {

    public:

    TipTheta(const Selector<double>& intipNe, const NodeSelector<vector<double>>& innodetree, int inu_idx, double inu_scale) : 
        tipNe(intipNe),
        nodetree(innodetree),
        u_idx(inu_idx), u_scale(inu_scale) {
    }

    ~TipTheta() {}

    double GetVal(int node_index) const {
        double ret = 4 * u_scale * exp(nodetree.GetVal(node_index)[u_idx]) * tipNe.GetVal(node_index);
        if (! ret)  {
            cerr << "in tiptheta: null value\n";
            cerr << u_scale << '\t' << nodetree.GetVal(node_index)[u_idx] << '\t' << tipNe.GetVal(node_index) << '\n';
            exit(1);
        }
        return ret;
    }

    private:

    const Selector<double>& tipNe;
    const NodeSelector<vector<double>>& nodetree;
    int u_idx;
    double u_scale;
};

class TipGamma {

    public:

    TipGamma(const Selector<double>& intipNe, double inA, double inalpha) :
        tipNe(intipNe),
        A(inA), alpha(inalpha) {
    }

    ~TipGamma() {}

    void SetParameters(double inA, double inalpha)  {
        A = inA;
        alpha = inalpha;
    }

    double GetVal(int node_index) const {
        double ret = A * exp(-alpha *  log(tipNe.GetVal(node_index)));
        if (! ret)  {
            cerr << "tip gamma: null value\n";
            cerr << tipNe.GetVal(node_index) << '\t' << A << '\t' << alpha << '\n';
            cerr << alpha * log(tipNe.GetVal(node_index)) << '\n';
            exit(1);
        }
        return ret;
    }

    private:

    const Selector<double>& tipNe;
    double A, alpha;
};

class pNpS  {

    public:

    pNpS(const Tree& intree, const IntegerData& from, const Selector<double>& inerror, 
            const TipTheta& intiptheta, const TipGamma& intipgamma) :
        tree(intree),
        Ntaxa(from.GetNtaxa()),
        Ks(Ntaxa, 0),
        Kn(Ntaxa, 0),
        Ls(Ntaxa, 0),
        Ln(Ntaxa, 0),
        error(inerror),
        tiptheta(intiptheta),
        tipgamma(intipgamma) {
            RegisterWithData(from);
    }

    ~pNpS() {
    }

    int GetNtaxa() const    {
        return Ntaxa;
    }

    void RegisterWithData(const IntegerData& data)   {
        int k = 0;
        int n = 0;
        RecursiveRegisterWithData(tree.GetRoot(), data, k, n);
		cerr << " polymorphism data: " << n-k << " out of " << n << " missing\n";
    }

    void RecursiveRegisterWithData(const Link* from, const IntegerData& data, int& k, int& n)   {
		if(from->isLeaf()){
			n++;
			int tax = data.GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			if (tax != -1)	{
                int idx = from->GetNode()->GetIndex();
				Ks[idx] = data.GetState(tax, 0);
				Ls[idx] = data.GetState(tax, 1);
				Kn[idx] = data.GetState(tax, 2);
				Ln[idx] = data.GetState(tax, 3);
                // cerr << from->GetNode()->GetName() << '\t' << double(Ks[idx]) / double(Ls[idx]) << '\t' << double(Kn[idx]) / double(Ln[idx]) / (double(Ks[idx]) / double(Ls[idx])) << '\n';
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
        double f = error.GetVal(tax);
        double theta = tiptheta.GetVal(tax);
        double gamma = tipgamma.GetVal(tax);
        ret += -Ls[tax]*f*theta + Ks[tax]*log(Ls[tax]*f*theta);
        // ret += -Ln[tax]*f*theta*gamma + Kn[tax]*log(Ln[tax]*f*theta*gamma);
        if (std::isinf(ret))    {
            cerr << "in pnps get log prob: inf\n";
            cerr << f << '\t' << theta << '\t' << gamma << '\n';
            exit(1);
        }
        return ret;
    }

    const Tree& tree;
    int Ntaxa;
    vector<int> Ks;
    vector<int> Kn;
    vector<int> Ls;
    vector<int> Ln;
    const Selector<double>& error;
    const TipTheta& tiptheta;
    const TipGamma& tipgamma;
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
    double u_scale;
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

    // short-term fluctuations of Ne
    double tipNe_invshape;
    TipNe* tipNe;
    GammaSuffStat tipNe_hypersuffstat;

    // uncertainty about number of callable positions
    double callable_invshape;
    IIDGamma* callable_error; 
    GammaSuffStat callable_hypersuffstat;

    // dN/dS = A N_e ^ -alpha 
    double macro_A;
    double macro_alpha;

    // pN/pS = A N_e ^ -alpha
    double micro_A;
    double micro_alpha;

    TipTheta* tiptheta;
    TipGamma* tipgamma;

    double nuom;
    double nuom2;
    double nuom3;
    ChronoGammaWhiteNoise* wnom;
    PoissonSuffStatBranchArray* wnomsuffstatbrancharray;

    dSOmegaPathSuffStatBranchArray* dsompathsuffstatarray;
    MultivariateNormalSuffStat* browniansuffstat;

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

        // check that
        root_age = inroot_age;
        u_scale = 1.0;
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
        tipNe = new TipNe(*process, Ne_idx, tipNe_invshape);

        micro_A = 1.0;
        micro_alpha = 1.0;

        tiptheta = new TipTheta(*tipNe, *process, u_idx, u_scale);
        tipgamma = new TipGamma(*tipNe, micro_A, micro_alpha);

        callable_invshape = 1.0;
        callable_error = new IIDGamma(Ntaxa, 1.0, 1.0);
        for (int i=0; i<Ntaxa; i++) {
            (*callable_error)[i] = 1.0;
        }

        pnps = new pNpS(*tree, *polydata, *callable_error, *tiptheta, *tipgamma);

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
        tipgamma->SetParameters(micro_A, micro_alpha);
        tipNe->SetInvShape(1.0/tipNe_invshape);
        callable_error->SetShape(1.0 / callable_invshape);
        callable_error->SetScale(1.0 / callable_invshape);
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
        total += MacroHyperLogPrior();
        total += MicroHyperLogPrior();
        if (std::isinf(total))  {
            cerr << "log prior is inf before tip\n";
            exit(1);
        }
        total += TipNeHyperLogPrior();
        total += TipNeLogPrior();
        if (std::isinf(total))  {
            cerr << "log prior is before callable\n";
            exit(1);
        }
        total += CallableHyperLogPrior();
        total += CallableLogPrior();
        if (std::isinf(total))  {
            cerr << "log prior is inf before wn\n";
            exit(1);
        }
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

    double TipNeHyperLogPrior() const   {
        return - tipNe_invshape;
    }

    double TipNeLogPrior() const    {
        return tipNe->GetLogProb();
    }

    double MacroHyperLogPrior() const   {
        return - macro_A - 10*macro_alpha;
    }

    double MicroHyperLogPrior() const   {
        return  - micro_A - 10*micro_alpha;
    }

    double CallableHyperLogPrior() const    {
        return -callable_invshape;
    }

    double CallableLogPrior() const {
        return callable_error->GetLogProb();
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
        return TipNeHyperLogPrior() + TipNeHyperSuffStatLogProb();
    }

    double TipNeHyperSuffStatLogProb()  const   {
        double alpha = 1.0 / tipNe_invshape;
        return tipNe_hypersuffstat.GetLogProb(alpha, alpha);
    }

    // when moving tip callable_errors or tip Ne
    double PnPsLogProb(int tax) const   {
        return pnps->GetLogProb(tax);
    }

    // when moving callable error invshape
    double CallableHyperLogProb() const {
        return CallableHyperLogPrior() + CallableHyperSuffStatLogProb();
    }

    double CallableHyperSuffStatLogProb() const {
        double alpha = 1.0 / callable_invshape;
        return callable_hypersuffstat.GetLogProb(alpha, alpha);
    }

    // when moving Brownian long-term Ne
    double TipNeLogProb(const Link* from) const {
        if (! from->isLeaf())   {
            return 0;
        }
        return tipNe->GetLogProb(from->GetNode()->GetIndex());
    }

    // when  moving Brownian u: changes theta
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

    // when moving micro_A and micro_alpha
    double MicroHyperLogProb() const    {
        return MicroHyperLogPrior() + pnps->GetLogProb();
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
            MoveTipNe();
            MoveTipNeHyper();
	    /*
            MoveCallableErrors(1.0, 10);
            MoveCallableErrors(0.1, 10);
            MoveCallableErrorsHyper();
	    */
            MoveMacroHyper();
            MoveMicroHyper();
        }
    }

    void MoveTipNe()    {
        tipNe->ScalingMove(1.0, 10, [] (int i) {}, [this] (int i) {return PnPsLogProb(i);});
        tipNe->ScalingMove(0.1, 10, [] (int i) {}, [this] (int i) {return PnPsLogProb(i);});
    }

    void MoveTipNeHyper()   {
        tipNe_hypersuffstat.Clear();
        tipNe->AddSuffStat(tipNe_hypersuffstat);
        ScalingMove(tipNe_invshape, 1.0, 10, &FastCoevolModel::TipNeHyperLogProb,
                    &FastCoevolModel::NoUpdate, this);
        ScalingMove(tipNe_invshape, 0.1, 10, &FastCoevolModel::TipNeHyperLogProb,
                    &FastCoevolModel::NoUpdate, this);
        UpdateMicro();

    }

    void MoveCallableErrors(double tuning, int nrep)   {
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<Ntaxa; i++) {
                double logprob1 = PnPsLogProb(i);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*callable_error)[i] *= e;
                double logprob2 = PnPsLogProb(i);
                double loghastings = m;
                double deltalogprob = logprob2 - logprob1 + loghastings;
                int accept = (log(Random::Uniform()) < deltalogprob);
                if (! accept)   {
                    (*callable_error)[i] /= e;
                }
            }
        }
    }

    void MoveCallableErrorsHyper()  {
        callable_hypersuffstat.Clear();
        callable_hypersuffstat.AddSuffStat(*callable_error);
        ScalingMove(callable_invshape, 1.0, 10, &FastCoevolModel::CallableHyperLogProb,
                    &FastCoevolModel::NoUpdate, this);
        UpdateMicro();
    }

    void MoveMicroHyper()    {
        // move micro_A and alpha
        ScalingMove(micro_A, 1.0, 10, &FastCoevolModel::MicroHyperLogProb, &FastCoevolModel::UpdateMicro, this);
        ScalingMove(micro_alpha, 1.0, 10, &FastCoevolModel::MicroHyperLogProb, &FastCoevolModel::UpdateMicro, this);
    }

    void MoveMacroHyper()    {
        // move macro_A and alpha
        ScalingMove(macro_A, 1.0, 10, &FastCoevolModel::MacroHyperLogProb, &FastCoevolModel::UpdateMacro, this);
        // ScalingMove(macro_alpha, 1.0, 10, &FastCoevolModel::MacroHyperLogProb, &FastCoevolModel::UpdateMacro, this);
    }

    // Times and Rates

    void MoveTimes()    {
        if (wndsmode)   {
           if (wndsmode == 1)   {
                chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );
           }
           else {
                chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from) + wnom->GetBranchLogProb(from);} );
           }
           ResampleWNdS();
        }
        else    {
            chronogram->MoveTimes([this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from);} );
        }
    }

    void MoveBrownianProcess()  {

        // 0 : long-term Ne
        // 1 : u
        // 2 : gen time
        //
        // long-term Ne impacts: log prob of short term Ne, value of dN/dS

        // gen time impacts: value of dS
        // u impacts: value of dS, log prob of pS 

        if (wndsmode)   {
            process->SingleNodeMove(1, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from) + PnPsLogProb(from);} );
            process->SingleNodeMove(1, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from) + PnPsLogProb(from);} );

            process->SingleNodeMove(2, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );
            process->SingleNodeMove(2, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbdSIntegrated(from);} );

            ResampleWNdS();
        }
        else    {
            process->SingleNodeMove(1, 0.01, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from) + PnPsLogProb(from);} );
            process->SingleNodeMove(1, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from) + PnPsLogProb(from);} );
            process->SingleNodeMove(1, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from) + PnPsLogProb(from);} );

            /*
            if (Random::Uniform() < 0.3)    {
                dsacc1 += process->FilterMove(0, 10, 0, 1,
                        [this] (const Link* from) {BranchUpdate(from);},
                        [this] (const Link* from) {return BranchSuffStatLogProb(from);} );
                dsacc01 += process->FilterMove(0, 10, 0, 0.1, 
                        [this] (const Link* from) {BranchUpdate(from);},
                        [this] (const Link* from) {return BranchSuffStatLogProb(from);} );
                dsacc001 += process->FilterMove(0, 10, 0, 0.01, 
                        [this] (const Link* from) {BranchUpdate(from);},
                        [this] (const Link* from) {return BranchSuffStatLogProb(from);} );
                dsntry ++;
            }
            */
        }

        if (wnommode)   {
            process->SingleNodeMove(0, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbOmIntegrated(from) + TipNeLogProb(from);} );
            process->SingleNodeMove(0, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProbOmIntegrated(from) + TipNeLogProb(from);} );

            ResampleWNOm();
        }
        else    {
            process->SingleNodeMove(0, 0.01, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from) + TipNeLogProb(from);} );
            process->SingleNodeMove(0, 0.1, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from) + TipNeLogProb(from);} );
            process->SingleNodeMove(0, 1.0, [this](const Link* from) {NodeUpdate(from);}, [this](const Link* from) {return NodeLogProb(from) + TipNeLogProb(from);} );

            /*
            if (Random::Uniform() < 0.3)    {
                omacc1 += process->FilterMove(1, 10, 0, 1, 
                        [this] (const Link* from) {BranchUpdate(from);},
                        [this] (const Link* from) {return BranchSuffStatLogProb(from);} );
                omacc01 += process->FilterMove(1, 10, 0, 0.1, 
                        [this] (const Link* from) {BranchUpdate(from);},
                        [this] (const Link* from) {return BranchSuffStatLogProb(from);} );
                omacc001 += process->FilterMove(1, 10, 0, 0.01, 
                        [this] (const Link* from) {BranchUpdate(from);},
                        [this] (const Link* from) {return BranchSuffStatLogProb(from);} );
                omntry ++;
            }
            */
        }

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
        os << "\tmacro_A\tmacro_alpha";
        os << "\tmicro_A\tmicro_alpha";
	os << "\ttipNeinvshape";
	// os << "\tcallinvshape\tcallmean";
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
        os << '\t' << tipNe->GetMeanLog() / log(10.0);
        os << '\t' << process->GetMean(1) / log(10.0);
        os << '\t' << macro_A << '\t' << macro_alpha;
        os << '\t' << micro_A << '\t' << micro_alpha;
	os << '\t' << tipNe_invshape;
	// os << '\t' << callable_invshape << '\t' << callable_error->GetMean();
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
    }
};
