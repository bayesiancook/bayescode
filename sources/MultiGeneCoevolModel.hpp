
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "CoevolModel.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"

class MultiGeneCoevolModel : public MultiGeneProbModel {

  private:

    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;
    std::vector<CodonSequenceAlignment*> alivector;
    const ContinuousData* contdata;

    string datafile;
    string treefile;
    string contdatafile;
    string rootfile;

    GeneticCodeType codetype;

    int Ntaxa;
    int Nbranch;

    int nucmode;

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

    // Nucleotide rates

    // shared nuc rates
    GTRSubMatrix *nucmatrix;
    NucPathSuffStat nucpathsuffstat;

    // gene-specific nuc rates
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    IIDDirichlet *nucrelratearray;
    DirichletSuffStat nucrelratesuffstat;

    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    IIDDirichlet *nucstatarray;
    DirichletSuffStat nucstatsuffstat;

    // each gene defines its own CoevolModel
    std::vector<CoevolModel *> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneCoevolModel(string indatafile, string incontdatafile, string intreefile, string inrootfile, GeneticCodeType incodetype, int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {

        nucmode = 2;

        datafile = indatafile;
        treefile = intreefile;
        contdatafile = incontdatafile;
        rootfile = inrootfile;
        codetype = incodetype;

        AllocateAlignments(datafile);

        refcodondata = new CodonSequenceAlignment(refdata, true, codetype);
        taxonset = refdata->GetTaxonSet();
        Ntaxa = refdata->GetNtaxa();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

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

        if (!myid) {
            cerr << "number of taxa : " << Ntaxa << '\n';
            cerr << "number of branches : " << Nbranch << '\n';
            cerr << "tree and data fit together\n";
        }
    }

    void Allocate() {

        // Chronogram
        chronogram = new Chronogram(*tree);

        // Brownian process
        kappa.assign(Ncont+L, 1.0);
        df = 0;
        sigma = new InverseWishart(kappa, df);

        process = new MultivariateBrownianTreeProcess(*chronogram, *sigma, rootmean, rootvar);
        for (int i=0; i<Ncont; i++)	{
            process->SetAndClamp(*contdata, L+i, i);
        }

        branchlength = new MVBranchExpoLengthArray(*process, *chronogram, dSindex);
        branchomega = new MVBranchExpoMeanArray(*process, omindex);

        dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);
        browniansuffstat = new MultivariateNormalSuffStat(process->GetDim());

        cerr << "total length : " << branchlength->GetTotalLength() << '\n';
        cerr << "mean omega   : " << branchomega->GetMean() << '\n';

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;

        if (nucmode == 2) {
            nucrelratearray = new IIDDirichlet(1, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(1, nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = new GTRSubMatrix(Nnuc, (*nucrelratearray)[0], (*nucstatarray)[0], true);
        } else {
            nucrelratearray =
                new IIDDirichlet(GetLocalNgene(), nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray =
                new IIDDirichlet(GetLocalNgene(), nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = 0;
        }

        // Gene processes 

        lnL = 0;
        GeneLogPrior = 0;

        if (!GetMyid()) {
            geneprocess.assign(0, (CoevolModel*)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (CoevolModel*)0);

            ifstream is(datafile.c_str());
            string tmp;
            is >> tmp;
            if (tmp == "ALI")   {
                int ngene;
                is >> ngene;
                if (ngene != GetNgene())    {
                    cerr << "error when reading alignments from cat file: non matching number of genes\n";
                    exit(1);
                }
                alivector.assign(GetLocalNgene(), (CodonSequenceAlignment*) 0);
                int index = 0;
                for (int gene=0; gene<GetNgene(); gene++)   {
                    string name;
                    is >> name;
                    FileSequenceAlignment tmp(is);
                    if ((index < GetLocalNgene()) && (name == GeneName[index]))    {
                        if (GetLocalGeneName(index) != name)    {
                            cerr << "error: non matching gene name\n";
                            exit(1);
                        }
                        if (alivector[index]) {
                            cerr << "error: alignment already allocated\n";
                            exit(1);
                        }
                        alivector[index] = new CodonSequenceAlignment(&tmp, true, codetype);
                        index++;
                    }
                }
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    if (! alivector[gene])  {
                        cerr << "error: alignment not allocated\n";
                        exit(1);
                    }
                    geneprocess[gene] = new CoevolModel(alivector[gene], contdata, tree, rootmean, rootvar);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    geneprocess[gene] = new CoevolModel(GetLocalGeneName(gene), contdatafile, treefile, rootfile, "None", "None", codetype);
                }
            }

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->SetCoevolMode(2);
                geneprocess[gene]->Allocate();
            }
        }
    }

    void SetNucMode(int innucmode)   {
        nucmode = innucmode;
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

    void FastUpdate() {
        branchlength->Update();
        branchomega->Update();
        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
        if (nucmode == 2)   {
            UpdateNucMatrix();
        }
    }

    void MasterUpdate() override {

        FastUpdate();

        if (nprocs > 1) {
            MasterSendNucRatesHyperParameters();
            if (nucmode == 2) {
                MasterSendGlobalNucRates();
            } else {
                MasterSendGeneNucRates();
            }
            MasterSendCoevol();
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {

        SlaveReceiveNucRatesHyperParameters();
        if (nucmode == 2) {
            SlaveReceiveGlobalNucRates();
        } else {
            SlaveReceiveGeneNucRates();
        }
        SlaveReceiveCoevol();
        GeneUpdate();
        SlaveSendLogProbs();
    }

    void GeneUpdate() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->Update();
        }
    }

    void MasterPostPred(string name) override {
        FastUpdate();
        if (nprocs > 1) {
            MasterSendNucRatesHyperParameters();
            if (nucmode == 2) {
                MasterSendGlobalNucRates();
            } else {
                MasterSendGeneNucRates();
            }
            MasterSendCoevol();
        }
    }

    void SlavePostPred(string name) override {
        SlaveReceiveNucRatesHyperParameters();
        if (nucmode == 2) {
            SlaveReceiveGlobalNucRates();
        } else {
            SlaveReceiveGeneNucRates();
        }
        SlaveReceiveCoevol();
        GenePostPred(name);
    }

    void GenePostPred(string name) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->PostPred(name + GetLocalGeneName(gene));
        }
    }

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    //-------------------
    // Traces and Monitors
    //-------------------

    void PrintEntries(ostream& os) const   {
        os << "dS\n";
        os << "dN/dS\n";
        for (int i=0; i<GetNcont(); i++)    {
            os << contdata->GetCharacterName(i) << '\n';
        }
    }

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
        os << "\tstatent";
        os << "\trrent";
        if (nucmode != 2) {
            os << "\tstdevrr\tcenter\thyperinvconc";
            os << "\tstdevstat\tcenter\thyperinvconc";
        }
        os << '\n';
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

        os << '\t' << nucstatarray->GetMeanEntropy();
        os << '\t' << nucrelratearray->GetMeanEntropy();
        if (nucmode != 2) {
            os << '\t' << sqrt(GetVarNucRelRate()) << '\t' << Random::GetEntropy(nucrelratehypercenter)
               << '\t' << nucrelratehyperinvconc;
            os << '\t' << sqrt(GetVarNucStat()) << '\t' << Random::GetEntropy(nucstathypercenter)
               << '\t' << nucstathyperinvconc;
        }
        os << '\n';
        os.flush();
    }

    // Nucleotide rates

    double GetVarNucRelRate() const {
        if (nucmode == 2) {
            cerr << "error in getvarnucrelrate\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nrr; j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < Ngene; g++) {
                double tmp = (*nucrelratearray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= Ngene;
            var /= Ngene;
            var -= mean * mean;
            tot += var;
        }
        tot /= Nrr;
        return tot;
    }

    double GetVarNucStat() const {
        if (nucmode == 2) {
            cerr << "error in getvarnucstat\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nnuc; j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < Ngene; g++) {
                double tmp = (*nucstatarray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= Ngene;
            var /= Ngene;
            var -= mean * mean;
            tot += var;
        }
        tot /= Nnuc;
        return tot;
    }

    void Monitor(ostream &os) const override {}

    void MasterFromStream(istream &is) override {

        is >> *chronogram;
        is >> kappa;
        is >> *sigma;
        is >> *process;

        is >> nucrelratehypercenter;
        is >> nucrelratehyperinvconc;
        is >> nucstathypercenter;
        is >> nucstathyperinvconc;
        is >> *nucrelratearray;
        is >> *nucstatarray;
    }

    void MasterToStream(ostream &os) const override {

        os << *chronogram << '\t';
        os << kappa << '\t';
        os << *sigma << '\t';
        os << *process << '\n';

        os << nucrelratehypercenter << '\t';
        os << nucrelratehyperinvconc << '\t';
        os << nucstathypercenter << '\t';
        os << nucstathyperinvconc << '\t';
        os << *nucrelratearray << '\t';
        os << *nucstatarray << '\t';
    }

    //-------------------
    // Updates
    //-------------------

    void UpdateNucMatrix() {
        nucmatrix->CopyStationary((*nucstatarray)[0]);
        nucmatrix->CorruptMatrix();
    }

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
        // gene contributions
        double total = GeneLogPrior;

        // nuc rates
        if (nucmode == 2) {
            total += GlobalNucRatesLogPrior();
        } else if (nucmode == 1) {
            total += GeneNucRatesHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        total += ChronoLogPrior();
        total += KappaLogPrior();
        total += SigmaLogPrior();
        total += BrownianProcessLogPrior();

        return total;
    }

    // Nucleotide rates

    double GlobalNucRatesLogPrior() const {
        return nucrelratearray->GetLogProb() + nucstatarray->GetLogProb();
    }

    // exponential of mean 1 for nucrelrate and nucstat hyper inverse
    // concentration
    double GeneNucRatesHyperLogPrior() const {
        double total = 0;
        if (nucmode == 1) {
            total -= nucrelratehyperinvconc;
            total -= nucstathyperinvconc;
        }
        return total;
    }

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

    double GetLogLikelihood() const { return lnL; }

    //-------------------
    // Nuc Rates Log Probs
    //-------------------

    // Nucleotide rates

    // suff stat for global nuc rates, as a function of nucleotide matrix
    // (which itself depends on nucstat and nucrelrate)
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    // suff stat for gene-specific nuc rates, as a function of nucrate
    // hyperparameters
    double NucRatesHyperSuffStatLogProb() const {
        double total = 0;
        total += nucrelratesuffstat.GetLogProb(nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total += nucstatsuffstat.GetLogProb(nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    // log prob for moving nuc rates hyper params
    double NucRatesHyperLogProb() const {
        return GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates
    double NucRatesLogProb() const { return GlobalNucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    //-------------------
    // Brownian Log Probs
    //-------------------

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
    // Moves
    //-------------------

    // all methods starting with Master are called only by master
    // for each such method, there is a corresponding method called by slave, and
    // starting with Slave
    //
    // all methods starting with Gene are called only be slaves, and do some work
    // across all genes allocated to that slave

    void MasterMove() override {
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                MasterReceiveNucPathSuffStat();
                MoveNucRates();
                MasterSendGlobalNucRates();
            } else if (nucmode == 1) {
                MasterReceiveNucRatesHyperSuffStat();
                MoveNucRatesHyperParameters();
                MasterSendNucRatesHyperParameters();
            }

            MasterReceivedSOmegaSuffStat();
            MoveCoevol();
            MasterSendCoevol();

        }

        // collect current state
        if (nucmode != 2) {
            MasterReceiveGeneNucRates();
        }
        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {
        GeneResampleSub(1.0);

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {

            GeneCollectPathSuffStat();

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                SlaveSendNucPathSuffStat();
                SlaveReceiveGlobalNucRates();
            } else if (nucmode == 1) {
                GeneMoveNucRates();
                SlaveSendNucRatesHyperSuffStat();
                SlaveReceiveNucRatesHyperParameters();
            }

            SlaveSenddSOmegaSuffStat();
            SlaveReceiveCoevol();
        }

        // collect current state
        if (nucmode != 2) {
            SlaveSendGeneNucRates();
        }
        SlaveSendLogProbs();
    }

    void GeneResampleSub(double frac) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void GeneCollectPathSuffStat()  {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectPathSuffStat();
        }
    }

    void GeneMoveNucRates() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveNuc();
        }
    }

    /*
    void MoveGeneParameters(int nrep)   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveParameters(nrep);
            if (nucmode != 2) {
                geneprocess[gene]->GetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
            }
        }
    }
    */

    // Nucleotide rates

    void MoveNucRatesHyperParameters() {
        ProfileMove(nucrelratehypercenter, 1.0, 1, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10, &MultiGeneCoevolModel::NucRatesHyperLogProb,
                    &MultiGeneCoevolModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveNucRates() {
        vector<double> &nucrelrate = (*nucrelratearray)[0];
        ProfileMove(nucrelrate, 0.1, 1, 10, &MultiGeneCoevolModel::NucRatesLogProb,
                    &MultiGeneCoevolModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &MultiGeneCoevolModel::NucRatesLogProb,
                    &MultiGeneCoevolModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &MultiGeneCoevolModel::NucRatesLogProb,
                    &MultiGeneCoevolModel::UpdateNucMatrix, this);

        vector<double> &nucstat = (*nucstatarray)[0];
        ProfileMove(nucstat, 0.1, 1, 10, &MultiGeneCoevolModel::NucRatesLogProb,
                    &MultiGeneCoevolModel::UpdateNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 10, &MultiGeneCoevolModel::NucRatesLogProb,
                    &MultiGeneCoevolModel::UpdateNucMatrix, this);
    }

    // Coevol

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
    // MPI send / receive
    //-------------------

    // Nucleotide Rates

    void MasterSendGlobalNucRates() {
        MasterSendGlobal(nucrelratearray->GetVal(0), nucstatarray->GetVal(0));
    }

    void SlaveReceiveGlobalNucRates() {
        SlaveReceiveGlobal((*nucrelratearray)[0], (*nucstatarray)[0]);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[0], (*nucstatarray)[0]);
        }
    }

    void MasterSendGeneNucRates() {
        MasterSendGeneArray(*nucrelratearray, *nucstatarray);
    }

    void SlaveReceiveGeneNucRates() {
        SlaveReceiveGeneArray(*nucrelratearray, *nucstatarray);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
        }
    }

    void SlaveSendGeneNucRates() {
        SlaveSendGeneArray(*nucrelratearray, *nucstatarray);
    }

    void MasterReceiveGeneNucRates() {
        MasterReceiveGeneArray(*nucrelratearray, *nucstatarray);
    }

    void MasterSendNucRatesHyperParameters() {
        MasterSendGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        MasterSendGlobal(nucstathypercenter, nucstathyperinvconc);
    }

    void SlaveReceiveNucRatesHyperParameters() {
        SlaveReceiveGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        SlaveReceiveGlobal(nucstathypercenter, nucstathyperinvconc);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter, nucrelratehyperinvconc,
                                                          nucstathypercenter, nucstathyperinvconc);
        }
    }

    void SlaveSendNucRatesHyperSuffStat() {
        nucrelratesuffstat.Clear();
        nucrelratearray->AddSuffStat(nucrelratesuffstat);
        SlaveSendAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        nucstatarray->AddSuffStat(nucstatsuffstat);
        SlaveSendAdditive(nucstatsuffstat);
    }

    void MasterReceiveNucRatesHyperSuffStat() {
        nucrelratesuffstat.Clear();
        MasterReceiveAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        MasterReceiveAdditive(nucstatsuffstat);
    }

    void SlaveSendNucPathSuffStat() {
        nucpathsuffstat.Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
        }

        SlaveSendAdditive(nucpathsuffstat);
    }

    void MasterReceiveNucPathSuffStat() {
        nucpathsuffstat.Clear();
        MasterReceiveAdditive(nucpathsuffstat);
    }

    void MasterSendCoevol() {
        MasterSendGlobal(*chronogram, *process);
    }

    void SlaveReceiveCoevol()   {
        SlaveReceiveGlobal(*chronogram, *process);
        for (int gene=0; gene<GetLocalNgene(); gene++) {
            geneprocess[gene]->SetCoevol(*chronogram, *process);
        }
    }

    void SlaveSenddSOmegaSuffStat()   {
        dsompathsuffstatarray->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectdSOmegaPathSuffStat();
            dsompathsuffstatarray->Add(geneprocess[gene]->GetdSOmegaPathSuffStatBranchArray());
        }
        SlaveSendAdditive(*dsompathsuffstatarray);
    }

    void MasterReceivedSOmegaSuffStat()   {
        dsompathsuffstatarray->Clear();
        MasterReceiveAdditive(*dsompathsuffstatarray);
    }

    // log probs

    void SlaveSendLogProbs() {
        GeneLogPrior = 0;
        lnL = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
    }

    void MasterReceiveLogProbs() {
        GeneLogPrior = 0;
        MasterReceiveAdditive(GeneLogPrior);
        lnL = 0;
        MasterReceiveAdditive(lnL);
    }
};
