#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GeneBranchCodonArray.hpp"
#include "RecursiveNewick.hpp"

class GeneBranchNucModel : public ProbModel {

  private:

    int integrated_move;
    int syn_devmode;
    int nuc_devmode;

    const TaxonSet* taxonset;
    Tree *tree;
    int Ntaxa;
    int Nbranch;
    int Ngene;

    // branch lengths

    GeneBranchGammaEffects* syn_model;

    // nucleotide models

    GeneBranchGammaEffects* rho_AC_model;
    GeneBranchGammaEffects* rho_AG_model;
    GeneBranchGammaEffects* rho_CG_model;
    GeneBranchGammaEffects* rho_CA_model;
    GeneBranchGammaEffects* rho_CT_model;

    NucRatesGeneBranchArray* nucrates;
    NucMatrixGeneBranchArray* nucmat;

    PathSuffStatGeneBranchArray *pathss;
    NucPathSuffStatGeneBranchArray *nucss;
    LengthPathSuffStatGeneBranchArray *synss;

    vector<string> gene_names;

  public:

    GeneBranchNucModel(string datafile, string treefile, string taxonfile, int insyn_devmode, int innuc_devmode) {

        syn_devmode = insyn_devmode;
        nuc_devmode = innuc_devmode;

        integrated_move = 1;

        taxonset = new TaxonSet(taxonfile);

        tree = new Tree(treefile);
        tree->RegisterWith(taxonset);
        tree->SetIndices();
        Ntaxa = tree->GetSize();
        Nbranch = tree->GetNbranch();

        ifstream is(datafile.c_str());
        is >> Ngene;
        cerr << "number of genes: " << Ngene << '\n';
        gene_names.assign(Ngene,"");

        cerr << "allocate\n";
        Allocate();

        cerr << "read suff stat\n";
        for (int gene=0; gene<Ngene; gene++)    {
            is >> gene_names[gene];
            cerr << gene_names[gene] << '\n';
            RelativePathSuffStatNodeArray tmp(*tree, Nnuc);
            is >> tmp;
            pathss->Get(gene, tmp);
        }

        /*
        vector<double> bl = pathss->GetEmpiricalBranchLengths();
        double tot = 0;
        for (int j=0; j<Nbranch; j++)   {
            tot += bl[j];
            cerr << j << '\t' << bl[j] << '\n';
        }
        cerr << "total : " << tot << '\n';
        ofstream os("tmp");
        ToNewick(os, *tree, bl, bl);
        os.close();
        */
    }

    void Allocate() {

        syn_model = new GeneBranchGammaEffects(Ngene, Nbranch, syn_devmode, 1, 0);

        rho_AC_model = new GeneBranchGammaEffects(Ngene, Nbranch, nuc_devmode, 0, 1);
        rho_AG_model = new GeneBranchGammaEffects(Ngene, Nbranch, nuc_devmode, 0, 1);
        rho_CA_model = new GeneBranchGammaEffects(Ngene, Nbranch, nuc_devmode, 0, 1);
        rho_CG_model = new GeneBranchGammaEffects(Ngene, Nbranch, nuc_devmode, 0, 1);
        rho_CT_model = new GeneBranchGammaEffects(Ngene, Nbranch, nuc_devmode, 0, 1);
            
        nucrates = new NucRatesGeneBranchArray(Ngene, Nbranch, 
                *rho_AC_model, *rho_AG_model, *rho_CA_model, *rho_CG_model, *rho_CT_model);

        nucmat = new NucMatrixGeneBranchArray(Ngene, Nbranch, *nucrates, false);
        pathss = new PathSuffStatGeneBranchArray(Ngene, Nbranch, Nnuc);
        nucss = new NucPathSuffStatGeneBranchArray(Ngene, Nbranch);
        synss = new LengthPathSuffStatGeneBranchArray(Ngene, Nbranch);
    }

    int GetNtaxa() const    {
        return Ntaxa;
    }

    int GetNbranch() const  {
        return Nbranch;
    }

    int GetNgene() const    {
        return Ngene;
    }

    string GetGeneName(int gene) const  {
        return gene_names[gene];
    }

    const Tree* GetTree() const {
        return tree;
    }

    const GeneBranchGammaEffects* GetSynModel() const   {
        return syn_model;
    }

    void TraceHeader(ostream &os) const override {
        os << "logprior\tlogl\tlength";
        syn_model->TraceHeader(os, "syn");
        rho_AC_model->TraceHeader(os, "AC");
        rho_AG_model->TraceHeader(os, "AG");
        rho_CA_model->TraceHeader(os, "CA");
        rho_CG_model->TraceHeader(os, "CG");
        rho_CT_model->TraceHeader(os, "CT");
        os << '\n';
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t' << GetLogLikelihood();
        os << '\t' << syn_model->GetMeanBranchTotal() * nucmat->GetMeanRate();
        syn_model->Trace(os);
        rho_AC_model->Trace(os);
        rho_AG_model->Trace(os);
        rho_CA_model->Trace(os);
        rho_CG_model->Trace(os);
        rho_CT_model->Trace(os);
        os << '\n';
    }

    void Monitor(ostream &os) const override {
    }

    void ToStream(ostream &os) const override {
        syn_model->ToStream(os);
        rho_AC_model->ToStream(os);
        rho_AG_model->ToStream(os);
        rho_CA_model->ToStream(os);
        rho_CG_model->ToStream(os);
        rho_CT_model->ToStream(os);
        os << '\n';
    }

    void FromStream(istream &is) override {
        syn_model->FromStream(is);
        rho_AC_model->FromStream(is);
        rho_AG_model->FromStream(is);
        rho_CA_model->FromStream(is);
        rho_CG_model->FromStream(is);
        rho_CT_model->FromStream(is);
    }

    void Update() override {
        syn_model->Update();
        rho_AC_model->Update();
        rho_AG_model->Update();
        rho_CA_model->Update();
        rho_CG_model->Update();
        rho_CT_model->Update();
        nucrates->Update();
        nucmat->UpdateMatrices();
    }

    void NoUpdate() {}

    // Log Priors

    double GetLogPrior() const {
        double total = 0;
        total += syn_model->GetLogPrior();
        total += rho_AC_model->GetLogPrior();
        total += rho_AG_model->GetLogPrior();
        total += rho_CA_model->GetLogPrior();
        total += rho_CG_model->GetLogPrior();
        total += rho_CT_model->GetLogPrior();
        return total;
    }

    // likelihood
    
    double GetLogLikelihood() const {
        double tot = 0;
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                tot += GetLogLikelihood(i,j);
            }
        }
        return tot;
    }

    double GetLogLikelihood(int gene, int branch) const {
        return pathss->GetVal(gene,branch).GetLogProb(
                nucmat->GetVal(gene,branch), syn_model->GetVal(gene,branch));
    }

    double GetLogProb() const override  {
        return GetLogPrior() + GetLogLikelihood();
    }

    double Move() override  {
        for (int rep=0; rep<10; rep++)  {
            CollectNucPathSuffStat();
            MoveNuc(3,1);
            CollectSynPathSuffStat();
            MoveSyn(3,1);
        }
        return 1.0;
    }

    void CollectNucPathSuffStat()   {
        nucss->Clear();
        nucss->AddSuffStat(*pathss, *syn_model);
    }

    void CollectSynPathSuffStat()  {
        synss->Clear();
        synss->AddSuffStat(*pathss, *nucmat);
    }

    double MoveNuc(int nrep, int nsmallrep) {

        auto get_AC_ss = [this] (int gene, int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(gene, branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(0,1) + suffstat.GetPairCount(3,2),
                    suffstat.GetPairBeta(0,1) + suffstat.GetPairBeta(3,2));
        };

        auto get_AG_ss = [this] (int gene, int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(gene, branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(0,2) + suffstat.GetPairCount(3,1),
                    suffstat.GetPairBeta(0,2) + suffstat.GetPairBeta(3,1));
        };

        auto get_CA_ss = [this] (int gene, int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(gene, branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(1,0) + suffstat.GetPairCount(2,3),
                    suffstat.GetPairBeta(1,0) + suffstat.GetPairBeta(2,3));
        };

        auto get_CG_ss = [this] (int gene, int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(gene, branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(1,2) + suffstat.GetPairCount(2,1),
                    suffstat.GetPairBeta(1,2) + suffstat.GetPairBeta(2,1));
        };

        auto get_CT_ss = [this] (int gene, int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(gene, branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(1,3) + suffstat.GetPairCount(2,0),
                    suffstat.GetPairBeta(1,3) + suffstat.GetPairBeta(2,0));
        };

        for (int rep=0; rep<nrep; rep++) {
            if (integrated_move)    {
                rho_AC_model->IntegratedMove(1.0, nsmallrep, get_AC_ss);
                rho_AG_model->IntegratedMove(1.0, nsmallrep, get_AG_ss);
                rho_CA_model->IntegratedMove(1.0, nsmallrep, get_CA_ss);
                rho_CG_model->IntegratedMove(1.0, nsmallrep, get_CG_ss);
                rho_CT_model->IntegratedMove(1.0, nsmallrep, get_CT_ss);
            }
            else    {
                rho_AC_model->NonIntegratedMove(1.0, nsmallrep, get_AC_ss);
                rho_AG_model->NonIntegratedMove(1.0, nsmallrep, get_AG_ss);
                rho_CA_model->NonIntegratedMove(1.0, nsmallrep, get_CA_ss);
                rho_CG_model->NonIntegratedMove(1.0, nsmallrep, get_CG_ss);
                rho_CT_model->NonIntegratedMove(1.0, nsmallrep, get_CT_ss);
            }
            rho_AC_model->MoveHyper(1.0, 1);
            rho_AG_model->MoveHyper(1.0, 1);
            rho_CA_model->MoveHyper(1.0, 1);
            rho_CG_model->MoveHyper(1.0, 1);
            rho_CT_model->MoveHyper(1.0, 1);
        }
        Update();
        return 1.0;
    }

    double MoveSyn(int nrep, int nsmallrep) {

        auto get_syn_ss = [this] (int gene, int branch) {
            return this->synss->GetVal(gene,branch);
        };

        for (int rep=0; rep<nrep; rep++) {
            if (integrated_move)    {
                syn_model->IntegratedMove(1.0, nsmallrep, get_syn_ss);
            }
            else    {
                syn_model->NonIntegratedMove(1.0, nsmallrep, get_syn_ss);
            }
            syn_model->MoveHyper(1.0, 1);
        }
        Update();
        return 1.0;
    }

    void AddNucStats(double& meanAC, double& meanAG, double& meanCA, double& meanCG, double& meanCT,
            double& geneAC, double& geneAG, double& geneCA, double& geneCG, double& geneCT,
            double& branchAC, double& branchAG, double& branchCA, double& branchCG, double& branchCT,
            double& devAC, double& devAG, double& devCA, double& devCG, double& devCT)  {

        rho_AC_model->AddStats(meanAC, geneAC, branchAC, devAC);
        rho_AG_model->AddStats(meanAG, geneAG, branchAG, devAG);
        rho_CA_model->AddStats(meanCA, geneCA, branchCA, devCA);
        rho_CG_model->AddStats(meanCG, geneCG, branchCG, devCG);
        rho_CT_model->AddStats(meanCT, geneCT, branchCT, devCT);
    }

    void AddEmpiricalBranchLengthsFromNucPathSuffStat(vector<double>& into)    {
        CollectNucPathSuffStat();
        vector<double> tmp = nucss->GetEmpiricalBranchLengths();
        for (int j=0; j<Nbranch; j++)   {
            into[j] += tmp[j];
        }
    }
};

