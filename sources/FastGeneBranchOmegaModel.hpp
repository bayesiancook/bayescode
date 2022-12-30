
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "dSOmegaPathSuffStatGeneBranchArray.hpp"
#include "RecursiveNewick.hpp"
#include "GeneBranchGammaEffects.hpp"

class FastGeneBranchOmegaModel : public ProbModel {

  private:

    int integrated_move;
    int syn_devmode;
    int om_devmode;

    Tree *tree;
    // const TaxonSet *taxonset;
    int Ntaxa;
    int Nbranch;
    int Ngene;

    GeneBranchGammaEffects* syn_model;
    GeneBranchGammaEffects* om_model;

    dSOmegaPathSuffStatGeneBranchArray *dsomss;
    vector<string> gene_names;

  public:

    FastGeneBranchOmegaModel(string datafile, string treefile, int insyn_devmode, int inom_devmode) {

        syn_devmode = insyn_devmode;
        om_devmode = inom_devmode;
        integrated_move = 1;

        // taxonset = new TaxonSet(taxonfile);

        tree = new Tree(treefile);
        // tree->RegisterWith(taxonset);
        tree->SetIndices();
        Ntaxa = tree->GetSize();
        Nbranch = tree->GetNbranch();

        ifstream is(datafile.c_str());
        is >> Ngene;
        gene_names.assign(Ngene,"");

        syn_model = new GeneBranchGammaEffects(Ngene, Nbranch, syn_devmode, 1, 0);
        om_model = new GeneBranchGammaEffects(Ngene, Nbranch, syn_devmode, 0, 1);

        dsomss = new dSOmegaPathSuffStatGeneBranchArray(*tree, Ngene);

        // read suff stats
        for (int gene=0; gene<Ngene; gene++)    {
            is >> gene_names[gene];

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

            (*dsomss)[gene].Add(treedscount, treedsbeta, treedncount, treednbeta);
        }
    }

    void SetIntegratedMove(int in)  {
        integrated_move = in;
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

    const GeneBranchGammaEffects* GetOmegaModel() const {
        return om_model;
    }

    void TraceHeader(ostream &os) const override {
        os << "logprior\tlogl";
        syn_model->TraceHeader(os, "syn");
        om_model->TraceHeader(os, "om");
        os << '\n';
    }


    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t' << GetLogLikelihood();
        syn_model->Trace(os);
        om_model->Trace(os);
        os << '\n';
    }

    void Monitor(ostream &os) const override {
    }

    void ToStream(ostream &os) const override {
        syn_model->ToStream(os);
        om_model->ToStream(os);
        os << '\n';
    }

    void FromStream(istream &is) override {
        syn_model->FromStream(is);
        om_model->FromStream(is);
    }

    void Update() override {
        syn_model->Update();
        om_model->Update();
    }

    void NoUpdate() {}

    // Log Priors

    double GetLogPrior() const {
        double total = 0;
        total += syn_model->GetLogPrior();
        total += om_model->GetLogPrior();
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

        const dSOmegaPathSuffStat &suffstat = dsomss->GetVal(gene).GetVal(branch);

        double syncount = suffstat.GetSynCount();
        double synbeta = suffstat.GetSynBeta();
        double nonsyncount = suffstat.GetNonSynCount();
        double nonsynbeta = suffstat.GetNonSynBeta();

        double syn = syn_model->GetVal(gene, branch);
        double om = om_model->GetVal(gene, branch);

        return syncount * log(syn) - synbeta * syn 
            + nonsyncount * log(syn*om) - nonsynbeta * (syn*om);

    }

    double GetLogProb() const override  {
        return GetLogPrior() + GetLogLikelihood();
    }

    double Move() override  {
        int nrep = 10;
        int nsmallrep = 3;

        auto get_syn_ss = [this] (int gene, int branch) {
            const dSOmegaPathSuffStat &suffstat = this->dsomss->GetVal(gene).GetVal(branch);
            double om = this->om_model->GetVal(gene,branch);
            return MeanPoissonSuffStat(
                    suffstat.GetSynCount() + suffstat.GetNonSynCount(),
                    suffstat.GetSynBeta() + om*suffstat.GetNonSynBeta());
        };

        auto get_om_ss = [this] (int gene, int branch)  {
            const dSOmegaPathSuffStat &suffstat = this->dsomss->GetVal(gene).GetVal(branch);
            double syn = this->syn_model->GetVal(gene,branch);
            return MeanPoissonSuffStat(
                    suffstat.GetNonSynCount(),
                    syn*suffstat.GetNonSynBeta());
        };

        for (int rep=0; rep<nrep; rep++) {
            if (integrated_move)    {
                syn_model->IntegratedMove(1.0, nsmallrep, get_syn_ss);
            }
            else    {
                syn_model->NonIntegratedMove(1.0, nsmallrep, get_syn_ss);
            }
            syn_model->MoveHyper(1.0, 1);
            if (integrated_move)    {
                om_model->IntegratedMove(1.0, nsmallrep, get_om_ss);
            }
            else    {
                om_model->NonIntegratedMove(1.0, nsmallrep, get_om_ss);
            }
            om_model->MoveHyper(1.0, 1);
        }
        Update();
        return 1.0;
    }

    void AddSynDevPostProbsTo(vector<vector<double>>& array) const   {
        auto get_syn_ss = [this] (int gene, int branch) {
            const dSOmegaPathSuffStat &suffstat = this->dsomss->GetVal(gene).GetVal(branch);
            double om = this->om_model->GetVal(gene,branch);
            return MeanPoissonSuffStat(
                    suffstat.GetSynCount() + suffstat.GetNonSynCount(),
                    suffstat.GetSynBeta() + om*suffstat.GetNonSynBeta());
        };
        syn_model->AddDevPostProbsTo(array, get_syn_ss);
    }
    
    void AddOmegaDevPostProbsTo(vector<vector<double>>& array) const   {
        auto get_om_ss = [this] (int gene, int branch)  {
            const dSOmegaPathSuffStat &suffstat = this->dsomss->GetVal(gene).GetVal(branch);
            double syn = this->syn_model->GetVal(gene,branch);
            return MeanPoissonSuffStat(
                    suffstat.GetNonSynCount(),
                    syn*suffstat.GetNonSynBeta());
        };
        om_model->AddDevPostProbsTo(array, get_om_ss);
    }
};

