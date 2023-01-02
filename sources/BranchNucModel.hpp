#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GeneBranchCodonArray.hpp"
#include "RecursiveNewick.hpp"

class BranchNucModel : public ProbModel {

  private:

    int integrated_move;
    int syn_devmode;

    const TaxonSet* taxonset;
    Tree *tree;
    int Ntaxa;
    int Nbranch;
    int Ngene;

    // branch lengths

    GeneBranchGammaEffects* syn_model;

    // nucleotide models
    double meanAC, varAC;
    double meanAG, varAG;
    double meanCA, varCA;
    double meanCG, varCG;
    double meanCT, varCT;

    IIDGamma* rho_AC_model;
    IIDGamma* rho_AG_model;
    IIDGamma* rho_CG_model;
    IIDGamma* rho_CA_model;
    IIDGamma* rho_CT_model;

    NucRatesBranchArray* nucrates;
    NucMatrixBranchArray* nucmat;
    BidimRowHomogeneousSelector<StrandSymmetricIrreversibleSubMatrix>* nucmatselector;

    PathSuffStatGeneBranchArray *pathss;
    LengthPathSuffStatGeneBranchArray *synss;

    NucPathSuffStatBranchArray *nucss;
    NucRateSuffStatBranchArray* ss_AC;
    NucRateSuffStatBranchArray* ss_AG;
    NucRateSuffStatBranchArray* ss_CA;
    NucRateSuffStatBranchArray* ss_CG;
    NucRateSuffStatBranchArray* ss_CT;

    vector<string> gene_names;

  public:

    BranchNucModel(string datafile, string treefile, string taxonfile, int insyn_devmode) {

        syn_devmode = insyn_devmode;

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
    }

    void Allocate() {

        syn_model = new GeneBranchGammaEffects(Ngene, Nbranch, syn_devmode, 1, 0);

        meanAC = varAC = 1.0;
        meanAG = varAG = 1.0;
        meanCA = varCA = 1.0;
        meanCG = varCG = 1.0;
        meanCT = varCT = 1.0;

        rho_AC_model = new IIDGamma(Nbranch, 1.0, 1.0);
        rho_AG_model = new IIDGamma(Nbranch, 1.0, 1.0);
        rho_CA_model = new IIDGamma(Nbranch, 1.0, 1.0);
        rho_CG_model = new IIDGamma(Nbranch, 1.0, 1.0);
        rho_CT_model = new IIDGamma(Nbranch, 1.0, 1.0);
            
        nucrates = new NucRatesBranchArray(Nbranch, 
                *rho_AC_model, *rho_AG_model, *rho_CA_model, *rho_CG_model, *rho_CT_model);

        nucmat = new NucMatrixBranchArray(Nbranch, *nucrates, false);
        nucmatselector = new BidimRowHomogeneousSelector<StrandSymmetricIrreversibleSubMatrix>(Ngene, *nucmat);

        pathss = new PathSuffStatGeneBranchArray(Ngene, Nbranch, Nnuc);
        synss = new LengthPathSuffStatGeneBranchArray(Ngene, Nbranch);
        nucss = new NucPathSuffStatBranchArray(Nbranch);

        ss_AC = new NucRateSuffStatBranchArray(Nbranch);
        ss_AG = new NucRateSuffStatBranchArray(Nbranch);
        ss_CA = new NucRateSuffStatBranchArray(Nbranch);
        ss_CG = new NucRateSuffStatBranchArray(Nbranch);
        ss_CT = new NucRateSuffStatBranchArray(Nbranch);
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
        os << "logprior\tlogl";
        syn_model->TraceHeader(os, "syn");
        os << "\tmeanAC\tvarAC";
        os << "\tmeanAG\tvarAG";
        os << "\tmeanCA\tvarCA";
        os << "\tmeanCG\tvarCG";
        os << "\tmeanCT\tvarCT";
        os << '\n';
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t' << GetLogLikelihood();
        syn_model->Trace(os);
        os << '\t' << meanAC << '\t' << varAC;
        os << '\t' << meanAG << '\t' << varAG;
        os << '\t' << meanCA << '\t' << varCA;
        os << '\t' << meanCG << '\t' << varCG;
        os << '\t' << meanCT << '\t' << varCT;
        os << '\n';
    }

    void Monitor(ostream &os) const override {
    }

    void ToStream(ostream &os) const override {
        syn_model->ToStream(os);
        os << '\t' << meanAC << '\t' << varAC;
        os << '\t' << meanAG << '\t' << varAG;
        os << '\t' << meanCA << '\t' << varCA;
        os << '\t' << meanCG << '\t' << varCG;
        os << '\t' << meanCT << '\t' << varCT;
        os << '\t' << *rho_AC_model;
        os << '\t' << *rho_AG_model;
        os << '\t' << *rho_CA_model;
        os << '\t' << *rho_CG_model;
        os << '\t' << *rho_CT_model;
        os << '\n';
    }

    void FromStream(istream &is) override {
        syn_model->FromStream(is);
        is >> meanAC >> varAC;
        is >> meanAG >> varAG;
        is >> meanCA >> varCA;
        is >> meanCG >> varCG;
        is >> meanCT >> varCT;
        is >> *rho_AC_model;
        is >> *rho_AG_model;
        is >> *rho_CA_model;
        is >> *rho_CG_model;
        is >> *rho_CT_model;
    }

    void Update(IIDGamma* array, double mean, double invshape)    {
        double alpha = 1.0 / invshape;
        double beta = alpha / mean;
        array->SetShape(alpha);
        array->SetScale(beta);
    }

    void Update() override {
        syn_model->Update();

        Update(rho_AC_model, meanAC, varAC);
        Update(rho_AG_model, meanAG, varAG);
        Update(rho_CA_model, meanCA, varCA);
        Update(rho_CG_model, meanCG, varCG);
        Update(rho_CT_model, meanCT, varCT);

        nucrates->Update();
        nucmat->UpdateMatrices();
    }

    void NoUpdate() {}

    // Log Priors

    double GetLogPrior() const {
        double total = 0;
        total += syn_model->GetLogPrior();
        total += GetNucHyperLogPrior();
        total += GetNucLogPrior();
        return total;
    }

    double GetNucHyperLogPrior() const   {
        double total = 0;
        total -= meanAC + varAC;
        total -= meanAG + varAG;
        total -= meanCA + varCA;
        total -= meanCG + varCG;
        total -= meanCT + varCT;
        return total;
    }

    double GetNucLogPrior() const   {
        double total = 0;
        total += rho_AC_model->GetLogProb();
        total += rho_AG_model->GetLogProb();
        total += rho_CG_model->GetLogProb();
        total += rho_CA_model->GetLogProb();
        total += rho_CT_model->GetLogProb();
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
                nucmat->GetVal(branch), syn_model->GetVal(gene,branch));
    }

    double GetLogProb() const override  {
        return GetLogPrior() + GetLogLikelihood();
    }

    double Move() override  {
        for (int rep=0; rep<10; rep++)  {
            CollectNucPathSuffStat();
            CollectNucRateSuffStat();
            MoveNuc(3);
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

    void CollectNucRateSuffStat()   {

        auto get_AC_ss = [this] (int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(0,1) + suffstat.GetPairCount(3,2),
                    suffstat.GetPairBeta(0,1) + suffstat.GetPairBeta(3,2));
        };

        ss_AC->Clear();
        ss_AC->AddSuffStat(get_AC_ss);

        auto get_AG_ss = [this] (int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(0,2) + suffstat.GetPairCount(3,1),
                    suffstat.GetPairBeta(0,2) + suffstat.GetPairBeta(3,1));
        };

        ss_AG->Clear();
        ss_AG->AddSuffStat(get_AG_ss);

        auto get_CA_ss = [this] (int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(1,0) + suffstat.GetPairCount(2,3),
                    suffstat.GetPairBeta(1,0) + suffstat.GetPairBeta(2,3));
        };

        ss_CA->Clear();
        ss_CA->AddSuffStat(get_CA_ss);

        auto get_CG_ss = [this] (int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(1,2) + suffstat.GetPairCount(2,1),
                    suffstat.GetPairBeta(1,2) + suffstat.GetPairBeta(2,1));
        };

        ss_CG->Clear();
        ss_CG->AddSuffStat(get_CG_ss);

        auto get_CT_ss = [this] (int branch) {
            const MeanNucPathSuffStat &suffstat = this->nucss->GetVal(branch);
            return MeanPoissonSuffStat(
                    suffstat.GetPairCount(1,3) + suffstat.GetPairCount(2,0),
                    suffstat.GetPairBeta(1,3) + suffstat.GetPairBeta(2,0));
        };

        ss_CT->Clear();
        ss_CT->AddSuffStat(get_CT_ss);
    }

    double MoveNuc(int nrep)    {

        auto no_update = [] () {};

        auto logprob_AC = [this] () {
            return this->ss_AC->GetMarginalLogProbMeanInvShape(this->meanAC, this->varAC);
        };

        auto logprob_AG = [this] () {
            return this->ss_AG->GetMarginalLogProbMeanInvShape(this->meanAG, this->varAG);
        };

        auto logprob_CA = [this] () {
            return this->ss_CA->GetMarginalLogProbMeanInvShape(this->meanCA, this->varCA);
        };

        auto logprob_CG = [this] () {
            return this->ss_CG->GetMarginalLogProbMeanInvShape(this->meanCG, this->varCG);
        };

        auto logprob_CT = [this] () {
            return this->ss_CT->GetMarginalLogProbMeanInvShape(this->meanCT, this->varCT);
        };

        ScalingMove(meanAC, 0.3, nrep, logprob_AC, no_update);
        ScalingMove(varAC, 0.3, nrep, logprob_AC, no_update);
        Update(rho_AC_model, meanAC, varAC);
        rho_AC_model->GibbsResample(*ss_AC);

        ScalingMove(meanAG, 0.3, nrep, logprob_AG, no_update);
        ScalingMove(varAG, 0.3, nrep, logprob_AG, no_update);
        Update(rho_AG_model, meanAC, varAG);
        rho_AG_model->GibbsResample(*ss_AG);

        ScalingMove(meanCA, 0.3, nrep, logprob_CA, no_update);
        ScalingMove(varCA, 0.3, nrep, logprob_CA, no_update);
        Update(rho_CA_model, meanAC, varCA);
        rho_CA_model->GibbsResample(*ss_CA);

        ScalingMove(meanCG, 0.3, nrep, logprob_CG, no_update);
        ScalingMove(varCG, 0.3, nrep, logprob_CG, no_update);
        Update(rho_CG_model, meanAC, varCG);
        rho_CG_model->GibbsResample(*ss_CG);

        ScalingMove(meanCT, 0.3, nrep, logprob_CT, no_update);
        ScalingMove(varCT, 0.3, nrep, logprob_CT, no_update);
        Update(rho_CT_model, meanAC, varCT);
        rho_CT_model->GibbsResample(*ss_CT);

        Update();
        return 1.0;
    }

    template <class LogProbF, class UpdateF>
    double ScalingMove(double &x, double tuning, int nrep, LogProbF logprobf, UpdateF updatef)  {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -logprobf();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            x *= e;
            updatef();
            deltalogprob += logprobf();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x /= e;
                updatef();
            }
            ntot++;
        }
        return nacc / ntot;
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

    void AddNucStatsTo(double& mAC, double& mAG, double& mCA, double& mCG, double& mCT,
		    double& vAC, double& vAG, double& vCA, double& vCG, double& vCT)	{
	    mAC += meanAC;
	    mAG += meanAG;
	    mCA += meanCA;
	    mCG += meanCG;
	    mCT += meanCT;
	    vAC += varAC;
	    vAG += varAG;
	    vCA += varCA;
	    vCG += varCG;
	    vCT += varCT;
    }

    void AddNucBranchRatesTo(vector<vector<double>>& array) const   {
        for (int j=0; j<Nbranch; j++)   {
            array[j][0] += rho_AC_model->GetVal(j);
            array[j][1] += rho_AG_model->GetVal(j);
            array[j][2] += rho_CA_model->GetVal(j);
            array[j][3] += rho_CG_model->GetVal(j);
            array[j][4] += rho_CT_model->GetVal(j);
        }
    }
};

