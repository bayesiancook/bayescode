
#include "SingleOmegaModel.hpp"
#include "Parallel.hpp"
#include "MultiGeneProbModel.hpp"

class MultiGeneSingleOmegaModel : public MultiGeneProbModel {

    private:

	Tree* tree;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	double omegahypermean;
	double omegahyperinvshape;
	IIDGamma* omegaarray;
	GammaSuffStat omegahypersuffstat;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

    NucPathSuffStat nucpathsuffstat;
    std::vector<SingleOmegaModel*> geneprocess;

    double lnL;
    double GeneLogPrior;

    public:

    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneSingleOmegaModel(string datafile, string intreefile, int inmyid, int innprocs) : MultiGeneProbModel(inmyid,innprocs) {

        AllocateAlignments(datafile);
        treefile = intreefile;

        refcodondata = new CodonSequenceAlignment(refdata, true);
        taxonset = refdata->GetTaxonSet();
        Ntaxa = refdata->GetNtaxa();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        if (! myid) {
            std::cerr << "number of taxa : " << Ntaxa << '\n';
            std::cerr << "number of branches : " << Nbranch << '\n';
            std::cerr << "-- Tree and data fit together\n";
        }
    }

    void Allocate() {

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree,1.0,lambda);
        lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

        nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));

        nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omegaarray = new IIDGamma(GetLocalNgene(),1.0,1.0);

        lnL = 0;
        GeneLogPrior = 0;

        if (! GetMyid())    {
            geneprocess.assign(0,(SingleOmegaModel*) 0);
        }
        else    {
            geneprocess.assign(GetLocalNgene(),(SingleOmegaModel*) 0);

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene] = new SingleOmegaModel(GetLocalGeneName(gene),treefile);
            }
        }
    }

    void Unfold()   {

        if (! GetMyid())    {
            MasterSendGlobalBranchLengths();
            MasterSendGlobalNucRates();
            MasterSendOmegaHyperParameters();
            MasterSendOmega();
            MasterReceiveLogProbs();
        }
        else    {

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->Allocate();
            }

            SlaveReceiveGlobalBranchLengths();
            SlaveReceiveGlobalNucRates();
            SlaveReceiveOmegaHyperParameters();
            SlaveReceiveOmega();

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->UpdateMatrices();
                geneprocess[gene]->Unfold();
            }

            SlaveSendLogProbs();
        }
    }

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    //-------------------
    // Traces and Monitors
    //-------------------

    void TraceHeader(ostream& os) const {

        os << "#logprior\tlnL\tlength\t";
        os << "meanomega\t";
        os << "varomega\t";
        os << "omegahypermean\tinvshape\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream& os) const {
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << omegaarray->GetMean() << '\t';
        os << omegaarray->GetVar() << '\t';
        os << omegahypermean << '\t' << omegahyperinvshape << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
		os.flush();
    }

	void Monitor(ostream& os) const {}
	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

    //-------------------
    // Updates
    //-------------------

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
		double total = 0;
		total += BranchLengthsHyperLogPrior();
		total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += OmegaHyperLogPrior();
		total += OmegaLogPrior();
		return total;
    }

	double BranchLengthsHyperLogPrior() const {
		return -lambda / 10;
	}

	double BranchLengthsLogPrior() const {
		return branchlength->GetLogProb();
	}

    double OmegaHyperLogPrior() const {
        double total = 0;
        total -= omegahypermean;
        total -= omegahyperinvshape;
        return total;
    }

    double OmegaLogPrior()   const {
        return omegaarray->GetLogProb();
    }

    double NucRatesLogPrior() const {
        return 0;
    }

    double GetLogLikelihood() const {
        return lnL;
    }

    //-------------------
    // Suff Stat Log Probs
    //-------------------

    // suff stat for moving branch lengths hyperparameter (lambda)
	double BranchLengthsHyperSuffStatLogProb() const {
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    // suff stats for moving nuc rates
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

    // suff stats for moving omega hyper parameters
    double OmegaHyperSuffStatLogProb() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return omegahypersuffstat.GetLogProb(alpha,beta);
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // log prob for moving branch lengths hyperparameter (lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + NucRatesSuffStatLogProb();
    }

    // log prob for moving omega hyperparameters
    double OmegaHyperLogProb() const {
        return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb();
    }

    //-------------------
    // Moves
    //-------------------

    void MasterMove() override {

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            MasterReceiveOmega();
            MoveOmegaHyperParameters();
            MasterSendOmegaHyperParameters();

            MasterReceiveLengthSuffStat();
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();
            MasterSendGlobalBranchLengths();

            MasterReceiveNucPathSuffStat();
            MoveNucRates();
            MasterSendGlobalNucRates();
        }

        MasterReceiveOmega();
        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {

        GeneResampleSub(1.0);

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            GeneCollectPathSuffStat();
            MoveGeneOmegas();
            SlaveSendOmega();
            SlaveReceiveOmegaHyperParameters();

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalBranchLengths();

            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalNucRates();
        }

        SlaveSendOmega();
        SlaveSendLogProbs();
    }

    void GeneResampleSub(double frac)  {

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void GeneCollectPathSuffStat()  {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectPathSuffStat();
        }
    }

    void MoveGeneOmegas()  {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->MoveOmega();
            (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
        }
    }

    void ResampleBranchLengths()    {
		branchlength->GibbsResample(*lengthsuffstatarray);
    }

	void MoveBranchLengthsHyperParameter()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);

        ScalingMove(lambda,1.0,10,&MultiGeneSingleOmegaModel::BranchLengthsHyperLogProb,&MultiGeneSingleOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&MultiGeneSingleOmegaModel::BranchLengthsHyperLogProb,&MultiGeneSingleOmegaModel::NoUpdate,this);

		branchlength->SetScale(lambda);
    }

    void MoveOmegaHyperParameters()  {

		omegahypersuffstat.Clear();
		omegaarray->AddSuffStat(omegahypersuffstat);

        ScalingMove(omegahypermean,1.0,10,&MultiGeneSingleOmegaModel::OmegaHyperLogProb,&MultiGeneSingleOmegaModel::NoUpdate,this);
        ScalingMove(omegahypermean,0.3,10,&MultiGeneSingleOmegaModel::OmegaHyperLogProb,&MultiGeneSingleOmegaModel::NoUpdate,this);
        ScalingMove(omegahyperinvshape,1.0,10,&MultiGeneSingleOmegaModel::OmegaHyperLogProb,&MultiGeneSingleOmegaModel::NoUpdate,this);
        ScalingMove(omegahyperinvshape,0.3,10,&MultiGeneSingleOmegaModel::OmegaHyperLogProb,&MultiGeneSingleOmegaModel::NoUpdate,this);

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    void MoveNucRates()    {

        ProfileMove(nucrelrate,0.1,1,3,&MultiGeneSingleOmegaModel::NucRatesLogProb,&MultiGeneSingleOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.03,3,3,&MultiGeneSingleOmegaModel::NucRatesLogProb,&MultiGeneSingleOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucrelrate,0.01,3,3,&MultiGeneSingleOmegaModel::NucRatesLogProb,&MultiGeneSingleOmegaModel::UpdateNucMatrix,this);

        ProfileMove(nucstat,0.1,1,3,&MultiGeneSingleOmegaModel::NucRatesLogProb,&MultiGeneSingleOmegaModel::UpdateNucMatrix,this);
        ProfileMove(nucstat,0.01,1,3,&MultiGeneSingleOmegaModel::NucRatesLogProb,&MultiGeneSingleOmegaModel::UpdateNucMatrix,this);
    }

    //-------------------
    // MPI send / receive
    //-------------------

    // global branch lengths
    
    void MasterSendGlobalBranchLengths() {
        MasterSendGlobal(*branchlength);
    }

    void SlaveReceiveGlobalBranchLengths()   {
        SlaveReceiveGlobal(*branchlength);
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetBranchLengths(*branchlength);
        }
    }

    // global nuc rates

    void MasterSendGlobalNucRates()   {
        MasterSendGlobal(nucrelrate,nucstat);
    }

    void SlaveReceiveGlobalNucRates()   {

        SlaveReceiveGlobal(nucrelrate,nucstat);
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
        }
    }

    // omega (and hyperparameters)

    void SlaveSendOmega()   {
        SlaveSendGeneArray(*omegaarray);
    }

    void MasterReceiveOmega()    {
        MasterReceiveGeneArray(*omegaarray);
    }

    void MasterSendOmega()  {
        MasterSendGeneArray(*omegaarray);
    }
    
    void SlaveReceiveOmega()    {
        SlaveReceiveGeneArray(*omegaarray);
    }

    // omega hyperparameters

    void MasterSendOmegaHyperParameters()   {
        MasterSendGlobal(omegahypermean,omegahyperinvshape);
    }

    void SlaveReceiveOmegaHyperParameters() {
        SlaveReceiveGlobal(omegahypermean,omegahyperinvshape);
        for (int gene=0; gene<GetLocalNgene(); gene++)    {
            geneprocess[gene]->SetOmegaHyperParameters(omegahypermean,omegahyperinvshape);
        }
    }

    // branch length suff stat

    void SlaveSendLengthSuffStat()  {
        lengthsuffstatarray->Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectLengthSuffStat();
            lengthsuffstatarray->Add(*geneprocess[gene]->GetLengthSuffStatArray());
        }
        SlaveSendAdditive(*lengthsuffstatarray);
    }

    void MasterReceiveLengthSuffStat()  {

        lengthsuffstatarray->Clear();
        MasterReceiveAdditive(*lengthsuffstatarray);
    }

    // nuc rate suff stat

    void SlaveSendNucPathSuffStat()  {
        nucpathsuffstat.Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
        }
        SlaveSendAdditive(nucpathsuffstat);
    }

    void MasterReceiveNucPathSuffStat()  {
        nucpathsuffstat.Clear();
        MasterReceiveAdditive(nucpathsuffstat);
    }

    // log probs

    void SlaveSendLogProbs()   {

        GeneLogPrior = 0;
        lnL = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
    }

    void MasterReceiveLogProbs()    {

        GeneLogPrior = 0;
        MasterReceiveAdditive(GeneLogPrior);
        lnL = 0;
        MasterReceiveAdditive(lnL);
    }
};

