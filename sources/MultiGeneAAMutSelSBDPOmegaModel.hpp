
#include "AAMutSelSBDPOmegaModel.hpp"
#include "Parallel.hpp"
#include "MultiGeneProbModel.hpp"

class MultiGeneAAMutSelSBDPOmegaModel : public MultiGeneProbModel {

    private:

	Tree* tree;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
    // gathering gene-specific omegas
    // and implementing the hyperprior over them
	double omegahypermean;
	double omegahyperinvshape;
	IIDGamma* omegaarray;
    // suff stats for gamma-distributed omega's across genes
	GammaSuffStat omegahypersuffstat;

    // branch lengths shared across genes
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
	GammaSuffStat lambdasuffstat;

    std::vector<AAMutSelSBDPOmegaModel*> geneprocess;

    double lnL;
    double GeneLogPrior;
    double MeanNcluster;
    double MeanAAEntropy;

    int Ncat;

    public:

    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneAAMutSelSBDPOmegaModel(string datafile, string intreefile, int inNcat, int inmyid, int innprocs) : MultiGeneProbModel(inmyid,innprocs) {

        AllocateAlignments(datafile);
        treefile = intreefile;
        Ncat = inNcat;

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

        // lambda and branchlengths  will be resampled at the global level
        lambda = 10;
        branchlength = new BranchIIDGamma(*tree,1.0,lambda);

        // will collect suff stats from genes
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        // these 2 hyperparameters will be resampled at the global level
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;

        // these are just copies from gene-specific omega's
		omegaarray = new IIDGamma(GetLocalNgene(),1.0,1.0);

        lnL = 0;
        GeneLogPrior = 0;

        if (! GetMyid())    {
            geneprocess.assign(0,(AAMutSelSBDPOmegaModel*) 0);
        }
        else    {
            geneprocess.assign(GetLocalNgene(),(AAMutSelSBDPOmegaModel*) 0);

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene] = new AAMutSelSBDPOmegaModel(GetLocalGeneName(gene),treefile,Ncat);
                geneprocess[gene]->SetFixBL(1);
            }
        }
    }

    void Unfold()   {

        if (! GetMyid())    {
            MasterSendGlobalBranchLengths();
            MasterSendOmegaHyperParameters();
            MasterSendOmega();
            MasterReceiveLogProbs();
            MasterReceiveSummaryStats();
        }
        else    {

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->Allocate();
            }

            SlaveReceiveGlobalBranchLengths();
            SlaveReceiveOmegaHyperParameters();
            SlaveReceiveOmega();

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->UpdateMatrices();
                geneprocess[gene]->Unfold();
            }

            SlaveSendLogProbs();
            SlaveSendSummaryStats();
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
        os << "meanncluster\t";
        os << "meanstatent\n";
    }

    void Trace(ostream& os) const {
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << omegaarray->GetMean() << '\t';
        os << omegaarray->GetVar() << '\t';
        os << MeanNcluster << '\t';
        os << MeanAAEntropy << '\n';

		os.flush();
    }

    void TraceGeneOmegas(ostream& os) const {
        os << *omegaarray;
        os.flush();
    }

	void Monitor(ostream& os) const {}
	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

    //-------------------
    // Updates
    //-------------------

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
		double total = 0;
		total += BranchLengthsHyperLogPrior();
		total += BranchLengthsLogPrior();
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

            MasterReceiveLengthSuffStat();
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();
            MasterSendGlobalBranchLengths();

            MasterReceiveOmega();
            MoveOmegaHyperParameters();
            MasterSendOmegaHyperParameters();
        }

        MasterReceiveLogProbs();
        MasterReceiveSummaryStats();
    }

    // slave move
    void SlaveMove() override {

        GeneResampleSub(1.0);

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalBranchLengths();

            GeneMoveParameters(1);

            SlaveSendOmega();
            SlaveReceiveOmegaHyperParameters();
        }

        SlaveSendLogProbs();
        SlaveSendSummaryStats();
    }

    void GeneResampleSub(double frac)  {

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void GeneMoveParameters(int nrep)  {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->MoveParameters(1);
        }
    }

    void ResampleBranchLengths()    {
		branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

	void MoveBranchLengthsHyperParameter()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);

        ScalingMove(lambda,1.0,10,&MultiGeneAAMutSelSBDPOmegaModel::BranchLengthsHyperLogProb,&MultiGeneAAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&MultiGeneAAMutSelSBDPOmegaModel::BranchLengthsHyperLogProb,&MultiGeneAAMutSelSBDPOmegaModel::NoUpdate,this);

		branchlength->SetScale(lambda);
    }

    void MoveOmegaHyperParameters()  {

		omegahypersuffstat.Clear();
		omegaarray->AddSuffStat(omegahypersuffstat);

        ScalingMove(omegahypermean,1.0,10,&MultiGeneAAMutSelSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(omegahypermean,0.3,10,&MultiGeneAAMutSelSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(omegahyperinvshape,1.0,10,&MultiGeneAAMutSelSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelSBDPOmegaModel::NoUpdate,this);
        ScalingMove(omegahyperinvshape,0.3,10,&MultiGeneAAMutSelSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelSBDPOmegaModel::NoUpdate,this);

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
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

    // omega (and hyperparameters)

    void SlaveSendOmega()   {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
        }
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
        lengthpathsuffstatarray->Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectLengthSuffStat();
            lengthpathsuffstatarray->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
        }
        SlaveSendAdditive(*lengthpathsuffstatarray);
    }

    void MasterReceiveLengthSuffStat()  {
        lengthpathsuffstatarray->Clear();
        MasterReceiveAdditive(*lengthpathsuffstatarray);
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

    // mean number of clusters across genes (dirichlet processes)

    void SlaveSendSummaryStats()   {

        MeanNcluster = 0;
        MeanAAEntropy = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            MeanNcluster += geneprocess[gene]->GetNcluster();
            MeanAAEntropy += geneprocess[gene]->GetMeanAAEntropy();
        }
        SlaveSendAdditive(MeanNcluster);
        SlaveSendAdditive(MeanAAEntropy);
    }

    void MasterReceiveSummaryStats()    {

        MeanNcluster = 0;
        MeanAAEntropy = 0;
        MasterReceiveAdditive(MeanNcluster);
        MasterReceiveAdditive(MeanAAEntropy);
        MeanNcluster /= GetLocalNgene();
        MeanAAEntropy /= GetLocalNgene();
    }
};

