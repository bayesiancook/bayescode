
#include "AAMutSelHyperSBDPOmegaModel.hpp"
#include "Parallel.hpp"
#include "MultiGeneProbModel.hpp"

class MultiGeneAAMutSelHyperSBDPOmegaModel : public MultiGeneProbModel {

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

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

    // mixture components
    // set of Ncat Dirichlet densities
    // centers
    vector<double> aacenterhypercenter;
    double aacenterhyperinvconc;
    IIDDirichlet* componentaacenterarray;
    // concentrations
    double aaconchypermean;
    double aaconchyperinvshape;
    IIDGamma* componentaaconcentrationarray;
    // and associated suffstatarray
    DirichletSuffStatArray* aahypersuffstatarray;

    // now, a mixture model drawing from this set of Ncat components
    // weights:
    double kappa;
    StickBreakingProcess* weight;
    OccupancySuffStat* occupancy;

    std::vector<AAMutSelHyperSBDPOmegaModel*> geneprocess;

    double lnL;
    double GeneLogPrior;
    double MeanStatEnt;
    double MeanAAConc;
    double MeanAACenterEnt;

    int Ncat;

    int fixomega;

    public:

    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneAAMutSelHyperSBDPOmegaModel(string datafile, string intreefile, int inNcat, int infixomega, int inmyid, int innprocs) : MultiGeneProbModel(inmyid,innprocs) {

        AllocateAlignments(datafile);
        treefile = intreefile;
        Ncat = inNcat;
        fixomega = infixomega;

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

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omegaarray = new IIDGamma(GetLocalNgene(),1.0,1.0);
        if (fixomega) {
            for (int i=0; i<GetLocalNgene(); i++)   {
                (*omegaarray)[i] = 1.0;
            }
        }

        lnL = 0;
        GeneLogPrior = 0;
        MeanStatEnt = 0;
        MeanAAConc = 0;
        MeanAACenterEnt = 0;

        // mixture components
        aacenterhypercenter.assign(Naa,1.0/Naa);
        aacenterhyperinvconc = 1.0/Naa;
        componentaacenterarray = new IIDDirichlet(Ncat,aacenterhypercenter,1.0/aacenterhyperinvconc);
        componentaacenterarray->SetUniform();

        aaconchypermean = Naa;
        aaconchyperinvshape = 1.0;
        double alpha = 1.0 / aaconchyperinvshape;
        double beta = alpha / aaconchypermean;
        componentaaconcentrationarray = new IIDGamma(Ncat,alpha,beta);
        for (int k=0; k<Ncat; k++)  {
            (*componentaaconcentrationarray)[k] = 20.0;
        }

        // suff stats for mixture components
        aahypersuffstatarray = new DirichletSuffStatArray(Ncat,Naa);

        // mixture weights
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat,kappa);
        occupancy = new OccupancySuffStat(Ncat);

        if (! GetMyid())    {
            geneprocess.assign(0,(AAMutSelHyperSBDPOmegaModel*) 0);
        }
        else    {
            geneprocess.assign(GetLocalNgene(),(AAMutSelHyperSBDPOmegaModel*) 0);

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene] = new AAMutSelHyperSBDPOmegaModel(GetLocalGeneName(gene),treefile,Ncat);
                geneprocess[gene]->SetFixAAHyperMix(1);
                geneprocess[gene]->SetFixOmega(fixomega);
            }
        }
    }

    void Unfold()   {

        if (! GetMyid())    {
            MasterSendGlobalBranchLengths();
            // MasterSendGlobalNucRates();
            if (! fixomega) {
                MasterSendOmegaHyperParameters();
                MasterSendOmega();
            }
            MasterSendAAHyperMixture();
            MasterReceiveLogProbs();
        }
        else    {

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->Allocate();
            }

            SlaveReceiveGlobalBranchLengths();
            // SlaveReceiveGlobalNucRates();
            if (! fixomega) {
                SlaveReceiveOmegaHyperParameters();
                SlaveReceiveOmega();
            }
            SlaveReceiveAAHyperMixture();

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
        os << "ncluster\t";
        os << "aastatent\t";
        os << "hypercenter\t";
        os << "hyperstatalpha\n";
    }

    void Trace(ostream& os) const {
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << 3*branchlength->GetTotalLength() << '\t';
        os << omegaarray->GetMean() << '\t';
        os << omegaarray->GetVar() << '\t';
        os << GetNcluster() << '\t';
        os << MeanStatEnt << '\t';
        os << MeanAAConc << '\t';
		os << MeanAACenterEnt << '\n';
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
		double total = GeneLogPrior;
		total += BranchLengthsHyperLogPrior();
		total += BranchLengthsLogPrior();
        // total += NucRatesLogPrior();
        if (! fixomega) {
            total += OmegaHyperLogPrior();
            total += OmegaLogPrior();
        }
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AAHyperLogPrior();
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

    double StickBreakingHyperLogPrior() const   {
        return -kappa/10;
    }

    double StickBreakingLogPrior() const    {
        return weight->GetLogProb(kappa);
    }

    double AAHyperLogPrior() const {
        double total = 0;
        total += componentaacenterarray->GetLogProb();
        total += componentaaconcentrationarray->GetLogProb();
        return total;
    }

    double AAHyperLogPrior(int k) const {
        double total = 0;
        total += componentaacenterarray->GetLogProb(k);
        total += componentaaconcentrationarray->GetLogProb(k);
        return total;
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

    double AAHyperSuffStatLogProb(int k) const   {
        return aahypersuffstatarray->GetVal(k).GetLogProb(componentaacenterarray->GetVal(k),componentaaconcentrationarray->GetVal(k));
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

    // for moving aa hyper params (aacenter and aainvconc)
    // for component k of the mixture
    double AAHyperLogProb(int k) const   {
        return AAHyperLogPrior(k) + AAHyperSuffStatLogProb(k);
    }

    // for moving kappa
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    // Moves
    //-------------------

    void MasterMove() override {

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            MasterReceiveAAHyperSuffStat();
            MoveAAHyperMixture();
            MasterSendAAHyperMixture();

            if (! fixomega) {
                MasterReceiveOmega();
                MoveOmegaHyperParameters();
                MasterSendOmegaHyperParameters();
            }

            MasterReceiveLengthSuffStat();
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();
            MasterSendGlobalBranchLengths();
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

            MoveGeneAA();
            SlaveSendAAHyperSuffStat();
            SlaveReceiveAAHyperMixture();

            if (! fixomega) {
                MoveGeneOmegas();
                SlaveSendOmega();
                SlaveReceiveOmegaHyperParameters();
            }

            MoveGeneNucRates();

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalBranchLengths();
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

    void MoveGeneAA()   {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->MoveAA(3);
        }
    }

    void MoveGeneNucRates() {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->MoveNucRates();
        }
    }

    /*
    void MoveAAHyperMixtureHyperParameters()    {
    }
    */

    void MoveAAHyperMixture()   {
        for (int rep=0; rep<3; rep++)  {
            MoveAAHyperMixtureComponents();
            // LabelSwitchingMove();
            ResampleWeights();
            MoveKappa();
            // MoveAAHyperMixtureHyperParameters();
        }
    }

    void MoveAAHyperMixtureComponents() {

        for (int k=0; k<Ncat; k++)  {
            for (int rep=0; rep<10; rep++)    {
                MoveAAHyperCenters(1.0,1);
                MoveAAHyperCenters(1.0,3);
                MoveAAHyperCenters(0.3,3);
                MoveAAHyperConcentrations(1.0);
                MoveAAHyperConcentrations(0.3);
            }
        }
    }


    double MoveAAHyperCenters(double tuning, int n) {
		double nacc = 0;
		double ntot = 0;
        vector<double> bk(Naa,0);
        for (int k=0; k<Ncat; k++)  {
            vector<double>& aa = (*componentaacenterarray)[k];
            bk = aa;
            double deltalogprob = -AAHyperLogProb(k);
            double loghastings = Random::ProfileProposeMove(aa,Naa,tuning,n);
            deltalogprob += loghastings;
            deltalogprob += AAHyperLogProb(k);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                aa = bk;
            }
            ntot++;
        }
		return nacc/ntot;
	}

    double MoveAAHyperConcentrations(double tuning)  {
		double nacc = 0;
		double ntot = 0;
        for (int k=0; k<Ncat; k++)  {
            double& c = (*componentaaconcentrationarray)[k];
            double bk = c;
            double deltalogprob = -AAHyperLogProb(k);
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            c *= e;
            deltalogprob += m;
            deltalogprob += AAHyperLogProb(k);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                c = bk;
            }
            ntot++;
        }
		return nacc/ntot;
    }

    /*
    void LabelSwitchingMove()   {
        MoveOccupiedCompAlloc(5);
        MoveAdjacentCompAlloc(5);
    }

    double MoveOccupiedCompAlloc(int k0)	{

        const vector<double>& V = weight->GetBetaVariates();
        const vector<double>& w = weight->GetArray();

        int nrep = (int) (k0 * kappa);
        ResampleWeights();
        double total = 0.0;
        int Nocc = GetNcluster();
        if (Nocc != 1)	{
            for (int i=0; i<nrep; i++)	{
                int occupiedComponentIndices[Nocc];
                int j=0;
                for (int k=0; k<Ncat; k++)	{
                    if ((*occupancy)[k] != 0)	{
                        occupiedComponentIndices[j] = k;
                        j++;
                    }
                }
                if (j != Nocc)	{
                    cerr << "error in MoveOccupiedCompAlloc.\n";
                    exit(1);
                }
                int indices[2];
                Random::DrawFromUrn(indices,2,Nocc);
                int cat1 = occupiedComponentIndices[indices[0]];
                int cat2 = occupiedComponentIndices[indices[1]];
                double logMetropolis = ((*occupancy)[cat2] - (*occupancy)[cat1]) * log(w[cat1] / w[cat2]);
                int accepted = (log(Random::Uniform()) < logMetropolis);
                if (accepted)	{
                    total += 1.0;
                    componentaacenterarray->Swap(cat1,cat2);
                    componentaaconcentrationarray->Swap(cat1,cat2);
                    sitealloc->SwapComponents(cat1,cat2);
                    occupancy->Swap(cat1,cat2);
                }
            }
            return total /= nrep;
        }
        return 0;
    }

    double MoveAdjacentCompAlloc(int k0)	{

        ResampleWeights();
        int nrep = (int) (k0 * kappa);
        
        double total = 0;

        const vector<double>& V = weight->GetBetaVariates();

        for (int i=0; i<nrep; i++)	{
            int cat1 = (int)(Random::Uniform() * (Ncat-2));  
            int cat2 = cat1 + 1;
            double logMetropolis = ((*occupancy)[cat1] * log(1 - V[cat2])) - ((*occupancy)[cat2] * log(1-V[cat1]));
            int accepted = (log(Random::Uniform()) < logMetropolis);
            if (accepted)	{
                total += 1.0;
                componentaacenterarray->Swap(cat1,cat2);
                componentaaconcentrationarray->Swap(cat1,cat2);
                sitealloc->SwapComponents(cat1,cat2);
                weight->SwapComponents(cat1,cat2);
                occupancy->Swap(cat1,cat2);
            }
        }

        return total /= nrep;
    }
    */

    void ResampleWeights()  {
        weight->GibbsResample(*occupancy);
    }

    void MoveKappa()    {
        ScalingMove(kappa,1.0,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::StickBreakingHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);
        ScalingMove(kappa,0.3,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::StickBreakingHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);
    }

    int GetNcluster() const {

        int n = 0;
        for (int i=0; i<Ncat; i++)  {
            if (occupancy->GetVal(i))    {
                n++;
            }
        }
        return n;
    }

    void ResampleBranchLengths()    {
		branchlength->GibbsResample(*lengthsuffstatarray);
    }

	void MoveBranchLengthsHyperParameter()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);

        ScalingMove(lambda,1.0,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::BranchLengthsHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::BranchLengthsHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);

		branchlength->SetScale(lambda);
    }

    void MoveOmegaHyperParameters()  {

		omegahypersuffstat.Clear();
		omegaarray->AddSuffStat(omegahypersuffstat);

        ScalingMove(omegahypermean,1.0,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);
        ScalingMove(omegahypermean,0.3,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);
        ScalingMove(omegahyperinvshape,1.0,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);
        ScalingMove(omegahyperinvshape,0.3,10,&MultiGeneAAMutSelHyperSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelHyperSBDPOmegaModel::NoUpdate,this);

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

    // aa hyper mixture

    void MasterSendAAHyperMixture() {
        MasterSendGlobal(*componentaacenterarray,*componentaaconcentrationarray);
        MasterSendGlobal(*weight);
    }

    void SlaveReceiveAAHyperMixture()   {
        SlaveReceiveGlobal(*componentaacenterarray,*componentaaconcentrationarray);
        SlaveReceiveGlobal(*weight);
        for (int gene=0; gene<GetLocalNgene(); gene++)    {
            geneprocess[gene]->SetAAHyperMixture(*componentaacenterarray,*componentaaconcentrationarray,*weight);
        }
    }

    // aa hyper suff stat

    void SlaveSendAAHyperSuffStat() {
        aahypersuffstatarray->Clear();
        occupancy->Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectAAHyperSuffStat();
            aahypersuffstatarray->Add(*geneprocess[gene]->GetAAHyperSuffStatArray());
            geneprocess[gene]->UpdateOccupancies();
            occupancy->Add(*geneprocess[gene]->GetOccupancies());
        }
        SlaveSendAdditive(*aahypersuffstatarray);
        SlaveSendAdditive(*occupancy);
    }

    void MasterReceiveAAHyperSuffStat() {

        aahypersuffstatarray->Clear();
        occupancy->Clear();
        MasterReceiveAdditive(*aahypersuffstatarray);
        MasterReceiveAdditive(*occupancy);
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

    // log probs

    void SlaveSendLogProbs()   {

        GeneLogPrior = 0;
        lnL = 0;
        MeanStatEnt = 0;
        MeanAAConc = 0;
        MeanAACenterEnt = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
            MeanStatEnt += geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanAAEntropy();
            MeanAAConc += geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanComponentAAConcentration();
            MeanAACenterEnt += geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanComponentAAEntropy();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
        SlaveSendAdditive(MeanStatEnt);
        SlaveSendAdditive(MeanAAConc);
        SlaveSendAdditive(MeanAACenterEnt);
    }

    void MasterReceiveLogProbs()    {

        GeneLogPrior = 0;
        lnL = 0;
        MasterReceiveAdditive(GeneLogPrior);
        MasterReceiveAdditive(lnL);
        MeanStatEnt = 0;
        MeanAAConc = 0;
        MeanAACenterEnt = 0;
        MasterReceiveAdditive(MeanStatEnt);
        MasterReceiveAdditive(MeanAAConc);
        MasterReceiveAdditive(MeanAACenterEnt);
        MeanStatEnt /= GetTotNsite();
        MeanAAConc /= GetTotNsite();
        MeanAACenterEnt /= GetTotNsite();
    }
};

