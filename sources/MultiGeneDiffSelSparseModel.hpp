
// this is a multigene version of singleomegamodel
//
// - branch lengths are shared across genes, and are iid Exponential of rate lambda
// - nucleotide relative exchangeabilities and stationaries are also shared across genes (uniform Dirichlet)
// - the array of gene-specific omega's are iid gamma with hyperparameters omegahypermean and omegahyperinvshape
//
// the sequence of MCMC moves is as follows:
// - genes resample substitution histories, gather path suff stats and move their omega's
// - master receives the array of omega's across genes, moves their hyperparameters and then broadcast the new value of these hyperparams
// - master collects branch length suff stats across genes, moves branch lengths and broadcasts their new value
// - master collects nuc path suffstats across genes, moves nuc rates and broadcasts their new value

#include "DiffSelSparseModel.hpp"
#include "Parallel.hpp"
#include "MultiGeneProbModel.hpp"
#include "IIDMultiBernBeta.hpp"
#include "Chrono.hpp"

class MultiGeneDiffSelSparseModel : public MultiGeneProbModel {

    private:

	Tree* tree;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

    // branch lengths are shared across genes
    // iid expo (gamma of shape 1 and scale lambda)
    // where lambda is a hyperparameter
	double lambda;
	BranchIIDGamma* branchlength;
	
    // across genes:
    int Ncond;
    int Nlevel;
    int codonmodel;
    
    vector<double> shiftprobhypermean;
    vector<double> shiftprobhyperinvconc;
    vector<double> pi;
    double pihypermean;
    double pihyperinvconc;

    // IIDMultiBernBeta* shiftprobarray;
    vector<int> totcount;
    IIDMultiCount* shiftcountarray;

    // suffstats for paths, as a function of branch lengths
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
    // suff stats for branch lengths, as a function of lambda
	GammaSuffStat hyperlengthsuffstat;

    // each gene defines its own DiffSelSparseModel
    std::vector<DiffSelSparseModel*> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    public:

    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneDiffSelSparseModel(string datafile, string intreefile, int inNcond, int inNlevel, int incodonmodel, int inmyid, int innprocs) : MultiGeneProbModel(inmyid,innprocs)  {

        codonmodel = incodonmodel;
        Ncond = inNcond;
        Nlevel = inNlevel;

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

        GeneLogPrior = 0;
        lnL = 0;
    }

    void Allocate() {

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree,1.0,lambda);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        shiftprobhypermean.assign(Ncond-1,0.5);
        shiftprobhyperinvconc.assign(Ncond-1,0.5);
        pihypermean = 0.1;
        pihyperinvconc = 0.1;
        pi.assign(Ncond-1,0.1);
        // shiftprobarray = new IIDMultiBernBeta(GetLocalNgene(),pi,shiftprobhypermean,shiftprobhyperinvconc);
        totcount.assign(GetLocalNgene(),0);
        for (int gene=0; gene<GetLocalNgene(); gene++)  {
            totcount[gene] = GetLocalGeneNsite(gene) * Naa;
        }
        shiftcountarray = new IIDMultiCount(totcount,pi,shiftprobhypermean,shiftprobhyperinvconc);

        lnL = 0;
        GeneLogPrior = 0;

        if (! GetMyid())    {
            geneprocess.assign(0,(DiffSelSparseModel*) 0);
        }
        else    {
            geneprocess.assign(GetLocalNgene(),(DiffSelSparseModel*) 0);

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene] = new DiffSelSparseModel(GetLocalGeneName(gene),treefile,Ncond,Nlevel,codonmodel);
                geneprocess[gene]->SetBLMode(2);
                geneprocess[gene]->SetNucMode(0);
                geneprocess[gene]->SetFitnessHyperMode(3);
                geneprocess[gene]->Allocate();
            }
        }
    }

    void FastUpdate()   {
        branchlength->SetScale(lambda);
    }

    void MasterUpdate() override {
        FastUpdate();
        if (nprocs > 1) {
            MasterSendGlobalBranchLengths();
            MasterSendShiftProbHyperParameters();
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {
        SlaveReceiveGlobalBranchLengths();
        SlaveReceiveShiftProbHyperParameters();
        GeneUpdate();
        SlaveSendLogProbs();
    }

    void GeneUpdate()	{
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->Update();
        }
    }

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    int GetNbranch() const {
        return tree->GetNbranch();
    }

    const Tree* GetTree() const {
        return tree;
    }

    //-------------------
    // Traces and Monitors
    //-------------------

    void TraceHeader(ostream& os) const {

        os << "#logprior\tlnL\tlength\t";
        for (int k=1;k<Ncond;k++)   {
            os << "pi" << k << '\t';
            os << "mean\t";
            os << "invconc\t";
        }
        os << '\n';
    }

    void Trace(ostream& os) const {
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        for (int k=1;k<Ncond;k++)   {
            os << pi[k-1] << '\t';
            os << shiftprobhypermean[k-1] << '\t';
            os << shiftprobhyperinvconc[k-1] << '\t';
        }
        os << '\n';
		os.flush();
    }

	void Monitor(ostream& os) const {}

    /*
	void MasterToStream(ostream& os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << shiftprobhypermean << '\n';
        os << shiftprobhyperinvconc << '\n';
        os << pi << '\n';

        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Status stat;
            int size;
            MPI_Recv(&size,1,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPIBuffer buffer(size);
            MPI_Recv(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            os << size << '\n';
            buffer.ToStream(os);
        }
    }

    void SlaveToStream() const override {

        int size = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            size += geneprocess[gene]->GetMPISize();
        }
        MPIBuffer buffer(size);
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            buffer << *geneprocess[gene];
        }
        MPI_Send(&size,1,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }

	void MasterFromStream(istream& is) {
        is >> lambda;
        is >> *branchlength;
        is >> shiftprobhypermean;
        is >> shiftprobhyperinvconc;
        is >> pi;

        for (int proc=1; proc<GetNprocs(); proc++)  {
            int size;
            is >> size;
            MPI_Send(&size,1,MPI_INT,proc,TAG1,MPI_COMM_WORLD);
            MPIBuffer buffer(size);
            buffer.FromStream(is);
            MPI_Send(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
        }
    }

    void SlaveFromStream() override {

        int checksize = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            checksize += geneprocess[gene]->GetMPISize();
        }
        MPI_Status stat;
        int size;
        MPI_Recv(&size,1,MPI_INT,0,TAG1,MPI_COMM_WORLD,&stat);
        if (size != checksize)  {
            cerr << "error in SlaveFromStream: non matching buffer size\n";
            exit(1);
        }
        MPIBuffer buffer(size);
        MPI_Recv(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            buffer >> *geneprocess[gene];
        }
    }
    */

    double GetSlaveMoveTime() const {
        return moveTime;
    }

    double GetSlaveMapTime() const   {
        return mapTime;
    }

    double GetMasterMoveTime() const {
        return movechrono.GetTime();
    }

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
        total += ShiftProbHyperLogPrior();
        total += GeneLogPrior;
        if (std::isnan(total))   {
            cerr << "GetLogPrior is nan\n";
            exit(1);
        }
		return total;
    }

	double BranchLengthsHyperLogPrior() const {
		return -lambda / 10;
	}

	double BranchLengthsLogPrior() const {
		return branchlength->GetLogProb();
	}

    double ShiftProbHyperLogPrior() const   {
        double total = 0;
        for (int k=1; k<Ncond; k++) {
            total += GetPiLogPrior(k);
            total += GetShiftProbHyperMeanLogPrior(k);
            total += GetShiftProbHyperInvConcLogPrior(k);
        }
        return total;
    }

    double GetPiLogPrior(int k) const {
        double alpha = pihypermean / pihyperinvconc;
        double beta = (1-pihypermean) / pihyperinvconc;
        return Random::logBetaDensity(pi[k-1],alpha,beta);
    }

    double GetShiftProbHyperMeanLogPrior(int k) const {
        return 0;
    }

    double GetShiftProbHyperInvConcLogPrior(int k) const {
        return -shiftprobhyperinvconc[k-1];
    }

    double GetCountLogProb(int k) const {
        return shiftcountarray->GetMarginalLogProb(k-1);
    }

    double GetCountLogProb() const {
        double total = 0;
        for (int k=1; k<Ncond; k++) {
            total += GetCountLogProb(k);
        }
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
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // log prob for moving branch lengths hyperparameter (lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //-------------------
    // Moves
    //-------------------

    // all methods starting with Master are called only by master
    // for each such method, there is a corresponding method called by slave, and starting with Slave
    //
    // all methods starting with Gene are called only be slaves, and do some work across all genes allocated to that slave

    void MasterMove() override {

        int nrep = 1;
		// int nrep = 3;

		for (int rep=0; rep<nrep; rep++)	{

            MasterReceiveShiftCounts();
            movechrono.Start();
            MoveShiftProbHyperParameters(3);
            movechrono.Stop();
            MasterSendShiftProbHyperParameters();

            MasterReceiveLengthSuffStat();
            movechrono.Start();
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();
            movechrono.Stop();
            MasterSendGlobalBranchLengths();
        }

        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {

        movechrono.Start();
        mapchrono.Start();
        GeneResampleSub(1.0);
        mapchrono.Stop();
        movechrono.Stop();

		int nrep = 1;
		// int nrep = 3;

		for (int rep=0; rep<nrep; rep++)	{

            movechrono.Start();
            GeneMove();
            movechrono.Stop();

            SlaveSendShiftCounts();
            SlaveReceiveShiftProbHyperParameters();
            movechrono.Start();
            GeneResampleShiftProbs();
            movechrono.Stop();

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalBranchLengths();
        }

        SlaveSendLogProbs();
    }

    void GeneCollectShiftCounts()    {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            for (int k=1; k<Ncond; k++) {
                (*shiftcountarray)[gene][k-1] = geneprocess[gene]->GetNshift(k);
            }
        }
    }

    void GeneResampleShiftProbs()    {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleShiftProb();
        }
    }

    void GeneMove() {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->MoveParameters(1,10);
            // geneprocess[gene]->MoveParameters(1,20);
        }
    }

    void GeneResampleSub(double frac)  {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void ResampleBranchLengths()    {
		branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

	void MoveBranchLengthsHyperParameter()	{

		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);

        ScalingMove(lambda,1.0,10,&MultiGeneDiffSelSparseModel::BranchLengthsHyperLogProb,&MultiGeneDiffSelSparseModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&MultiGeneDiffSelSparseModel::BranchLengthsHyperLogProb,&MultiGeneDiffSelSparseModel::NoUpdate,this);

		branchlength->SetScale(lambda);
    }

    void MoveShiftProbHyperParameters(int nrep) {
        MovePi(0.3,nrep);
        MoveShiftProbHyperMean(0.3,nrep);
        MoveShiftProbHyperInvConc(0.3,nrep);
    }

    double MovePi(double tuning, int nrep)   {
        double nacc = 0;
        double ntot = 0;
        for (int k=1; k<Ncond; k++) {
            for (int rep=0; rep<nrep; rep++)    {
                double deltalogprob = - GetPiLogPrior(k) - GetCountLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double bk = pi[k-1];
                pi[k-1] += m;
                while ((pi[k-1] < 0) || (pi[k-1] > 1))  {
                    if (pi[k-1] < 0)    {
                        pi[k-1] = - pi[k-1];
                    }
                    if (pi[k-1] > 1)    {
                        pi[k-1] = 2 - pi[k-1];
                    }
                }
                deltalogprob += GetPiLogPrior(k) + GetCountLogProb(k);

                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc)    {
                    nacc++;
                }
                else    {
                    pi[k-1] = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveShiftProbHyperMean(double tuning, int nrep)  {
        double nacc = 0;
        double ntot = 0;
        for (int k=1; k<Ncond; k++) {
            for (int rep=0; rep<nrep; rep++)    {
                double deltalogprob = - GetShiftProbHyperMeanLogPrior(k) - GetCountLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double bk = shiftprobhypermean[k-1];
                shiftprobhypermean[k-1] += m;
                while ((shiftprobhypermean[k-1] < 0) || (shiftprobhypermean[k-1] > 1))  {
                    if (shiftprobhypermean[k-1] < 0)    {
                        shiftprobhypermean[k-1] = - shiftprobhypermean[k-1];
                    }
                    if (shiftprobhypermean[k-1] > 1)    {
                        shiftprobhypermean[k-1] = 2 - shiftprobhypermean[k-1];
                    }
                }
                deltalogprob += GetShiftProbHyperMeanLogPrior(k) + GetCountLogProb(k);

                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc)    {
                    nacc++;
                }
                else    {
                    shiftprobhypermean[k-1] = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveShiftProbHyperInvConc(double tuning, int nrep)  {
        double nacc = 0;
        double ntot = 0;
        for (int k=1; k<Ncond; k++) {
            for (int rep=0; rep<nrep; rep++)    {
                double deltalogprob = - GetShiftProbHyperInvConcLogPrior(k) - GetCountLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                double bk = shiftprobhyperinvconc[k-1];
                shiftprobhyperinvconc[k-1] *= e;
                deltalogprob += GetShiftProbHyperInvConcLogPrior(k) + GetCountLogProb(k);
                deltalogprob += m;

                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc)    {
                    nacc++;
                }
                else    {
                    shiftprobhyperinvconc[k-1] = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
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

        moveTime = movechrono.GetTime();
        mapTime = mapchrono.GetTime();
        SlaveSendAdditive(moveTime);
        SlaveSendAdditive(mapTime);
    }

    void MasterReceiveLogProbs()    {
        GeneLogPrior = 0;
        MasterReceiveAdditive(GeneLogPrior);
        lnL = 0;
        MasterReceiveAdditive(lnL);
        moveTime = 0;
        mapTime = 0;
        MasterReceiveAdditive(moveTime);
        MasterReceiveAdditive(mapTime);
        moveTime /= (GetNprocs()-1);
        mapTime /= (GetNprocs()-1);
    }

    void MasterReceiveShiftCounts()   {
        MasterReceiveGeneArray(*shiftcountarray);
    }

    void SlaveSendShiftCounts()   {
        GeneCollectShiftCounts();
        SlaveSendGeneArray(*shiftcountarray);
    }

    void MasterSendShiftProbHyperParameters()   {
        MasterSendGlobal(shiftprobhypermean,shiftprobhyperinvconc);
        MasterSendGlobal(pi);
    }

    void SlaveReceiveShiftProbHyperParameters() {
        SlaveReceiveGlobal(shiftprobhypermean,shiftprobhyperinvconc);
        SlaveReceiveGlobal(pi);
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetShiftProbHyperParameters(pi,shiftprobhypermean,shiftprobhyperinvconc);
        }
    }
};

