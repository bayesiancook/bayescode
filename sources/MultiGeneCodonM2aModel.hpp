
#include "CodonM2aModel.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"

#include "Chrono.hpp"


class MultiGeneCodonM2aModel : public MultiGeneMPIModule	{

    private:

	Tree* tree;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

    double puromhypermean;
	double puromhyperinvconc;
	IIDBeta* puromarray;
	BetaSuffStat puromsuffstat;

    double dposomhypermean;
    double dposomhyperinvshape;
    IIDGamma* dposomarray;
    GammaSuffStat dposomsuffstat;

    double purwhypermean;
    double purwhyperinvconc;
    IIDBeta* purwarray;
    BetaSuffStat purwsuffstat;

    double poswhypermean;
    double poswhyperinvconc;
    IIDBernoulliBeta* poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double pi;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    NucPathSuffStat nucpathsuffstat;

    std::vector<CodonM2aModel*> geneprocess;

    double totlnL;
    double* lnL;

    int burnin;

    Chrono timepercycle;
    Chrono omegachrono,hyperchrono,mastersampling;

    public:

	CodonStateSpace* GetCodonStateSpace()   {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    MultiGeneCodonM2aModel(string datafile, string intreefile, double inpihypermean, double inpihyperinvconc, int inmyid, int innprocs) : MultiGeneMPIModule(inmyid,innprocs) {

        burnin = 0;

        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;

        AllocateAlignments(datafile);
        treefile = intreefile;

        // all datafiles have all taxa (with missing data if needed) in same order
        // makes it easier to register tree with data, etc.

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

        timepercycle.Reset();
        mastersampling.Reset();
        // Allocate();
    }

    void Allocate() {

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree,1.0,lambda);
        lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelrate.assign(Nrr,0);
        double totrr = 0;
        for (int k=0; k<Nrr; k++)	{
            nucrelrate[k] = Random::sExpo();
            totrr += nucrelrate[k];
        }
        for (int k=0; k<Nrr; k++)	{
            nucrelrate[k] /= totrr;
        }

        nucstat.assign(Nnuc,0);
        double totstat = 0;
        for (int k=0; k<Nnuc; k++)	{
            nucstat[k] = Random::sGamma(1.0);
            totstat += nucstat[k];
        }
        for (int k=0; k<Nnuc; k++)	{
            nucstat[k] /= totstat;
        }
        nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        puromhypermean = 0.5;
        puromhyperinvconc = 0.5;

        dposomhypermean = 1.0;
        dposomhyperinvshape = 1.0;

        purwhypermean = 0.5;
        purwhyperinvconc = 0.5;

        pi = pihypermean;
        if (! pi)   {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }
        else    {
            poswhypermean = 0.1;
            poswhyperinvconc = 1;
        }

        double puromalpha = puromhypermean / puromhyperinvconc;
        double purombeta = (1-puromhypermean) / puromhyperinvconc;
        puromarray = new IIDBeta(GetLocalNgene(),puromalpha,purombeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray = new IIDGamma(GetLocalNgene(),dposomalpha,dposombeta);

        double purwalpha = purwhypermean / purwhyperinvconc;
        double purwbeta = (1-purwhypermean) / purwhyperinvconc;
        purwarray = new IIDBeta(GetLocalNgene(),purwalpha,purwbeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        poswarray = new IIDBernoulliBeta(GetLocalNgene(),pi,poswalpha,poswbeta);

        if (myid)   {
            lnL = new double[GetLocalNgene()];
        }
        else    {
            lnL = 0;
        }


        if (! GetMyid())    {
            geneprocess.assign(0,(CodonM2aModel*) 0);
        }
        else    {
            geneprocess.assign(GetLocalNgene(),(CodonM2aModel*) 0);

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene] = new CodonM2aModel(GetLocalGeneName(gene),treefile,pi);
            }
        }
    }

    void Unfold()   {

        if (! GetMyid())    {
            MasterSendHyperParameters();
            MasterSendGlobalParameters();
            MasterSendMixture();
            MasterReceiveLogLikelihood();
        }
        else    {
            SlaveReceiveHyperParameters();
            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->Allocate();
            }

            SlaveReceiveGlobalParameters();
            SlaveReceiveMixture();

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->UpdateMatrices();
                geneprocess[gene]->Unfold();
            }
            SlaveSendLogLikelihood();
        }
    }

    void TraceHeader(ostream& os)   {

        os << "#time\t";
        os << "logprior\tlnL\tlength\t";
        os << "pi\t";
        os << "nposfrac\t";
        os << "purommean\tinvconc\t";
        os << "dposommean\tinvshape\t";
        os << "purwmean\tinvconc\t";
        os << "poswmean\tinvconc\t";
        os << "statent\t";
        os << "rrent\n";
        timepercycle.Start();
    }

    void Trace(ostream& os)    {
        timepercycle.Stop();
        os << timepercycle.GetTime() << '\t';
        // cerr << myid << '\t' << timepercycle.GetTime() << '\t' << mastersampling.GetTime() << '\t' << omegachrono.GetTime() << '\t' << hyperchrono.GetTime() << '\n';
        mastersampling.Reset();
        omegachrono.Reset();
        hyperchrono.Reset();
        timepercycle.Reset();
        timepercycle.Start();
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
        os << pi << '\t';
        os << GetNpos() << '\t';
        os << puromhypermean << '\t' << puromhyperinvconc << '\t';
        os << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
        os << purwhypermean << '\t' << purwhyperinvconc << '\t';
        os << poswhypermean << '\t' << poswhyperinvconc << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
		os.flush();
    }

    void TraceGlobalParameters(ostream& os)   {
        for (int j=0; j<Nbranch; j++)   {
            os << branchlength->GetVal(j) << '\t';
        }
        for (int j=0; j<Nrr; j++)   {
            os << nucrelrate[j] << '\t';
        }
        for (int j=0; j<Nnuc; j++)  {
            os << nucstat[j] << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TraceHyperParameters(ostream& os)    {
        os << puromhypermean << '\t';
        os << puromhyperinvconc << '\t';
        os << dposomhypermean << '\t';
        os << dposomhyperinvshape << '\t';
        os << purwhypermean << '\t';
        os << purwhyperinvconc << '\t';
        os << pi << '\t';
        os << poswhypermean << '\t';
        os << poswhyperinvconc << '\t';
        os << '\n';
        os.flush();
    }

    void TracePosWeight(ostream& os) {

        for (int gene=0; gene<Ngene; gene++)    {
            os << (*poswarray)[gene] << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TracePosOm(ostream& os) {

        for (int gene=0; gene<Ngene; gene++)    {
            os << 1 + (*dposomarray)[gene] << '\t';
        }
        os << '\n';
        os.flush();
    }

    void SlaveTracePostProbHeader(string name)    {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            ofstream os((name + "_" + GetLocalGeneName(gene) + ".sitepp").c_str());
        }
    }

    void SlaveTracePostProb(string name)    {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            ofstream os((name + "_" + GetLocalGeneName(gene) + ".sitepp").c_str());
            geneprocess[gene]->TracePostProb(os);
        }
    }

    int GetNpos()    {
        return GetNgene() - poswarray->GetNullSet();
    }

    double GetLogPrior()    {
		double total = 0;
		total += LambdaLogProb();
		total += LengthLogProb();
		total += HyperLogPrior();
        total += NucLogPrior();
		return total;
    }

    double NucLogPrior()    {
        return 0;
    }

	double LambdaLogProb()	{
		return -lambda / 10;
	}

	double LengthSuffStatLogProb()	{
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double HyperLogPrior()   {
        double total = 0;
        total -= puromhyperinvconc;
        total -= purwhyperinvconc;
        total -= 10*poswhyperinvconc;
        total -= dposomhypermean;
        total -= 10*dposomhyperinvshape;
        return total;
    }

    double HyperSuffStatLogProb()   {

        double total = 0;

        double puromalpha = puromhypermean / puromhyperinvconc;
        double purombeta = (1-puromhypermean) / puromhyperinvconc;
        total += puromsuffstat.GetLogProb(puromalpha,purombeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        total += dposomsuffstat.GetLogProb(dposomalpha,dposombeta);

        double purwalpha = purwhypermean / purwhyperinvconc;
        double purwbeta = (1-purwhypermean) / purwhyperinvconc;
        total += purwsuffstat.GetLogProb(purwalpha,purwbeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        total += poswsuffstat.GetLogProb(pi,poswalpha,poswbeta);

        if (isnan(total))   {
            cerr << "hyper suff stat log prob is nan\n";
            exit(1);
        }

        return total;
    }

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

    double GetLogLikelihood()   {
        return totlnL;
    }

	double GetTotalLength()	{
		double tot = 0;
		for (int j=1; j<Nbranch; j++)	{
			tot += branchlength->GetVal(j);
		}
		return tot;
	}

	double GetEntropy(const std::vector<double>& profile, int dim) const {
		double tot = 0;
		for (int i=0; i<dim; i++)	{
			tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
		}
		return tot;
	}

	void Monitor(ostream& os) {}
	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

    void MasterMove() {

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            if (rep)    {
                mastersampling.Start();
            }
            MasterReceiveLengthSuffStat();
            MasterResampleBranchLengths();
            MasterMoveLambda();

            MasterSendGlobalParameters();
            omegachrono.Start();
            MasterReceiveMixtureHyperSuffStat();
            omegachrono.Stop();
            hyperchrono.Start();
            MasterMoveMixtureHyperParameters();
            hyperchrono.Stop();
            MasterSendHyperParameters();

            MasterReceiveNucPathSuffStat();
            MasterMoveNuc();
            MasterSendGlobalParameters();
            if (rep)    {
                mastersampling.Stop();
            }
        }
        burnin++;
        MasterReceiveMixture();
        MasterReceiveLogLikelihood();
    }

    // slave move
    void SlaveMove() {

        Chrono mapping, sampling;

        mapping.Start();
        SlaveResampleSub();
        mapping.Stop();

		int nrep = 30;

        sampling.Start();
		for (int rep=0; rep<nrep; rep++)	{

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalParameters();

            SlaveCollectPathSuffStat();
            SlaveMoveOmega();
            SlaveSendMixtureHyperSuffStat();
            SlaveReceiveHyperParameters();
            SlaveSetArrays();

            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalParameters();

        }
        sampling.Stop();
        burnin++;
        SlaveSendMixture();
        SlaveSendLogLikelihood();
    }

    void MasterSendGlobalParameters() {

        int N = Nbranch + Nrr + Nnuc;
        double* array = new double[N];
        int i = 0;
        for (int j=0; j<Nbranch; j++)   {
            array[i++] = branchlength->GetVal(j);
        }
        for (int j=0; j<Nrr; j++)   {
            array[i++] = nucrelrate[j];
        }
        for (int j=0; j<Nnuc; j++)  {
            array[i++] = nucstat[j];
        }
        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveGlobalParameters()   {

        int N = Nbranch + Nrr + Nnuc;
        double* array = new double[N];
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        int i = 0;
        for (int j=0; j<Nbranch; j++)   {
            (*branchlength)[j] = array[i++];
        }
        for (int j=0; j<Nrr; j++)   {
            nucrelrate[j] = array[i++];
        }
        for (int j=0; j<Nnuc; j++)  {
            nucstat[j] = array[i++];
        }

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetBranchLengths(*branchlength);
            geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
        }
        delete[] array;
    }

    void MasterSendHyperParameters() {

        int N = 9;
        double* array = new double[N];
        int i = 0;
        array[i++] = puromhypermean;
        array[i++] = puromhyperinvconc;
        array[i++] = dposomhypermean;
        array[i++] = dposomhyperinvshape;
        array[i++] = purwhypermean;
        array[i++] = purwhyperinvconc;
        array[i++] = pi;
        array[i++] = poswhypermean;
        array[i++] = poswhyperinvconc;

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveHyperParameters()   {

        int N = 9;
        double* array = new double[N];
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        int i = 0;
        puromhypermean = array[i++];
        puromhyperinvconc = array[i++];
        dposomhypermean = array[i++];
        dposomhyperinvshape = array[i++];
        purwhypermean = array[i++];
        purwhyperinvconc = array[i++];
        pi = array[i++];
        poswhypermean = array[i++];
        poswhyperinvconc = array[i++];

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetMixtureHyperParameters(puromhypermean,puromhyperinvconc,dposomhypermean,dposomhyperinvshape,pi,purwhypermean,purwhyperinvconc,poswhypermean,poswhyperinvconc);
        }
        delete[] array;
    }

    void SlaveResampleSub()  {

        double frac = 1.0;
        /*
        if (burnin > 10)    {
            frac = 0.2;
        }
        */
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void SlaveSendLengthSuffStat()  {

        lengthsuffstatarray->Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectLengthSuffStat();
            lengthsuffstatarray->Add(*geneprocess[gene]->GetLengthSuffStatArray());
        }

        int* count = new int[Nbranch];
        double* beta = new double[Nbranch];
        lengthsuffstatarray->Push(count,beta);
        MPI_Send(count,Nbranch,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(beta,Nbranch,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] count;
        delete[] beta;
    }

    void MasterReceiveLengthSuffStat()  {

        int* count = new int[Nbranch];
        double* beta = new double[Nbranch];
        lengthsuffstatarray->Clear();
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nbranch,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Nbranch,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            lengthsuffstatarray->Add(count,beta);
        }
        delete[] count;
        delete[] beta;
    }

    void MasterResampleBranchLengths()    {
		branchlength->GibbsResample(*lengthsuffstatarray);
    }

	void MasterMoveLambda()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
		MoveLambda(1.0,10);
		MoveLambda(0.3,10);
		branchlength->SetScale(lambda);
    }

	double MoveLambda(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - LambdaLogProb() - LengthSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			lambda *= e;
			deltalogprob += LambdaLogProb() + LengthSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				lambda /= e;
			}
			ntot++;
		}
		return nacc/ntot;
	}

    void SlaveSendMixtureHyperSuffStat()  {

        int Nint = 1 + 1 + 1 + 2;
        int Ndouble = 2 + 2 + 2 + 2;
        int* count = new int[Nint];
        double* beta = new double[Ndouble];

        int i = 0;
        int d = 0;

		puromsuffstat.Clear();
		puromarray->AddSuffStat(puromsuffstat);
        count[i++] = puromsuffstat.GetN();
        beta[d++] = puromsuffstat.GetSumLog0();
        beta[d++] = puromsuffstat.GetSumLog1();

		dposomsuffstat.Clear();
		dposomarray->AddSuffStat(dposomsuffstat);
		// dposomarray->AddSuffStat(dposomsuffstat,*poswarray);
        count[i++] = dposomsuffstat.GetN();
        beta[d++] = dposomsuffstat.GetSum();
        beta[d++] = dposomsuffstat.GetSumLog();

		purwsuffstat.Clear();
		purwarray->AddSuffStat(purwsuffstat);
        count[i++] = purwsuffstat.GetN();
        beta[d++] = purwsuffstat.GetSumLog0();
        beta[d++] = purwsuffstat.GetSumLog1();

		poswsuffstat.Clear();
		poswarray->AddSuffStat(poswsuffstat);
        count[i++] = poswsuffstat.GetN0();
        count[i++] = poswsuffstat.GetN1();
        beta[d++] = poswsuffstat.GetSumLog0();
        beta[d++] = poswsuffstat.GetSumLog1();

        MPI_Send(count,Nint,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(beta,Ndouble,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

        delete[] count;
        delete[] beta;
    }

    void MasterReceiveMixtureHyperSuffStat()  {

		puromsuffstat.Clear();
		dposomsuffstat.Clear();
		purwsuffstat.Clear();
		poswsuffstat.Clear();

        int Nint = 1 + 1 + 1 + 2;
        int Ndouble = 2 + 2 + 2 + 2;
        int* count = new int[Nint];
        double* beta = new double[Ndouble];

        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nint,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Ndouble,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

            int i = 0;
            int d = 0;

            puromsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
            dposomsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
            purwsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
            poswsuffstat.AddNullSuffStat(count[i++]);
            poswsuffstat.AddPosSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
        }

        delete[] count;
        delete[] beta;
    }

    void MasterMoveMixtureHyperParameters()  {

		HyperSlidingMove(puromhypermean,1.0,10,0,1);
		HyperSlidingMove(puromhypermean,0.3,10,0,1);
		HyperScalingMove(puromhyperinvconc,1.0,10);
		HyperScalingMove(puromhyperinvconc,0.3,10);

		HyperSlidingMove(dposomhypermean,3.0,10,0.5,10);
		HyperSlidingMove(dposomhypermean,1.0,10,0.5,10);
		HyperSlidingMove(dposomhypermean,0.3,10,0.5,10);
		HyperScalingMove(dposomhyperinvshape,1.0,10);
		HyperScalingMove(dposomhyperinvshape,0.3,10);

		HyperSlidingMove(poswhypermean,1.0,10,0,1);
		HyperSlidingMove(poswhypermean,0.3,10,0,1);
		HyperScalingMove(poswhyperinvconc,1.0,10);
		HyperScalingMove(poswhyperinvconc,0.3,10);

		HyperSlidingMove(purwhypermean,1.0,10,0,1);
		HyperSlidingMove(purwhypermean,0.3,10,0,1);
		HyperScalingMove(purwhyperinvconc,1.0,10);
		HyperScalingMove(purwhyperinvconc,0.3,10);

        if (burnin > 10)    {
            if (pihyperinvconc)    {
                ResamplePi();
            }
        }
    }

    void SlaveSetArrays()    {

        double puromalpha = puromhypermean / puromhyperinvconc;
        double purombeta = (1-puromhypermean) / puromhyperinvconc;
        puromarray->SetAlpha(puromalpha);
        puromarray->SetBeta(purombeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray->SetShape(dposomalpha);
        dposomarray->SetScale(dposombeta);
        // dposomarray->PriorResample(*poswarray);
        // necessary after changing some dposom values
        /*
        for (int gene=0; gene<GetLocalNgene(); gene++)    {
            geneprocess[gene]->SetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
        }
        */

        double purwalpha = purwhypermean / purwhyperinvconc;
        double purwbeta = (1-purwhypermean) / purwhyperinvconc;
        purwarray->SetAlpha(purwalpha);
        purwarray->SetBeta(purwbeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        poswarray->SetPi(pi);
        poswarray->SetAlpha(poswalpha);
        poswarray->SetBeta(poswbeta);

    }

    void ResamplePi()   {

        int n0 = poswsuffstat.GetN0();
        int n1 = poswsuffstat.GetN1();
        if ((n0+n1) != Ngene)   {
            cerr << "error in resample pi\n";
            exit(1);
        }
        double pialpha = pihypermean / pihyperinvconc;
        double pibeta = (1-pihypermean) / pihyperinvconc;
        double a0 = Random::sGamma(pialpha + n0);
        double a1 = Random::sGamma(pibeta + n1);
        pi = a1 / (a0 + a1);
    }

	double HyperSlidingMove(double& x, double tuning, int nrep, double min = 0, double max = 0)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - HyperLogPrior() - HyperSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
		    x += m;
            if (max > min)  {
                while ((x < min) || (x > max))  {
                    if (x < min)    {
                        x = 2*min - x;
                    }
                    if (x > max)    {
                        x = 2*max - x;
                    }
                }
            }
			deltalogprob += HyperLogPrior() + HyperSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
			    x -= m;
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double HyperScalingMove(double& x, double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - HyperLogPrior() - HyperSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
		    x *= e;
			deltalogprob += HyperLogPrior() + HyperSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
			    x /= e;
			}
			ntot++;
		}
		return nacc/ntot;
	}

    void SlaveCollectPathSuffStat() {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectPathSuffStat();
        }
    }

    void SlaveSendNucPathSuffStat()  {

        nucpathsuffstat.Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectComponentPathSuffStat();
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat.Add(geneprocess[gene]->GetNucPathSuffStat());
        }

        int* count = new int[Nnuc+Nnuc*Nnuc];
        double* beta = new double[Nnuc*Nnuc];
        nucpathsuffstat.Push(count,beta);
        MPI_Send(count,Nnuc+Nnuc*Nnuc,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(beta,Nnuc*Nnuc,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] count;
        delete[] beta;
    }

    void MasterReceiveNucPathSuffStat()  {

        int* count = new int[Nnuc+Nnuc*Nnuc];
        double* beta = new double[Nnuc*Nnuc];
        nucpathsuffstat.Clear();
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nnuc+Nnuc*Nnuc,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Nnuc*Nnuc,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            nucpathsuffstat.Add(count,beta);
        }
        delete[] count;
        delete[] beta;
    }

    void MasterMoveNuc()    {

		MoveRR(0.1,1,3);
		MoveRR(0.03,3,3);
		MoveRR(0.01,3,3);

		MoveNucStat(0.1,1,3);
		MoveNucStat(0.01,1,3);
    }

    double NucPathSuffStatLogProb() {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

	double MoveRR(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Nrr];
		for (int rep=0; rep<nrep; rep++)	{
			for (int l=0; l<Nrr; l++)	{
				bk[l] = nucrelrate[l];
			}
			double deltalogprob = -NucPathSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucrelrate,Nrr,tuning,n);
			deltalogprob += loghastings;
            UpdateNucMatrix();
			deltalogprob += NucPathSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nrr; l++)	{
					nucrelrate[l] = bk[l];
				}
                UpdateNucMatrix();
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MoveNucStat(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Nnuc];
		for (int rep=0; rep<nrep; rep++)	{
			for (int l=0; l<Nnuc; l++)	{
				bk[l] = nucstat[l];
			}
			double deltalogprob = -NucPathSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucstat,Nnuc,tuning,n);
			deltalogprob += loghastings;
            UpdateNucMatrix();
			deltalogprob += NucPathSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nnuc; l++)	{
					nucstat[l] = bk[l];
				}
                UpdateNucMatrix();
			}
			ntot++;
		}
		return nacc/ntot;
	}

    void SlaveMoveOmega()  {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->MoveOmega();
            geneprocess[gene]->GetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
        }
    }

    void SlaveSendMixture()   {

        int ngene = GetLocalNgene();
        double* array = new double[4*ngene];
        for (int gene=0; gene<ngene; gene++)    {
            array[gene] = geneprocess[gene]->GetPurOm();
            array[1*ngene+gene] = geneprocess[gene]->GetDPosOm();
            array[2*ngene+gene] = geneprocess[gene]->GetPurW();
            array[3*ngene+gene] = geneprocess[gene]->GetPosW();
        }
        MPI_Send(array,4*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }

    void MasterReceiveMixture()    {

        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            double* array = new double[3*ngene];
            MPI_Recv(array,3*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            int index = 0;
            for (int gene=0; gene<Ngene; gene++)    {
                if (GeneAlloc[gene] == proc)    {
                    (*puromarray)[gene] = array[index];
                    (*dposomarray)[gene] = array[1*ngene+index];
                    (*purwarray)[gene] = array[2*ngene+index];
                    (*poswarray)[gene] = array[3*ngene+index];
                    index++;
                }
            }
            delete[] array;
        }
    }

    void SlaveReceiveMixture()   {

        int ngene = GetLocalNgene();
        double* array = new double[4*ngene];
        MPI_Status stat;
        MPI_Recv(array,4*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

        for (int gene=0; gene<ngene; gene++)    {
            (*puromarray)[gene] = array[gene];
            (*dposomarray)[gene] = array[1*ngene+gene];
            (*purwarray)[gene] = array[2*ngene+gene];
            (*poswarray)[gene] = array[3*ngene+gene];
            geneprocess[gene]->SetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
        }
        delete[] array;
    }

    void MasterSendMixture()    {

        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            double* array = new double[4*ngene];
            int index = 0;
            for (int gene=0; gene<Ngene; gene++)    {
                if (GeneAlloc[gene] == proc)    {
                    array[index] = (*puromarray)[gene];
                    array[1*ngene+index] = (*dposomarray)[gene];
                    array[2*ngene+index] = (*purwarray)[gene];
                    array[3*ngene+index] = (*poswarray)[gene];
                    index++;
                }
            }
            MPI_Send(array,4*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
            delete[] array;
        }
    }

    void SlaveSendLogLikelihood()   {

        totlnL = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            lnL[gene] = geneprocess[gene]->GetLogLikelihood();
            totlnL += lnL[gene];
        }
        MPI_Send(&totlnL,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }

    void MasterReceiveLogLikelihood()    {

        MPI_Status stat;
        totlnL = 0;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            double tmp;
            MPI_Recv(&tmp,1,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            totlnL += tmp;
        }
    }
};

