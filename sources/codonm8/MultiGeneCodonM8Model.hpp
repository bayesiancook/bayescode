
#include "CodonM8Model.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "IIDDirichlet.hpp"

#include "Chrono.hpp"


class MultiGeneCodonM8Model : public MultiGeneMPIModule	{

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

    double purifmeanhypermean;
	double purifmeanhyperinvconc;
	IIDBeta* purifmeanarray;
	BetaSuffStat purifmeansuffstat;

    double purifinvconchypermean;
    double purifinvconchyperinvshape;
	IIDGamma* purifinvconcarray;
	GammaSuffStat purifinvconcsuffstat;

    double dposomhypermean;
    double dposomhyperinvshape;
    IIDGamma* dposomarray;
    GammaSuffStat dposomsuffstat;

    double poswhypermean;
    double poswhyperinvconc;
    IIDBernoulliBeta* poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    vector<double> purifweighthypercenter;
    double purifweighthyperinvconc;
    IIDDirichlet* purifweightarray;
    DirichletSuffStat purifweightsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double pi;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    NucPathSuffStat nucpathsuffstat;

    std::vector<CodonM8Model*> geneprocess;

    double totlnL;
    double* lnL;

    int ncat;

    int burnin;

    Chrono timepercycle;
    Chrono omegachrono,hyperchrono,mastersampling;

    public:

	CodonStateSpace* GetCodonStateSpace()   {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    MultiGeneCodonM8Model(string datafile, string intreefile, int inncat, double inpihypermean, double inpihyperinvconc, int inmyid, int innprocs) : MultiGeneMPIModule(inmyid,innprocs), purifweightsuffstat(3) {

        burnin = 0;

        ncat = inncat;
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

        purifmeanhypermean = 0.5;
        purifmeanhyperinvconc = 0.5;
        purifinvconchypermean = 2.0;
        purifinvconchyperinvshape = 1.0;

        dposomhypermean = 1.0;
        dposomhyperinvshape = 1.0;

        pi = pihypermean;
        if (! pi)   {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }
        else    {
            poswhypermean = 0.1;
            poswhyperinvconc = 1;
        }

        purifweighthypercenter.assign(3,1.0/3);
        purifweighthyperinvconc = 1.0/3;

        double purifmeanalpha = purifmeanhypermean / purifmeanhyperinvconc;
        double purifmeanbeta = (1-purifmeanhypermean) / purifmeanhyperinvconc;
        purifmeanarray = new IIDBeta(GetLocalNgene(),purifmeanalpha,purifmeanbeta);

        double purifinvconcalpha = 1.0 / purifinvconchyperinvshape;
        double purifinvconcbeta = purifinvconcalpha / purifinvconchypermean;
        purifinvconcarray = new IIDGamma(GetLocalNgene(),purifinvconcalpha,purifinvconcbeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray = new IIDGamma(GetLocalNgene(),dposomalpha,dposombeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        poswarray = new IIDBernoulliBeta(GetLocalNgene(),pi,poswalpha,poswbeta);

        double purifweighthyperconc = 1.0 / purifweighthyperinvconc;
        purifweightarray = new IIDDirichlet(GetLocalNgene(),purifweighthypercenter,purifweighthyperconc);

        if (myid)   {
            lnL = new double[GetLocalNgene()];
        }
        else    {
            lnL = 0;
        }


        if (! GetMyid())    {
            geneprocess.assign(0,(CodonM8Model*) 0);
        }
        else    {
            geneprocess.assign(GetLocalNgene(),(CodonM8Model*) 0);

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene] = new CodonM8Model(GetLocalGeneName(gene),treefile,ncat,pi);
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
            // SlaveSetArrays();
            SlaveReceiveMixture();

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                // geneprocess[gene]->SetMixtureParameters((*purifmeanarray)[gene],(*purifinvconcarray)[gene],(*poswarray)[gene],(*dposomarray)[gene],(*purifweightarray)[gene]);
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
        os << "meanposfrac\t";
        os << "meanposom\t";
        os << "meanw0\tmeanw1\tpurifinvconc\t";
        os << "purifmean\tinvconc\tpurifinvconc\tinvshape\t";
        os << "poswmean\tinvconc\t";
        os << "dposommean\tinvshape\t";
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
        os << GetMeanPosFrac() << '\t';
        os << GetMeanPosOmega() << '\t';
        os << purifweighthypercenter[0] << '\t' << purifweighthypercenter[2] << '\t' << purifweighthyperinvconc << '\t';
        os << purifmeanhypermean << '\t' << purifmeanhyperinvconc << '\t' << purifinvconchypermean << '\t' << purifinvconchyperinvshape << '\t';
        os << poswhypermean << '\t' << poswhyperinvconc << '\t';
        os << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
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
        os << purifmeanhypermean << '\t';
        os << purifmeanhyperinvconc << '\t';
        os << purifinvconchypermean << '\t';
        os << purifinvconchyperinvshape << '\t';
        os << pi << '\t';
        os << poswhypermean << '\t';
        os << poswhyperinvconc << '\t';
        os << dposomhypermean << '\t';
        os << dposomhyperinvshape << '\t';
        os << purifweighthypercenter[0] << '\t';
        os << purifweighthypercenter[1] << '\t';
        os << purifweighthypercenter[2] << '\t';
        os << purifweighthyperinvconc;
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

    double GetMeanPosFrac() {
        return poswarray->GetPosMean();
    }

    double GetMeanPosOmega()    {

        double tot = 0;
        double totweight = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (poswarray->GetVal(gene)) {
                tot += poswarray->GetVal(gene) * dposomarray->GetVal(gene);
                totweight += poswarray->GetVal(gene);
            }
        }
        if (! totweight)    {
            return 1;
        }
        return tot / totweight + 1;
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
        total -= purifmeanhyperinvconc;
        total -= purifinvconchypermean;
        total -= purifinvconchyperinvshape;
        total -= 10*poswhyperinvconc;
        total -= dposomhypermean;
        total -= 10*dposomhyperinvshape;
        total -= purifweighthyperinvconc;
        return total;
    }

    double HyperSuffStatLogProb()   {
        double total = 0;

        double purifmeanalpha = purifmeanhypermean / purifmeanhyperinvconc;
        double purifmeanbeta = (1-purifmeanhypermean) / purifmeanhyperinvconc;
        total += purifmeansuffstat.GetLogProb(purifmeanalpha,purifmeanbeta);

        double purifinvconcalpha = 1.0 / purifinvconchyperinvshape;
        double purifinvconcbeta = purifinvconcalpha / purifinvconchypermean;
        total += purifinvconcsuffstat.GetLogProb(purifinvconcalpha,purifinvconcbeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        total += dposomsuffstat.GetLogProb(dposomalpha,dposombeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        total += poswsuffstat.GetLogProb(pi,poswalpha,poswbeta);

        double purifweighthyperconc = 1.0 / purifweighthyperinvconc;
        total += purifweightsuffstat.GetLogProb(purifweighthypercenter,purifweighthyperconc);

        if (isnan(total))   {
            cerr << "hyper suff stat log prob is nan\n";
            exit(1);
        }

        return total;
    }

    void PrintHyperSuffStat(ostream& os)   {

            os << "purifmean      : " << purifmeansuffstat.GetN() << '\t' << purifmeansuffstat.GetSumLog0() << '\t' << purifmeansuffstat.GetSumLog1() << '\n';
            os << "purifinvconc   : " << purifinvconcsuffstat.GetN() << '\t' << purifinvconcsuffstat.GetSum() << '\t' << purifinvconcsuffstat.GetSumLog() << '\n';
            os << "posw           : " << poswsuffstat.GetN0() << '\t' << poswsuffstat.GetN1() << '\t' << poswsuffstat.GetSumLog0() << '\t' << poswsuffstat.GetSumLog1() << '\n';
            os << "dposom         : " << dposomsuffstat.GetN() << '\t' << dposomsuffstat.GetSum() << '\t' << dposomsuffstat.GetSumLog() << '\n';
            os << "purifweight    : " << purifweightsuffstat.GetN() << '\t' << purifweightsuffstat.GetSumLog(0) << '\t' << purifweightsuffstat.GetSumLog(1) << '\t' << purifweightsuffstat.GetSumLog(2) << '\n';
    }

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

    double GetLogLikelihood()   {
        return totlnL;
        /*
        double tot = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            tot += lnL[gene];
        }
        return tot;
        */
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
            // MasterReceiveMixture();
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
            // SlaveSendMixture();
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

        int N = 13;
        double* array = new double[N];
        int i = 0;
        array[i++] = purifmeanhypermean;
        array[i++] = purifmeanhyperinvconc;
        array[i++] = purifinvconchypermean;
        array[i++] = purifinvconchyperinvshape;
        array[i++] = pi;
        array[i++] = poswhypermean;
        array[i++] = poswhyperinvconc;
        array[i++] = dposomhypermean;
        array[i++] = dposomhyperinvshape;
        array[i++] = purifweighthypercenter[0];
        array[i++] = purifweighthypercenter[1];
        array[i++] = purifweighthypercenter[2];
        array[i++] = purifweighthyperinvconc;

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveHyperParameters()   {

        int N = 13;
        double* array = new double[N];
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        int i = 0;
        purifmeanhypermean = array[i++];
        purifmeanhyperinvconc = array[i++];
        purifinvconchypermean = array[i++];
        purifinvconchyperinvshape = array[i++];
        pi = array[i++];
        poswhypermean = array[i++];
        poswhyperinvconc = array[i++];
        dposomhypermean = array[i++];
        dposomhyperinvshape = array[i++];
        purifweighthypercenter[0] = array[i++];
        purifweighthypercenter[1] = array[i++];
        purifweighthypercenter[2] = array[i++];
        purifweighthyperinvconc = array[i++];

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetMixtureHyperParameters(purifmeanhypermean,purifmeanhyperinvconc,purifinvconchypermean,purifinvconchyperinvshape,dposomhypermean,dposomhyperinvshape,pi,poswhypermean,poswhyperinvconc,purifweighthypercenter,purifweighthyperinvconc);
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

        int Nint = 1 + 1 + 2 + 1 + 1;
        int Ndouble = 2 + 2 + 2 + 2 + 3;
        int* count = new int[Nint];
        double* beta = new double[Ndouble];

        int i = 0;
        int d = 0;

		purifmeansuffstat.Clear();
		purifmeanarray->AddSuffStat(purifmeansuffstat);
        count[i++] = purifmeansuffstat.GetN();
        beta[d++] = purifmeansuffstat.GetSumLog0();
        beta[d++] = purifmeansuffstat.GetSumLog1();

		purifinvconcsuffstat.Clear();
		purifinvconcarray->AddSuffStat(purifinvconcsuffstat);
        count[i++] = purifinvconcsuffstat.GetN();
        beta[d++] = purifinvconcsuffstat.GetSum();
        beta[d++] = purifinvconcsuffstat.GetSumLog();

		poswsuffstat.Clear();
		poswarray->AddSuffStat(poswsuffstat);
        count[i++] = poswsuffstat.GetN0();
        count[i++] = poswsuffstat.GetN1();
        beta[d++] = poswsuffstat.GetSumLog0();
        beta[d++] = poswsuffstat.GetSumLog1();

		dposomsuffstat.Clear();
		dposomarray->AddSuffStat(dposomsuffstat);
		// dposomarray->AddSuffStat(dposomsuffstat,*poswarray);
        count[i++] = dposomsuffstat.GetN();
        beta[d++] = dposomsuffstat.GetSum();
        beta[d++] = dposomsuffstat.GetSumLog();

        purifweightsuffstat.Clear();
        purifweightarray->AddSuffStat(purifweightsuffstat);
        count[i++] = purifweightsuffstat.GetN();
        beta[d++] = purifweightsuffstat.GetSumLog(0);
        beta[d++] = purifweightsuffstat.GetSumLog(1);
        beta[d++] = purifweightsuffstat.GetSumLog(2);

        MPI_Send(count,Nint,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(beta,Ndouble,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

        delete[] count;
        delete[] beta;
    }

    void MasterReceiveMixtureHyperSuffStat()  {

		purifmeansuffstat.Clear();
		purifinvconcsuffstat.Clear();
		poswsuffstat.Clear();
		dposomsuffstat.Clear();
        purifweightsuffstat.Clear();

        int Nint = 1 + 1 + 2 + 1 + 1;
        int Ndouble = 2 + 2 + 2 + 2 + 3;
        int* count = new int[Nint];
        double* beta = new double[Ndouble];

        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nint,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Ndouble,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

            int i = 0;
            int d = 0;

            purifmeansuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
            purifinvconcsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
            poswsuffstat.AddNullSuffStat(count[i++]);
            poswsuffstat.AddPosSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
            dposomsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
            d+=2;
            purifweightsuffstat.AddSuffStat(beta+d,count[i++]);
            d += 3;
        }

        delete[] count;
        delete[] beta;
    }

    void MasterMoveMixtureHyperParameters()  {

		HyperSlidingMove(purifmeanhypermean,1.0,10,0,1);
		HyperSlidingMove(purifmeanhypermean,0.3,10,0,1);
		HyperScalingMove(purifmeanhyperinvconc,1.0,10);
		HyperScalingMove(purifmeanhyperinvconc,0.3,10);

		HyperScalingMove(purifinvconchypermean,1.0,10);
		HyperScalingMove(purifinvconchypermean,0.3,10);
		HyperScalingMove(purifinvconchyperinvshape,1.0,10);
		HyperScalingMove(purifinvconchyperinvshape,0.3,10);

		HyperSlidingMove(poswhypermean,1.0,10,0,1);
		HyperSlidingMove(poswhypermean,0.3,10,0,1);
		HyperScalingMove(poswhyperinvconc,1.0,10);
		HyperScalingMove(poswhyperinvconc,0.3,10);

		HyperSlidingMove(dposomhypermean,3.0,10,0.5,10);
		HyperSlidingMove(dposomhypermean,1.0,10,0.5,10);
		HyperSlidingMove(dposomhypermean,0.3,10,0.5,10);
		HyperScalingMove(dposomhyperinvshape,1.0,10);
		HyperScalingMove(dposomhyperinvshape,0.3,10);

        HyperProfileMove(purifweighthypercenter,1.0,1,10);
        HyperProfileMove(purifweighthypercenter,0.3,1,10);
        HyperScalingMove(purifweighthyperinvconc,1.0,10);
        HyperScalingMove(purifweighthyperinvconc,0.3,10);

        if (burnin > 10)    {
            if (pihyperinvconc)    {
                ResamplePi();
            }
        }
    }

    void SlaveSetArrays()    {

        double purifmeanalpha = purifmeanhypermean / purifmeanhyperinvconc;
        double purifmeanbeta = (1-purifmeanhypermean) / purifmeanhyperinvconc;
        purifmeanarray->SetAlpha(purifmeanalpha);
        purifmeanarray->SetBeta(purifmeanbeta);

        double purifinvconcalpha = 1.0 / purifinvconchyperinvshape;
        double purifinvconcbeta = purifinvconcalpha / purifinvconchypermean;
        purifinvconcarray->SetShape(purifinvconcalpha);
        purifinvconcarray->SetScale(purifinvconcbeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray->SetShape(dposomalpha);
        dposomarray->SetScale(dposombeta);
        // dposomarray->PriorResample(*poswarray);
        // necessary after changing some dposom values
        /*
        for (int gene=0; gene<GetLocalNgene(); gene++)    {
            geneprocess[gene]->SetMixtureParameters((*purifmeanarray)[gene],(*purifinvconcarray)[gene],(*poswarray)[gene],(*dposomarray)[gene],(*purifweightarray)[gene]);
        }
        */

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1-poswhypermean) / poswhyperinvconc;
        poswarray->SetPi(pi);
        poswarray->SetAlpha(poswalpha);
        poswarray->SetBeta(poswbeta);

        double purifweighthyperconc = 1.0 / purifweighthyperinvconc;
        purifweightarray->SetCenter(purifweighthypercenter);
        purifweightarray->SetConcentration(purifweighthyperconc);
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

	double HyperProfileMove(vector<double>& x, double tuning, int n, int nrep)	{

		double nacc = 0;
		double ntot = 0;
        vector<double> bk(x.size(),0);
		for (int rep=0; rep<nrep; rep++)	{
            bk = x;
			double deltalogprob = - HyperLogPrior() - HyperSuffStatLogProb();
            double loghastings = Random::ProfileProposeMove(x,x.size(),tuning,n);
			deltalogprob += HyperLogPrior() + HyperSuffStatLogProb();
			deltalogprob += loghastings;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
			    x = bk;
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
            geneprocess[gene]->GetMixtureParameters((*purifmeanarray)[gene],(*purifinvconcarray)[gene],(*poswarray)[gene],(*dposomarray)[gene],(*purifweightarray)[gene]);
        }
    }

    void SlaveSendMixture()   {

        int ngene = GetLocalNgene();
        double* array = new double[7*ngene];
        for (int gene=0; gene<ngene; gene++)    {
            array[gene] = geneprocess[gene]->GetPurifMean();
            array[ngene+gene] = geneprocess[gene]->GetPurifInvConc();
            array[2*ngene+gene] = geneprocess[gene]->GetPosW();
            array[3*ngene+gene] = geneprocess[gene]->GetDPosOm();
            array[4*ngene+gene] = geneprocess[gene]->GetPurifWeight(0);
            array[5*ngene+gene] = geneprocess[gene]->GetPurifWeight(1);
            array[6*ngene+gene] = geneprocess[gene]->GetPurifWeight(2);
        }
        MPI_Send(array,7*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }

    void MasterReceiveMixture()    {

        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            double* array = new double[7*ngene];
            MPI_Recv(array,7*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            int index = 0;
            for (int gene=0; gene<Ngene; gene++)    {
                if (GeneAlloc[gene] == proc)    {
                    (*purifmeanarray)[gene] = array[index];
                    (*purifinvconcarray)[gene] = array[ngene+index];
                    (*poswarray)[gene] = array[2*ngene+index];
                    (*dposomarray)[gene] = array[3*ngene+index];
                    (*purifweightarray)[gene][0] = array[4*ngene+index];
                    (*purifweightarray)[gene][1] = array[5*ngene+index];
                    (*purifweightarray)[gene][2] = array[6*ngene+index];
                    index++;
                }
            }
            delete[] array;
        }
    }

    void SlaveReceiveMixture()   {

        int ngene = GetLocalNgene();
        double* array = new double[7*ngene];
        MPI_Status stat;
        MPI_Recv(array,7*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

        for (int gene=0; gene<ngene; gene++)    {
            (*purifmeanarray)[gene] = array[gene];
            (*purifinvconcarray)[gene] = array[ngene+gene];
            (*poswarray)[gene] = array[2*ngene+gene];
            (*dposomarray)[gene] = array[3*ngene+gene];
            (*purifweightarray)[gene][0] = array[4*ngene+gene];
            (*purifweightarray)[gene][1] = array[5*ngene+gene];
            (*purifweightarray)[gene][2] = array[6*ngene+gene];
            geneprocess[gene]->SetMixtureParameters((*purifmeanarray)[gene],(*purifinvconcarray)[gene],(*poswarray)[gene],(*dposomarray)[gene],(*purifweightarray)[gene]);
        }
        delete[] array;
    }

    void MasterSendMixture()    {

        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            double* array = new double[7*ngene];
            int index = 0;
            for (int gene=0; gene<Ngene; gene++)    {
                if (GeneAlloc[gene] == proc)    {
                    array[index] = (*purifmeanarray)[gene];
                    array[ngene+index] = (*purifinvconcarray)[gene];
                    array[2*ngene+index] = (*poswarray)[gene];
                    array[3*ngene+index] = (*dposomarray)[gene];
                    array[4*ngene+index] = (*purifweightarray)[gene][0];
                    array[5*ngene+index] = (*purifweightarray)[gene][1];
                    array[6*ngene+index] = (*purifweightarray)[gene][2];
                    index++;
                }
            }
            MPI_Send(array,7*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
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

