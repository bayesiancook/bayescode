
#include "CodonM8Model.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDDirichlet.hpp"

class MultiGeneCodonM8Model : public MultiGeneMPIModule	{

    private:

	Tree* tree;
	FileSequenceAlignment* refdata;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

	double aalpha;
    double abeta;
	IIDGamma* alphaarray;
	GammaSuffStat alphasuffstat;

    double balpha;
    double bbeta;
	IIDGamma* betaarray;
	GammaSuffStat betasuffstat;

    double pi;
    double poswalpha;
    double poswbeta;
    IIDBernoulliBeta* poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double dposomalpha;
    double dposombeta;
    IIDGamma* dposomarray;
    GammaSuffStat dposomsuffstat;

    vector<double> purifcenter;
    double purifconcentration;
    IIDDirichlet* purifweightarray;
    DirichletSuffStat purifweightsuffstat;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    NucPathSuffStat nucpathsuffstat;

    std::vector<CodonM8Model*> geneprocess;

    double* lnL;

    int ncat;
    int withpos;

    public:

	CodonStateSpace* GetCodonStateSpace()   {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    MultiGeneCodonM8Model(string datafile, string intreefile, int inncat, int inwithpos, int inmyid, int innprocs) : MultiGeneMPIModule(inmyid,innprocs), purifweightsuffstat(3) {

        ncat = inncat;
        withpos = inwithpos;

        AllocateAlignments(datafile);
        treefile = intreefile;

        // all datafiles have all taxa (with missing data if needed) in same order
        // makes it easier to register tree with data, etc.

        string filename = genename[0];
        refdata = new FileSequenceAlignment(filename);
		refcodondata = new CodonSequenceAlignment(refdata, true);
        taxonset = refdata->GetTaxonSet();
        Ntaxa = refdata->GetNtaxa();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        std::cerr << "number of taxa : " << Ntaxa << '\n';
        std::cerr << "number of branches : " << Nbranch << '\n';
        std::cerr << "-- Tree and data fit together\n";

        Allocate();
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

        aalpha = abeta = 1.0;
        balpha = bbeta = 1.0;
        pi = 0.1;
        poswalpha = poswbeta = 1.0;
        dposomalpha = dposombeta = 1.0;

        purifconcentration = 3.0;
        purifcenter.assign(3,1.0/3);

		alphaarray = new IIDGamma(GetNgene(),aalpha,abeta);
		betaarray = new IIDGamma(GetNgene(),balpha,bbeta);
		dposomarray = new IIDGamma(GetNgene(),dposomalpha,dposombeta);
        poswarray = new IIDBernoulliBeta(GetNgene(),pi,poswalpha,poswbeta);
        purifweightarray = new IIDDirichlet(GetNgene(),purifcenter,purifconcentration);

        lnL = new double[Ngene];

        cerr << "gene processes\n";
        if (! GetMyid())    {
            geneprocess.assign(0,(CodonM8Model*) 0);
        }
        else    {
            geneprocess.assign(Ngene,(CodonM8Model*) 0);

            for (int gene=0; gene<GetNgene(); gene++)   {
                if (genealloc[gene] == myid)    {
                    geneprocess[gene] = new CodonM8Model(genename[gene],treefile,ncat,withpos);
                }
            }
        }
    }

    void Unfold()   {

        if (! GetMyid())    {
            MasterSendGlobalParameters();
            MasterSendMixture();
        }
        else    {
            SlaveReceiveGlobalParameters();
            SlaveReceiveMixture();

            for (int gene=0; gene<GetNgene(); gene++)   {
                if (genealloc[gene] == myid)    {
                    geneprocess[gene]->UpdateMatrices();
                    geneprocess[gene]->Unfold();
                }
            }
        }
    }

    void TraceHeader(ostream& os)   {

        os << "#logprior\tlnL\tlength\t";
        os << "pi\t";
        os << "nposfrac\t";
        os << "meanposfrac\t";
        os << "meanposom\t";
        os << "meanw0\tmeanw1\tpurifconc\t";
        os << "aalpha\tabeta\tbalpha\tbbeta\t";
        os << "poswalpha\tposwbeta\t";
        os << "dposomalpha\tdposombeta\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void MasterTrace(ostream& os)    {
		os << GetLogPrior() << '\t';
        MasterReceiveLogLikelihood();
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
        os << pi << '\t';
        os << GetNpos() << '\t';
        os << GetMeanPosFrac() << '\t';
        os << GetMeanPosOmega() << '\t';
        os << purifcenter[0] << '\t' << purifcenter[2] << '\t' << purifconcentration << '\t';
        os << aalpha << '\t' << abeta << '\t' << balpha << '\t' << bbeta << '\t';
        os << poswalpha << '\t' << poswbeta << '\t';
        os << dposomalpha << '\t' << dposombeta << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
		os.flush();
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
        return tot / totweight;
    }

    void SlaveTrace()   {
        SlaveSendLogLikelihood();
    }

    double GetLogPrior()    {
		double total = 0;
		total += LambdaLogProb();
		total += LengthLogProb();
		total += HyperLogPrior();
		return total;
    }

	double LambdaLogProb()	{
		return -lambda / 10;
	}

	double LengthSuffStatLogProb()	{
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double HyperLogPrior()   {
        return - (aalpha + abeta + balpha + bbeta + poswalpha + poswbeta + dposomalpha + dposombeta + purifconcentration) / 10;
    }

    double HyperSuffStatLogProb()   {
        double total = 0;
        total += alphasuffstat.GetLogProb(aalpha,abeta);
        total += betasuffstat.GetLogProb(balpha,bbeta);
        total += dposomsuffstat.GetLogProb(dposomalpha,dposombeta);
        total += poswsuffstat.GetLogProb(pi,poswalpha,poswbeta);
        total += purifweightsuffstat.GetLogProb(purifcenter,purifconcentration);
        return total;
    }

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

    double GetLogLikelihood()   {
        double tot = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            tot += lnL[gene];
        }
        return tot;
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

            MasterReceiveLengthSuffStat();
            MasterResampleBranchLengths();
            MasterMoveLambda();
            MasterSendGlobalParameters();

            MasterReceiveMixture();
            MasterMoveMixtureHyperParameters();
            MasterSendGlobalParameters();

            MasterReceiveNucPathSuffStat();
            MasterMoveNuc();
            MasterSendGlobalParameters();
        }
    }

    // slave move
    void SlaveMove() {

        SlaveResampleSub();

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalParameters();

            SlaveCollectPathSuffStat();
            SlaveMoveOmega();
            SlaveSendMixture();
            SlaveReceiveGlobalParameters();

            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalParameters();

        }
    }

    void MasterSendGlobalParameters() {

        int N = Nbranch + Nrr + Nnuc + 13;
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
        array[i++] = aalpha;
        array[i++] = abeta;
        array[i++] = balpha;
        array[i++] = bbeta;
        array[i++] = pi;
        array[i++] = poswalpha;
        array[i++] = poswbeta;
        array[i++] = dposomalpha;
        array[i++] = dposombeta;
        array[i++] = purifcenter[0];
        array[i++] = purifcenter[1];
        array[i++] = purifcenter[2];
        array[i++] = purifconcentration;

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveGlobalParameters()   {

        int N = Nbranch + Nrr + Nnuc + 13;
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

        aalpha = array[i++];
        abeta = array[i++];
        balpha = array[i++];
        bbeta = array[i++];
        pi = array[i++];
        poswalpha = array[i++];
        poswbeta = array[i++];
        dposomalpha = array[i++];
        dposombeta = array[i++];
        purifcenter[0] = array[i++];
        purifcenter[1] = array[i++];
        purifcenter[2] = array[i++];
        purifconcentration = array[i++];

        if (i != N) {
            cerr << "error when sending global params: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->SetBranchLengths(*branchlength);
                geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
                geneprocess[gene]->SetMixtureHyperParameters(aalpha,abeta,balpha,bbeta,pi,poswalpha,poswbeta,dposomalpha,dposombeta,purifcenter,purifconcentration);
            }
        }
        delete[] array;
    }

    void SlaveResampleSub()  {

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->ResampleSub(1.0);
            }
        }
    }

    void SlaveSendLengthSuffStat()  {

        lengthsuffstatarray->Clear();
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->CollectLengthSuffStat();
                lengthsuffstatarray->Add(*geneprocess[gene]->GetLengthSuffStatArray());
            }
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

    void MasterMoveMixtureHyperParameters()  {

		alphasuffstat.Clear();
		alphaarray->AddSuffStat(alphasuffstat);
		betasuffstat.Clear();
		betaarray->AddSuffStat(betasuffstat);
		poswsuffstat.Clear();
		poswarray->AddSuffStat(poswsuffstat);
		dposomsuffstat.Clear();
		dposomarray->AddSuffStat(dposomsuffstat);
        purifweightsuffstat.Clear();
        purifweightarray->AddSuffStat(purifweightsuffstat);

		HyperScalingMove(aalpha,1.0,10);
		HyperScalingMove(abeta,1.0,10);
		HyperScalingMove(aalpha,0.3,10);
		HyperScalingMove(abeta,0.3,10);

		HyperScalingMove(balpha,1.0,10);
		HyperScalingMove(bbeta,1.0,10);
		HyperScalingMove(balpha,0.3,10);
		HyperScalingMove(bbeta,0.3,10);

        ResamplePi();
		HyperScalingMove(poswalpha,1.0,10);
		HyperScalingMove(poswbeta,1.0,10);
		HyperScalingMove(poswalpha,0.3,10);
		HyperScalingMove(poswbeta,0.3,10);

		HyperScalingMove(dposomalpha,1.0,10);
		HyperScalingMove(dposombeta,1.0,10);
		HyperScalingMove(dposomalpha,0.3,10);
		HyperScalingMove(dposombeta,0.3,10);

        HyperScalingMove(purifconcentration,1.0,10);
        HyperScalingMove(purifconcentration,0.3,10);

        HyperProfileMove(purifcenter,1.0,1,10);
        HyperProfileMove(purifcenter,0.3,1,10);

        alphaarray->SetShape(aalpha);
        alphaarray->SetScale(abeta);
        betaarray->SetShape(balpha);
        betaarray->SetScale(bbeta);
        poswarray->SetPi(pi);
        poswarray->SetAlpha(poswalpha);
        poswarray->SetBeta(poswbeta);
        dposomarray->SetShape(dposomalpha);
        dposomarray->SetScale(dposombeta);
        purifweightarray->SetConcentration(purifconcentration);
        purifweightarray->SetCenter(purifcenter);
    }

    void ResamplePi()   {

        int n0 = poswsuffstat.GetN0();
        int n1 = poswsuffstat.GetN1();
        if ((n0+n1) != Ngene)   {
            cerr << "error in resample pi\n";
            exit(1);
        }
        double a0 = Random::sGamma(1.0 + n0);
        double a1 = Random::sGamma(1.0 + n1);
        pi = a1 / (a0 + a1);
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
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->CollectPathSuffStat();
            }
        }
    }

    void SlaveSendNucPathSuffStat()  {

        nucpathsuffstat.Clear();
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->CollectComponentPathSuffStat();
                geneprocess[gene]->CollectNucPathSuffStat();
                nucpathsuffstat.Add(geneprocess[gene]->GetNucPathSuffStat());
            }
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

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->MoveOmega();
            }
        }
    }

    void SlaveSendMixture()   {

        // alpha,beta,posw,dposom
        double* array = new double[7*Ngene];
        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                array[gene] = geneprocess[gene]->GetAlpha();
                array[Ngene+gene] = geneprocess[gene]->GetBeta();
                array[2*Ngene+gene] = geneprocess[gene]->GetPosW();
                array[3*Ngene+gene] = geneprocess[gene]->GetDPosOm();
                array[4*Ngene+gene] = geneprocess[gene]->GetPurifWeight(0);
                array[5*Ngene+gene] = geneprocess[gene]->GetPurifWeight(1);
                array[6*Ngene+gene] = geneprocess[gene]->GetPurifWeight(2);
            }
            else    {
                array[gene] = -1;
                array[Ngene+gene] = -1;
                array[2*Ngene+gene] = -1;
                array[3*Ngene+gene] = -1;
                array[4*Ngene+gene] = -1;
                array[5*Ngene+gene] = -1;
                array[6*Ngene+gene] = -1;
            }
        }
        MPI_Send(array,7*Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }

    void MasterReceiveMixture()    {

        double* array = new double[7*Ngene];
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(array,7*Ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            for (int gene=0; gene<Ngene; gene++)    {
                if (array[gene] != -1)    {
                    (*alphaarray)[gene] = array[gene];
                    (*betaarray)[gene] = array[Ngene+gene];
                    (*poswarray)[gene] = array[2*Ngene+gene];
                    (*dposomarray)[gene] = array[3*Ngene+gene];
                    (*purifweightarray)[gene][0] = array[4*Ngene+gene];
                    (*purifweightarray)[gene][1] = array[5*Ngene+gene];
                    (*purifweightarray)[gene][2] = array[6*Ngene+gene];
                }
            }
        }
        delete[] array;
    }

    void MasterSendMixture()  {

        double* array = new double[7*Ngene];
        for (int gene=0; gene<Ngene; gene++)    {
            array[gene] = (*alphaarray)[gene];
            array[Ngene+gene] = (*betaarray)[gene];
            array[2*Ngene+gene] = (*poswarray)[gene];
            array[3*Ngene+gene] = (*dposomarray)[gene];
            array[4*Ngene+gene] = (*purifweightarray)[gene][0];
            array[5*Ngene+gene] = (*purifweightarray)[gene][1];
            array[6*Ngene+gene] = (*purifweightarray)[gene][2];
        }
        MPI_Bcast(array,7*Ngene,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveMixture()    {

        double* array = new double[7*Ngene];
        MPI_Bcast(array,7*Ngene,MPI_DOUBLE,0,MPI_COMM_WORLD);

        for (int gene=0; gene<Ngene; gene++)    {
            (*alphaarray)[gene] = array[gene];
            (*betaarray)[gene] = array[Ngene+gene];
            (*poswarray)[gene] = array[2*Ngene+gene];
            (*dposomarray)[gene] = array[3*Ngene+gene];
            (*purifweightarray)[gene][0] = array[4*Ngene+gene];
            (*purifweightarray)[gene][1] = array[5*Ngene+gene];
            (*purifweightarray)[gene][2] = array[6*Ngene+gene];
        }

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                geneprocess[gene]->SetMixtureParameters((*alphaarray)[gene],(*betaarray)[gene],(*poswarray)[gene],(*dposomarray)[gene],(*purifweightarray)[gene]);
            }
        }
        delete[] array;
    }

    void SlaveSendLogLikelihood()   {

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene] == myid)    {
                lnL[gene] = geneprocess[gene]->GetLogLikelihood();
            }
            else    {
                lnL[gene] = 0;
            }
        }
        MPI_Send(lnL,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }

    void MasterReceiveLogLikelihood()    {

        double* array = new double[Ngene];
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(array,Ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            for (int gene=0; gene<Ngene; gene++)    {
                if (array[gene])    {
                    lnL[gene] = array[gene];
                }
            }
        }
        delete[] array;
    }

};

