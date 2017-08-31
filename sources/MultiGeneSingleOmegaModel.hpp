
#include "SingleOmegaModel.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"

class MultiGeneSingleOmegaModel : public MultiGeneMPIModule	{

    private:

	Tree* tree;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	double alpha;
	double beta;
	IIDGamma* omegaarray;
	GammaSuffStat alphabetasuffstat;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

    NucPathSuffStat nucpathsuffstat;
    std::vector<SingleOmegaModel*> geneprocess;

    double totlnL;
    double* lnL;

    public:

	CodonStateSpace* GetCodonStateSpace()   {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    MultiGeneSingleOmegaModel(string datafile, string intreefile, int inmyid, int innprocs) : MultiGeneMPIModule(inmyid,innprocs) {

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

        alpha = beta = 1.0;
		omegaarray = new IIDGamma(GetLocalNgene(),alpha,beta);

        if (myid)   {
            lnL = new double[GetLocalNgene()];
        }
        else    {
            lnL = 0;
        }

        cerr << "gene processes\n";
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
            MasterReceiveLogLikelihood();
        }
        else    {
            SlaveReceiveGlobalBranchLengths();
            SlaveReceiveGlobalNucRates();
            SlaveReceiveOmegaHyperParameters();
            SlaveReceiveOmega();

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->UpdateMatrices();
                geneprocess[gene]->Unfold();
            }
            SlaveSendLogLikelihood();
        }
    }


    void TraceHeader(ostream& os)   {

        os << "#logprior\tlnL\tlength\t";
        os << "meanomega\t";
        os << "varomega\t";
        os << "alpha\tbeta\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream& os)    {
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
        os << omegaarray->GetMean() << '\t';
        os << omegaarray->GetVar() << '\t';
        os << alpha << '\t' << beta << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
		os.flush();
    }

    double GetLogPrior()    {
		double total = 0;
		total += LambdaLogProb();
		total += LengthLogProb();
		total += OmegaLogProb();
        total += AlphaLogProb();
        total += BetaLogProb();
		return total;
    }

	double LambdaLogProb()	{
		return -lambda / 10;
	}

	double LengthSuffStatLogProb()	{
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double OmegaSuffStatLogProb()   {
        return alphabetasuffstat.GetLogProb(alpha,beta);
    }

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

    double AlphaLogProb()   {
        return -alpha/10;
    }

    double BetaLogProb()    {
        return -beta/10;
    }

    double OmegaLogProb()   {
        return omegaarray->GetLogProb();
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

            MasterReceiveLengthSuffStat();
            MasterResampleBranchLengths();
            MasterMoveLambda();
            MasterSendGlobalBranchLengths();

            MasterReceiveOmega();
            MasterMoveAlphaBeta();
            MasterSendOmegaHyperParameters();

            MasterReceiveNucPathSuffStat();
            MasterMoveNuc();
            MasterSendGlobalNucRates();
        }

        MasterReceiveOmega();
        MasterReceiveLogLikelihood();
    }

    // slave move
    void SlaveMove() {

        SlaveResampleSub();

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalBranchLengths();

            SlaveCollectPathSuffStat();
            SlaveMoveOmega();
            SlaveSendOmega();
            SlaveReceiveOmegaHyperParameters();

            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalNucRates();
        }

        SlaveSendOmega();
        SlaveSendLogLikelihood();
    }

    void SlaveResampleSub()  {

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleSub();
        }
    }

    void MasterSendGlobalBranchLengths() {

        double* array = new double[Nbranch];
        for (int j=0; j<Nbranch; j++)   {
            array[j] = branchlength->GetVal(j);
        }
        MPI_Bcast(array,Nbranch,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveGlobalBranchLengths()   {

        double* array = new double[Nbranch];
        MPI_Bcast(array,Nbranch,MPI_DOUBLE,0,MPI_COMM_WORLD);
        for (int j=0; j<Nbranch; j++)   {
            (*branchlength)[j] = array[j];
        }
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetBranchLengths(*branchlength);
        }
        delete[] array;
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

    void MasterMoveAlphaBeta()  {

		alphabetasuffstat.Clear();
		omegaarray->AddSuffStat(alphabetasuffstat);
		MoveAlpha(1.0,10);
		MoveBeta(1.0,10);
		MoveAlpha(0.3,10);
		MoveBeta(0.3,10);
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

	double MoveAlpha(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - AlphaLogProb() - OmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			alpha *= e;
			deltalogprob += AlphaLogProb() + OmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				alpha /= e;
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MoveBeta(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - BetaLogProb() - OmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			beta *= e;
			deltalogprob += BetaLogProb() + OmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				beta /= e;
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


    void MasterSendGlobalNucRates()   {

        int N = Nrr + Nnuc;
        double* array = new double[N];
        int i = 0;
        for (int j=0; j<Nrr; j++)   {
            array[i++] = nucrelrate[j];
        }
        for (int j=0; j<Nnuc; j++)  {
            array[i++] = nucstat[j];
        }
        if (i != N) {
            cerr << "error when sending global nuc rates: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveGlobalNucRates()   {

        int N = Nrr + Nnuc;
        double* array = new double[N];
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        int i = 0;
        for (int j=0; j<Nrr; j++)   {
            nucrelrate[j] = array[i++];
        }
        for (int j=0; j<Nnuc; j++)  {
            nucstat[j] = array[i++];
        }

        if (i != N) {
            cerr << "error when receiving global nuc rates: non matching vector size\n";
            cerr << i << '\t' << N << '\n';
            exit(1);
        }

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
        }
        delete[] array;
    }

    void SlaveMoveOmega()  {

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->MoveOmega();
            (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
        }
    }

    void SlaveSendOmega()   {

        int ngene = GetLocalNgene();
        double* array = new double[ngene];
        for (int gene=0; gene<ngene; gene++)   {
            array[gene] = (*omegaarray)[gene];
        }
        MPI_Send(array,ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }

    void MasterReceiveOmega()    {

        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            double* array = new double[ngene];
            MPI_Recv(array,ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            int index = 0;
            for (int gene=0; gene<Ngene; gene++)    {
                if (GeneAlloc[gene] == proc)    {
                    (*omegaarray)[gene] = array[index++];
                }
            }
            if (index != ngene) {
                cerr << "error: non matching number of genes in master receive omega\n";
                exit(1);
            }
            delete[] array;
        }
    }

    void MasterSendOmega()  {

        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            double* array = new double[4*ngene];
            int index = 0;
            for (int gene=0; gene<Ngene; gene++)    {
                if (GeneAlloc[gene] == proc)    {
                    array[index++] = (*omegaarray)[gene];
                }
            }
            if (index != ngene) {
                cerr << "error: non matching number of genes in master send omega\n";
                exit(1);
            }
            MPI_Send(array,ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
            delete[] array;
        }
    }

    void SlaveReceiveOmega()    {

        int ngene = GetLocalNgene();
        double* array = new double[ngene];
        MPI_Status stat;
        MPI_Recv(array,ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

        for (int gene=0; gene<ngene; gene++)    {
            (*omegaarray)[gene] = array[gene];
            geneprocess[gene]->SetOmega(array[gene]);
        }
        delete[] array;
    }

    void MasterSendOmegaHyperParameters()   {

        double* array = new double[2];
        array[0] = alpha;
        array[1] = beta;
        MPI_Bcast(array,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveOmegaHyperParameters() {

        double* array = new double[Ngene];
        MPI_Bcast(array,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
        alpha = array[0];
        beta = array[1];
        delete[] array;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
        for (int gene=0; gene<GetLocalNgene(); gene++)    {
            geneprocess[gene]->SetAlphaBeta(alpha,beta);
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

