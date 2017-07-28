/*
simple model: gene-specific omegas, iid from a gamma distribution
shared nuc rates and branch lengths across all genes
each slave allocates a small number of gene-specific models
master deals with branch lengths, nuc stat and shape and scale params of gamma distribution of omega's across genes

starting:
master broadcasts global params / slaves receive global params and update phyloprocess

move schedule:

slave resample sub

for i=1..Nrep

    slaves compile branch length suff stats and send them to master
    master collects branch length suff stat, resamples branch lengths (and hyperparams) and broadcasts new branch lengths

    slaves compile pathsuffstat, move omega and send to master 
    master moves alpha beta sends back to slaves

    slave compile nuc suff stats and send them to master
    master collects nuc suffstat across genes, resamples nuc rates and broadcasts them to slaves


- use SingleOmegaModel as a single class, usable under all conditions (with flags indicating which aspects are dependent / independent)
*/

#include "SingleOmegaModel.hpp"
#include "Parallel.hpp"

class MultiGeneSingleOmegaModel : public MultiGeneMPIModule	{

    private:

	Tree* tree;
	FileSequenceAlignment* refdata;
	const TaxonSet* taxonset;
	// CodonSequenceAlignment* codondata;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	double alpha;
	double beta;
	IIDGamma* omegaarray;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;

    NucPathSuffStat nucpathsuffstat;
    std::vector<SingleOmegaModel*> geneprocess;

    double* lnL;

    public:

    MultiGeneSingleOmegaModel(string datafile, string treefile, int inmyid, int innprocs) : MultiGeneMPIModule(inmyid,innprocs) {

        AllocateAlignments(datafile);

        // all datafiles have all taxa (with missing data if needed) in same order
        // makes it easier to register tree with data, etc.

        string datafile = genename[0];
        refdata = new FileSequenceAlignment(datafile);
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

        alpha = beta = 1.0;
		omegaarray = new IIDGamma(GetNgene(),alpha,beta);

        lnL = new double[Ngene];

        if (! GetMyid())    {
            geneprocess.assign(0,0);

            MasterSendGlobalParameters();
            MasterSendOmega();
        }
        else    {
            geneprocess.assign(Ngene,0);

            SlaveReceiveGlobalParameters();
            SlaveReceiveOmega();

            for (int gene=0; gene<GetNgene(); gene++)   {
                if (genealloc[gene])    {
                    geneprocess[gene] = new SingleOmegaModel(genename[gene],treefile,1);
                    geneprocess[gene]->SetBranchLengths(branchlength);
                    geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
                    geneprocess[gene]->SetAlphaBeta(alpha,beta);
                    geneprocess[gene]->SetOmega(omegarray->GetVal(gene));
                    geneprocess[gene]->UpdateMatrices();
                    geneprocess[gene]->Unfold();
                }
            }
        }
        Trace(cerr);
    }

    void Trace(ostream& os) {
        if (! GetMyid())    {
            MasterTrace(ostream& os);
        }
        else    {
            SlaveTrace();
        }
    }

    void TraceHeader(ostream& os)   {

        os << "#logprior\tlnL\tlength\t";
        os << "meanomega\t";
        os << "varomega\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void MasterTrace(ostream& os)    {
		os << GetLogPrior() << '\t';
        MasterReceiveLogLikelihood();
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
        os << omegaarray->GetMean() << '\t'
        os << omegaarray->GetVar() << '\t'
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
		os.flush();
    }

    void SlaveTrace()   {
        SlaveSendLogLikelihood();
    }

    double GetLogPrior()    {
		double total = 0;
		total += LambdaLogProb();
		total += LengthLogProb();
		total += OmegaLogProb();
        total += AlphaBetaLogProb();
		return total;
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

    void Move() {

        if (! GetMyid())    {
            MasterMove();
        }
        else    {
            SlaveMove();
        }
    }


    void MasterMove() {

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            MasterReceiveLengthSuffStat();
            MasterResampleBranchLengths();
            MasterMoveLambda();
            MasterSendGlobalParameters();

            MasterReceiveOmega();
            MasterMoveAlphaBeta();
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
            SlaveSendOmega();
            SlaveReceiveGlobalParameters();

            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalParameters();

        }
    }

    void MasterSendGlobalParameters() {

        // alpha
        // beta
        // branch lengths
        // nuc relrate
        // nucstat
        // 2 + Nbranch + 5 + 4

        int N = 2 + Nbranch + 5 + 4;
        double* array = new double[N];
        int i = 0;
        array[i++] = alpha;
        array[i++] = beta;
        for (int j=0; j<Nbranch; j++)   {
            array[i++] = branchlength->GetVal(j);
        }
        for (int j=0; j<Nrr; j++)   {
            array[i++] = nucrelrate[j];
        }
        for (int j=0; j<Nnuc; j++)  {
            array[i++] = nucstat[j];
        }
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveGlobalParameters()   {

        int N = 2 + Nbranch + 5 + 4;
        double* array = new double[N];
        MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        int i = 0;
        alpha = array[i++];
        beta = array[i++];
        for (int j=0; j<Nbranch; j++)   {
            (*branchlength)[j] = array[i++];
        }
        for (int j=0; j<Nrr; j++)   {
            nucrelrate[j] = array[i++];
        }
        for (int j=0; j<Nnuc; j++)  {
            nucstat[j] = array[i++];
        }
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene])    {
                geneprocess[gene]->SetBranchLengths(branchlength);
                geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
                geneprocess[gene]->SetAlphaBeta(alpha,beta);
            }
        }
        delete[] array;
    }

    void SlaveResampleSub()  {

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene])    {
                geneprocess[gene]->UpdateMatrices();
                geneprocess[gene]->ResampleSub();
            }
        }
    }

    void SlaveSendLengthSuffStat()  {

        lengthsuffstatarray->Clear();
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene])    {
                geneprocess[gene] = CollectLengthSuffStat();
                lengthsuffstatarray->Add(geneprocess[gene]->GetLengthSuffStat());
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

    void MasterMoveAlphaBeta()  {

    }

    void SlaveSendNucPathSuffStat()  {

        nucpathsuffstat.Clear();
        for (int gene=0; gene<GetNgene(); gene++)   {
            if (genealloc[gene])    {
                geneprocess[gene] = CollectNucPathSuffStat();
                nucpathsuffstat.Add(geneprocess[gene]->GetNucPathSuffStat());
            }
        }

        int* count = new int[Nnuc+Nnuc*Nnuc];
        double* beta = new double[Nnuc*Nnuc];
        nucpathsuffstat->Push(count,beta);
        MPI_Send(count,Nnuc+Nnuc*Nnuc,MPI_INT,0,TAG1,MPI_COMM_WORLD);
        MPI_Send(beta,Nnuc*Nnuc,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] count;
        delete[] beta;
    }

    void MasterReceiveNucPathSuffStat()  {

        int* count = new int[Nnuc+Nnuc*Nnuc];
        double* beta = new double[Nnuc*Nnuc];
        nucpathsuffstat->Clear();
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nnuc+Nnuc*Nnuc,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Nnuc*Nnuc,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            nucpathsuffstat->Add(count,beta);
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
            if (genealloc[gene])    {
                geneprocess[gene]->MoveOmega();
            }
        }
    }

    void SlaveSendOmega()   {

        double* array = new double[Ngene];
        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene])    {
                (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
            }
            else    {
                (*omegaarray)[gene] = -1;
            }
        }
        MPI_Send(array,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }

    void MasterReceiveOmega()    {

        double* array = new double[Ngene];
        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(array,Ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            for (int gene=0; gene<Ngene; gene++)    {
                if (array[gene] != -1)    {
                    (*omegaarray)[gene] = array[gene];
                }
            }
        }
        delete[] array;
    }

    void MasterSendOmega()  {

        double* array = new double[Ngene];
        omegaarray->Push(array);
        MPI_Bcast(array,Ngene,MPI_DOUBLE,0,MPI_COMM_WORLD);
        delete[] array;
    }

    void SlaveReceiveOmega()    {

        double* array = new double[Ngene];
        MPI_Bcast(array,Ngene,MPI_DOUBLE,0,MPI_COMM_WORLD);

        // this is not useful
        for (int gene=0; gene<Ngene; gene++)    {
            (*omegaarray)[gene] = array[gene];
        }

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene])    {
                geneprocess[gene]->SetOmega(array[gene]);
            }
        }
        delete[] array;
    }

    void SlaveSendLogLikelihood()   {

        for (int gene=0; gene<Ngene; gene++)    {
            if (genealloc[gene])    {
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

