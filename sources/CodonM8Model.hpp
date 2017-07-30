
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "CodonSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include "DiscBetaWithPos.hpp"
#include "MultinomialAllocationVector.hpp"

const int Nrr = Nnuc * (Nnuc-1) / 2;
const int Nstate = 61;

class CodonM8Model	{

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int ncat;
	int withpos;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	double aalpha;
    double abeta;
    double balpha;
    double bbeta;

    double pi;
    double poswalpha;
    double poswbeta;

    double dposomalpha;
    double dposombeta;

	double alpha;
	double beta;
	double posw;

	double dposom;
	DiscBetaWithPos* componentomegaarray;

	MultinomialAllocationVector* sitealloc;

	vector<vector<double> > sitepostprobarray;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	// array of matrices across components of the mixture
	MGOmegaCodonSubMatrixArray* componentcodonmatrixarray;

	// arrays of matrices across sites (such as determined by the site allocations to the mixture components)
	// two versions differing only by their exact type

	// used for collecting omega suffstats: need to have access to the *codon* matrix for each site
	ConstMixtureArray<MGOmegaCodonSubMatrix>* sitecodonmatrixarray;

	// used by PhyloProcess: has to be a ConstArray<SubMatrix>
	ConstMixtureArray<SubMatrix>* sitesubmatrixarray;

	PhyloProcess* phyloprocess;

	// suffstats

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;
	OmegaSuffStatArray* siteomegasuffstatarray;
	PathSuffStatArray* sitepathsuffstatarray;
	PathSuffStatArray* componentpathsuffstatarray;

    NucPathSuffStat nucpathsuffstat;

	public:

	CodonM8Model(string datafile, string treefile, int inncat, int inwithpos)	{

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);
		ncat = inncat;
		withpos = inwithpos;

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

		std::cerr << "-- Number of sites: " << Nsite << std::endl;

		taxonset = codondata->GetTaxonSet();

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

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogProb() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
		// Trace(cerr);
    }

	int GetNsite() {return codondata->GetNsite();}

    CodonStateSpace* GetCodonStateSpace()   {
		return (CodonStateSpace*) codondata->GetStateSpace();
    }

	void Allocate()	{

		lambda = 10.0;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);

        aalpha = abeta = 1.0;
        balpha = bbeta = 1.0;
        pi = 0.1;
        poswalpha = poswbeta = 1.0;
        dposomalpha = dposombeta = 1.0;

		alpha = beta = 1.0;
		posw = 0.1;
        if (! withpos)  {
            posw = 0;
        }
		dposom = 0.5;

		componentomegaarray = new DiscBetaWithPos(ncat,alpha,beta,posw,dposom+1);
		sitealloc = new MultinomialAllocationVector(GetNsite(),componentomegaarray->GetWeights());
		sitepostprobarray.assign(GetNsite(),vector<double>(ncat+1,0));

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

		componentcodonmatrixarray = new MGOmegaCodonSubMatrixArray((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,componentomegaarray);

		sitesubmatrixarray = new ConstMixtureArray<SubMatrix>(componentcodonmatrixarray,sitealloc);
		sitecodonmatrixarray = new ConstMixtureArray<MGOmegaCodonSubMatrix>(componentcodonmatrixarray,sitealloc);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitesubmatrixarray);

		lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);
		sitepathsuffstatarray = new PathSuffStatArray(GetNsite());
		componentpathsuffstatarray = new PathSuffStatArray(ncat+1);
		siteomegasuffstatarray = new OmegaSuffStatArray(GetNsite());
	}

    void SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
        for (int j=0; j<Nbranch; j++)   {
            (*branchlength)[j] = inbranchlength.GetVal(j);
        }
    }

    const PoissonSuffStatBranchArray* GetLengthSuffStatArray()  {
        return lengthsuffstatarray;
    }

    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        for (int j=0; j<Nrr; j++)   {
            nucrelrate[j] = innucrelrate[j];
        }
        for (int j=0; j<Nnuc; j++)  {
            nucstat[j] = innucstat[j];
        }
        UpdateMatrices();
    }

    void SetMixtureParameters(double inalpha, double inbeta, double inposw, double indposom)    {
        alpha = inalpha;
        beta = inbeta;
        posw = inposw;
        dposom = indposom;
        componentomegaarray->SetAlphaBeta(alpha,beta);
    }

    void SetMixtureHyperParameters(double inaalpha, double inabeta, double inbalpha, double inbbeta, double inpi, double inposwalpha, double inposwbeta, double indposomalpha, double indposombeta)  {
        aalpha = inaalpha;
        abeta = inabeta;
        balpha = inbalpha;
        bbeta = inbbeta;
        pi = inpi;
        poswalpha = inposwalpha;
        poswbeta = inposwbeta;
        dposomalpha = indposomalpha;
        dposombeta = indposombeta;
    }

    double GetAlpha() {
        return alpha;
    }

    double GetBeta()  {
        return beta;
    }

    double GetPosW()  {
        return posw;
    }

    double GetDPosOm()    {
        return dposom;
    }

    const NucPathSuffStat& GetNucPathSuffStat() {
        return nucpathsuffstat;
    }

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

	void UpdateCodonMatrices()	{
		componentcodonmatrixarray->UpdateCodonMatrices();
	}
		
	void UpdateMatrices()   {
		UpdateNucMatrix();
		UpdateCodonMatrices();
	}

	double PathSuffStatLogProb()	{
		return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
	}

	double OmegaSuffStatLogProb()    {
		return componentomegaarray->GetPostProbArray(*siteomegasuffstatarray,sitepostprobarray);
	}

	double PosOmegaLogProb()	{
        return dposomalpha*log(dposombeta) - Random::logGamma(dposomalpha) + (dposomalpha-1)*log(dposom) - dposombeta*dposom;
	}

	double PosWeightLogProb()   {
        if (posw)   {
            if (! pi)   {
                cerr << "in PosWeightLogProb: pi == 0 and posw > 0\n";
                exit(1);
            }
            return log(pi) + Random::logGamma(poswalpha+poswbeta) - Random::logGamma(poswalpha) - Random::logGamma(poswbeta) + (poswalpha-1)*log(posw) + (poswbeta-1)*log(1.0-posw);
        }
        else    {
            return log(1-pi);
        }
	}

    double PosSwitchLogProb()   {
        if (posw)   {
            return log(pi);
        }
        return log(1-pi);
    }

	double AlphaLogProb()   {
		return aalpha*log(abeta) - Random::logGamma(aalpha) + (aalpha-1)*log(alpha) - abeta*alpha;
	}

	double BetaLogProb()    {
        return balpha*log(bbeta) - Random::logGamma(balpha) + (balpha-1)*log(beta) - bbeta*beta;
	}

	double LambdaLogProb()	{
		return -lambda / 10;
	}

	double LengthSuffStatLogProb()	{
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

	void Move()	{

		ResampleSub();

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

			ResampleBranchLengths();
			MoveLambda();

			CollectPathSuffStat();

			MoveOmega();
			UpdateMatrices();
			MoveNuc();
		}
	}

    void ResampleSub()  {
        UpdateMatrices();
		phyloprocess->ResampleSub();
    }

    void CollectLengthSuffStat()    {

		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
    }

	void ResampleBranchLengths()	{

        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthsuffstatarray);
	}

	void MoveLambda()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
		MoveLambda(1.0,10);
		MoveLambda(0.3,10);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		sitepathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
	}

	void CollectComponentPathSuffStat()	{

		componentpathsuffstatarray->Clear();
		sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
	}

	void CollectOmegaSuffStat()	{

		siteomegasuffstatarray->Clear();
		siteomegasuffstatarray->AddSuffStat(*sitecodonmatrixarray,*sitepathsuffstatarray);
	}

	void MoveOmega() 	{

		CollectOmegaSuffStat();
		MoveAlpha(0.1,10);
		MoveBeta(0.1,10);
		MovePosOm(1,10);
		MovePosWeight(1,10);
        if (withpos == 1)    {
            SwitchPosWeight(10);
        }
		ResampleAlloc();
	}

    double NucPathSuffStatLogProb() {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

    void CollectNucPathSuffStat()   {
		UpdateMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
    }

	void MoveNuc()	{

		CollectComponentPathSuffStat();
        CollectNucPathSuffStat();

		MoveRR(0.1,1,3);
		MoveRR(0.03,3,3);
		MoveRR(0.01,3,3);

		MoveNucStat(0.1,1,3);
		MoveNucStat(0.01,1,3);

        UpdateMatrices();
	}

	void ResampleAlloc()	{
		OmegaSuffStatLogProb();
		sitealloc->GibbsResample(sitepostprobarray);
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

	double MoveAlpha(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - AlphaLogProb() - OmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			alpha *= e;
			componentomegaarray->SetAlphaBeta(alpha,beta);
			deltalogprob += AlphaLogProb() + OmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				alpha /= e;
				componentomegaarray->SetAlphaBeta(alpha,beta);
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
			componentomegaarray->SetAlphaBeta(alpha,beta);
			deltalogprob += BetaLogProb() + OmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				beta /= e;
				componentomegaarray->SetAlphaBeta(alpha,beta);
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MovePosOm(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - PosOmegaLogProb() - OmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			dposom *= e;
			componentomegaarray->SetPos(posw,dposom+1);
			deltalogprob += PosOmegaLogProb() + OmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				dposom /= e;
				componentomegaarray->SetPos(posw,dposom+1);
			}
			ntot++;
		}
		return nacc/ntot;
	}

    double DrawBetaPosWeight()    {
        double a = Random::sGamma(poswalpha);
        double b = Random::sGamma(poswbeta);
        double ret = a / (a+b);
        return ret;
    }

	double SwitchPosWeight(int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double bkposw = posw;
			double deltalogprob = - PosSwitchLogProb() - OmegaSuffStatLogProb();
            if (posw)   {
                posw = 0;
            }
            else    {
                posw = DrawBetaPosWeight();
            }
			componentomegaarray->SetPos(posw,dposom+1);
			deltalogprob += PosSwitchLogProb() + OmegaSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				posw = bkposw;
				componentomegaarray->SetPos(posw,dposom+1);
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MovePosWeight(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double bkposw = posw;
			double deltalogprob = - PosWeightLogProb() - OmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			posw += m;
			while ((posw<0) || (posw>1))    {
				if (posw <0)    {
					posw = -posw;
				}
				if (posw > 1)   {
					posw = 2- posw;
				}
			}
			componentomegaarray->SetPos(posw,dposom+1);
			deltalogprob += PosWeightLogProb() + OmegaSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				posw = bkposw;
				componentomegaarray->SetPos(posw,dposom+1);
			}
			ntot++;
		}
		return nacc/ntot;
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

	// summary statistics

	double GetLogPrior() {
		double total = 0;
		total += LambdaLogProb();
		total += LengthLogProb();
		total += AlphaLogProb();
		total += BetaLogProb();
		total += PosOmegaLogProb();
		total += PosWeightLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}

	double GetTotalLength()	{
		double tot = 0;
		for (int j=1; j<Nbranch; j++)	{
			tot += branchlength->GetVal(j);
		}
		return tot;
	}

	double GetMeanOmega()   {
		return (1-posw)*alpha/(alpha+beta) + posw*(1+dposom);
	}

	double GetEntropy(const std::vector<double>& profile, int dim) const {
		double tot = 0;
		for (int i=0; i<dim; i++)	{
			tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
		}
		return tot;
	}

	void TraceHeader(std::ostream& os)  {
		os << "#logprior\tlnL\tlength\t";
		os << "meanomega\tposom\tposw\t";
		os << "alpha\tbeta\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
		os << GetMeanOmega() << '\t';
		os << dposom+1 << '\t' << posw << '\t';
		os << alpha << '\t' << beta << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
	}

	void TracePostProb(ostream& os) {
		for (int i=0; i<GetNsite(); i++)    {
			os << sitepostprobarray[i][ncat] << '\t';
		}
		os << '\n';
	}

	void Monitor(ostream& os) {}

	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

};

