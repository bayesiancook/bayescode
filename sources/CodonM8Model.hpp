
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "CodonSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include "cdf.hpp"

const int Nrr = Nnuc * (Nnuc-1) / 2;
const int Nstate = 61;

class MultinomialAllocationVector : public SimpleArray<int> {

	public:

	MultinomialAllocationVector(int insize, const vector<double>& inweight) : SimpleArray<int>(insize), weight(inweight) {
            SampleAlloc();
    }

    ~MultinomialAllocationVector() {}

	void SampleAlloc()  {
        for (int i=0; i<GetSize(); i++) {
            (*this)[i] = Random::DrawFromDiscreteDistribution(weight);
        }
    }
	
	void AddSuffStat(vector<int>& occupancy) const {
        if (occupancy.size() != weight.size())  {
            cerr << "error: non matching size\n";
            exit(1);
        }
        for (unsigned int k=0 ; k<occupancy.size(); k++)  {
            occupancy[k] = 0;
        }
        for (int i=0; i<GetSize(); i++) {
            occupancy[GetVal(i)]++;
        }
    }

	private:
	const vector<double>& weight;
};

template<class T> class FiniteMixture : public virtual ConstArray<T>	{

	public:

	FiniteMixture(const Array<T>* incomponents, const Array<int>* inalloc) : ConstArray<T>(inalloc->GetSize()), components(incomponents), alloc(inalloc)	{
	}

	~FiniteMixture() {}

	const T& GetVal(int i) const override {
		return components->GetVal(alloc->GetVal(i));
	}

	private:
	const Array<T>* components;
	const Array<int>* alloc;
};

class DiscBetaWithPos : public SimpleArray<double>  {

    public:

    DiscBetaWithPos(int inncat, double inalpha, double inbeta, double inposw, double inposom) : SimpleArray<double>(inncat+1) , ncat(inncat), weight(inncat+1), alpha(inalpha), beta(inbeta), posw(inposw), posom(inposom) {
            ComputeDiscBeta();
            (*this)[ncat] = posom;
            ComputeWeights();
    }

    ~DiscBetaWithPos() {}

    void SetAlphaBeta(double inalpha, double inbeta)    {
        alpha = inalpha;
        beta = inbeta;
        ComputeDiscBeta();
    }
    
    void SetPos(double inposw, double inposom)  {
        posw = inposw;
        posom = inposom;
        (*this)[ncat] = posom;
        ComputeWeights();
    }

    const vector<double>& GetWeights() const {return weight;}

    double GetPostProbArray(const OmegaSuffStat& suffstat, vector<double>& postprob) const {

            double logp[ncat+1];
            double max = 0;
            for (int cat=0; cat<ncat+1; cat++)  {
                logp[cat] = suffstat.GetLogProb(GetVal(cat));
                if ((!cat) || (max < logp[cat]))    {
                    max = logp[cat];
                }
            }
            double tot = 0;
            for (int cat=0; cat<ncat+1; cat++)  {
                postprob[cat] = weight[cat] * exp(logp[cat] - max);
                tot += postprob[cat];
            }
            for (int cat=0; cat<ncat+1; cat++)  {
                postprob[cat] /= tot;
            }
            return log(tot) + max;
    }

    double GetPostProbArray(const OmegaSuffStatArray& suffstatarray, vector<vector<double> >& postprobarray)  const {

        double total = 0;
        for (int i=0; i<suffstatarray.GetSize(); i++)   {
            total += GetPostProbArray(suffstatarray.GetConstOmegaSuffStat(i),postprobarray[i]);
        }
        return total;
    }

    private:

    void ComputeDiscBeta()   {
        for (int cat=0; cat<ncat; cat++)  {
            (*this)[cat] = invbetaInc(alpha,beta,((double) (cat+0.5))/ncat);
        }
    }

    void ComputeWeights()   {
        for (int cat=0; cat<ncat; cat++)    {
            weight[cat] = (1-posw)/ncat;
        }
        weight[ncat] = posw;
    }

    int ncat;
    vector<double> weight;
    double alpha;
    double beta;
    double posw;
    double posom;
};

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
	
	double alpha;
	double beta;
    double posw;
    double dposom;
    DiscBetaWithPos* omegaarray;

	MultinomialAllocationVector* alloc;

    vector<vector<double> > postprobarray;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	MGOmegaHeterogeneousCodonSubMatrixArray* codonmatrixarray;
    FiniteMixture<SubMatrix>* submatrixarray;

    PhyloProcess* phyloprocess;

	// suffstats

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;
	OmegaSuffStatArray* omegasuffstatarray;
	PathSuffStatArray* pathsuffstatarray;
	// NucPathSuffStat nucsuffstat;	

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
		cerr << "-- unfold\n";
		phyloprocess->Unfold();
        cerr << "-- logprob\n";
		cerr << phyloprocess->GetLogProb() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
		Trace(cerr);
	}

    int GetNsite() {return data->GetNsite();}

	void Allocate()	{

		lambda = 10.0;
		branchlength = new BranchIIDGamma(tree,1.0,lambda);

		alpha = beta = 1.0;
        posw = 0.1;
        dposom = 0.5;

		omegaarray = new DiscBetaWithPos(ncat,alpha,beta,posw,dposom+1);
        alloc = new MultinomialAllocationVector(GetNsite(),omegaarray->GetWeights());
        postprobarray.assign(GetNsite(),vector<double>(ncat+1,0));

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

		codonmatrixarray = new MGOmegaHeterogeneousCodonSubMatrixArray((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,omegaarray);
        for (int k=0; k<=ncat; k++) {
            cerr << codonmatrixarray->GetVal(k).GetNstate() << '\n';
        }

        submatrixarray = new FiniteMixture<SubMatrix>(codonmatrixarray,alloc);

		phyloprocess = new PhyloProcess(tree,data,branchlength,0,submatrixarray);

		lengthsuffstatarray = new PoissonSuffStatBranchArray(tree);
		pathsuffstatarray = new PathSuffStatArray(GetNsite());
		omegasuffstatarray = new OmegaSuffStatArray(GetNsite());
	}

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

	void UpdateCodonMatrices()	{
		codonmatrixarray->UpdateCodonMatrices();
	}
		
	void UpdateMatrices()   {
		UpdateNucMatrix();
		UpdateCodonMatrices();
	}

	double GetPathSuffStatLogProb()	{
		return pathsuffstatarray->GetLogProb(codonmatrixarray);
	}

    double GetOmegaSuffStatLogProb()    {
        return 0;
    }

    // only for posom
	double PosOmegaLogProb()	{
		return -dposom;
	}

    double PosWeightLogProb()   {
        return 0;
    }

    double AlphaLogProb()   {
        return -alpha;
    }

    double BetaLogProb()    {
        return -beta;
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

		phyloprocess->ResampleSub();

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

			ResampleBranchLengths();
			MoveLambda();

			CollectPathSuffStat();

            MoveAlpha(1,10);
            MoveBeta(1,10);
            MovePosOm(1,10);
            MovePosWeight(1,10);

			MoveRR(0.1,1,3);
			MoveRR(0.03,3,3);
			MoveRR(0.01,3,3);

			MoveNucStat(0.1,1,3);
			MoveNucStat(0.01,1,3);
		}
	}

	void ResampleBranchLengths()	{

		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(lengthsuffstatarray);
		branchlength->GibbsResample(lengthsuffstatarray);
	}

	void MoveLambda()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
		MoveLambda(1.0,10);
		MoveLambda(0.3,10);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		pathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(pathsuffstatarray);
	}

	void CollectOmegaSuffStat()	{

		omegasuffstatarray->Clear();
		omegasuffstatarray->AddSuffStat(*codonmatrixarray,*pathsuffstatarray);
	}

	double MoveRR(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Nrr];
		for (int rep=0; rep<nrep; rep++)	{
			for (int l=0; l<Nrr; l++)	{
				bk[l] = nucrelrate[l];
			}
			double deltalogprob = -GetPathSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucrelrate,Nrr,tuning,n);
			deltalogprob += loghastings;
			UpdateMatrices();
			deltalogprob += GetPathSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nrr; l++)	{
					nucrelrate[l] = bk[l];
				}
                UpdateMatrices();
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
			double deltalogprob = -GetPathSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucstat,Nnuc,tuning,n);
			deltalogprob += loghastings;
            UpdateMatrices();
			deltalogprob += GetPathSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nnuc; l++)	{
					nucstat[l] = bk[l];
				}
                UpdateMatrices();
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MoveAlpha(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - AlphaLogProb() - GetOmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			alpha *= e;
            omegaarray->SetAlphaBeta(alpha,beta);
			deltalogprob += AlphaLogProb() + GetOmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				alpha /= e;
                omegaarray->SetAlphaBeta(alpha,beta);
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MoveBeta(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - BetaLogProb() - GetOmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			beta *= e;
            omegaarray->SetAlphaBeta(alpha,beta);
			deltalogprob += BetaLogProb() + GetOmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				beta /= e;
                omegaarray->SetAlphaBeta(alpha,beta);
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MovePosOm(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - PosOmegaLogProb() - GetOmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			dposom *= e;
            omegaarray->SetPos(posw,dposom+1);
			deltalogprob += PosOmegaLogProb() + GetOmegaSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				dposom /= e;
                omegaarray->SetPos(posw,dposom+1);
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
			double deltalogprob = - PosWeightLogProb() - GetOmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
            posw += m;
            while ((posw<0) && (posw>1))    {
                if (posw <0)    {
                    posw = -posw;
                }
                if (posw > 1)   {
                    posw = 2- posw;
                }
            }
            omegaarray->SetPos(posw,dposom+1);
			deltalogprob += PosWeightLogProb() + GetOmegaSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
                posw = bkposw;
                omegaarray->SetPos(posw,dposom+1);
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
        return 0;
    }

    double GetVarOmega()    {
        return 0;
    }

	double GetEntropy(const std::vector<double>& profile, int dim) const {
		double tot = 0;
		for (int i=0; i<dim; i++)	{
			tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
		}
		return tot;
	}

	void TraceHeader(std::ostream& os)  {
		os << "#logprior\tlnL\tlength\tlambda\t";
		os << "meanomega\tvaromega\talpha\tbeta\t";
        os << "posom\tposw\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
		os << lambda << '\t';
		os << GetMeanOmega() << '\t';
        os << GetVarOmega() << '\t';
        os << alpha << '\t' << beta << '\t';
        os << dposom+1 << '\t' << posw << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
	}

	void Monitor(ostream& os) {}

	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

};

