
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

#include "Chrono.hpp"

const int Nrr = Nnuc * (Nnuc-1) / 2;
const int Nstate = 61;

class CodonM8Model	{

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int ncat;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
    //
    // parameters of the distribution of omega across sites
    //

    // mean and inverse concentration of the discretized Beta for 0 < omega < 1
    double purifmean;
    double purifinvconc;

    // positive selection component: omega_pos > 1
    // here dposm = omega_pos - 1 
	double dposom;

    // purifweight[0] : weight of point mass at omega = 0
    // purifweight[1] : weight for 0<omega<1 (beta distribution)
    // purifweight[2] : weight of point mass at omega = 1
    vector<double> purifweight;

    // weight of positive selection component
	double posw;

	DiscBetaWithPos* componentomegaarray;
	MultinomialAllocationVector* sitealloc;
	vector<vector<double> > sitepostprobarray;

    // 
    // hyperparameters of the priors over the mixture parameters
    //

    // prior probability for the gene to be under positive selection (i.e. prior prob that posw > 0)
    double pi;

    // Beta prior for purifmean (with hypermean and hyper inverse concentration)
    double purifmeanhypermean;
	double purifmeanhyperinvconc;

    // Gamma prior for purifinvconc (with hyper mean and hyper inverse shape parameter)
    double purifinvconchypermean;
    double purifinvconchyperinvshape;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    // Beta prior for posw (assuming posw>0)
    double poswhypermean;
    double poswhyperinvconc;

    // Dirichlet prior for the weights for mass at 0, Beta and mass at 1
    // (with hyper center and inverse concentration)
    vector<double> purifweighthypercenter;
    double purifweighthyperinvconc;

    // nucleotide rate parameters
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

    int fixglobal;

	public:

	CodonM8Model(string datafile, string treefile, int inncat, double inpi)	{

        fixglobal = 0;
		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);
		ncat = inncat;
        pi = inpi;

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

        purifweight.assign(3,1.0/3);
        purifweighthypercenter.assign(3,1.0/3);

		// Allocate();
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

        purifmean = 0.5;
        if (ncat == 1)  {
            purifinvconc = 1.0;
        }
        else    {
            purifinvconc = 2.0;
        }
        purifmeanhypermean = 0.5;
        purifmeanhyperinvconc = 0.5;
        purifinvconchypermean = 0.5;
        purifinvconchyperinvshape = 1.0;

		dposom = 1.0;
        dposomhypermean = 0.5;
        dposomhyperinvshape = 0.5;

        if (! pi)   {
            posw = 0;
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }
        else    {
            posw = 0.1;
            poswhypermean = 0.5;
            poswhyperinvconc = 0.5;
        }

        purifweight.assign(3,1.0/3);
        purifweighthypercenter.assign(3,1.0/3);
        purifweighthyperinvconc = 1.0/3;

		componentomegaarray = new DiscBetaWithPos(ncat,purifmean,purifinvconc,posw,dposom+1,purifweight);
		sitealloc = new MultinomialAllocationVector(GetNsite(),componentomegaarray->GetWeights());
		sitepostprobarray.assign(GetNsite(),vector<double>(ncat+3,0));

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
		componentpathsuffstatarray = new PathSuffStatArray(ncat+3);
		siteomegasuffstatarray = new OmegaSuffStatArray(GetNsite());
	}

    double GetPurifMean()   {
        return purifmean;
    }

    double GetPurifInvConc()    {
        return purifinvconc;
    }

    double GetPosW()    {
        return posw;
    }

    double GetDPosOm()   {
        return dposom;
    }

    double GetPurifWeight(int k)    {
        return purifweight[k];
    }

    void SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
        for (int j=0; j<Nbranch; j++)   {
            (*branchlength)[j] = inbranchlength.GetVal(j);
        }
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

    void SetMixtureParameters(double inpurifmean, double inpurifinvconc, double inposw, double indposom, const vector<double>& inpurifweight)    {
        purifmean = inpurifmean;
        purifinvconc = inpurifinvconc;
        posw = inposw;
        dposom = indposom;
        purifweight = inpurifweight;
        componentomegaarray->SetParameters(purifmean,purifinvconc,posw,dposom+1,purifweight);
    }

    void SetMixtureHyperParameters(double inpurifmeanhypermean, double inpurifmeanhyperinvconc, double inpurifinvconchypermean, double inpurifinvconchyperinvshape, double indposomhypermean, double indposomhyperinvshape, double inpi, double inposwhypermean, double inposwhyperinvconc, const vector<double>& inpurifweighthypercenter, double inpurifweighthyperinvconc)  {

        purifmeanhypermean = inpurifmeanhypermean;
        purifmeanhyperinvconc = inpurifmeanhyperinvconc;
        purifinvconchypermean = inpurifinvconchypermean;
        purifinvconchyperinvshape = inpurifinvconchyperinvshape;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
        pi = inpi;
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;
        purifweighthypercenter = inpurifweighthypercenter;
        purifweighthyperinvconc = inpurifweighthyperinvconc;
    }

    // 
    // Matrices
    //

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

    //
    // Likelihood
    //

	double GetLogLikelihood()	{
        // return GetIntegratedLogLikelihood();
		return phyloprocess->GetLogProb();
	}

    double GetIntegratedLogLikelihood() {

        double total = 0;
        double logp[ncat+3];
        const vector<double>& w = componentomegaarray->GetWeights();
        double max = 0;
        for (int i=0; i<GetNsite(); i++) {
            int bkalloc = sitealloc->GetVal(i);

            for (int k=0; k<ncat+3; k++) {
                (*sitealloc)[i] = k;
                logp[k] = phyloprocess->SiteLogLikelihood(i);
                if ((!k) || (max<logp[k]))  {
                    max = logp[k];
                }
            }

            double p = 0;
            for (int k=0; k<ncat+3; k++) {
                p += w[k] * exp(logp[k]-max);
            }
            double logl = log(p) + max;
            total += logl;

            (*sitealloc)[i] = bkalloc;
        }
        return total;
    }

    //
    // Suff Stat and suffstatlogprobs
    //

    const PoissonSuffStatBranchArray* GetLengthSuffStatArray()  {
        return lengthsuffstatarray;
    }

	double LengthSuffStatLogProb()	{
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    const NucPathSuffStat& GetNucPathSuffStat() {
        return nucpathsuffstat;
    }

    double NucPathSuffStatLogProb() {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

	double PathSuffStatLogProb()	{
		return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
	}

	double OmegaSuffStatLogProb()    {
        componentomegaarray->SetParameters(purifmean,purifinvconc,posw,dposom+1,purifweight);
		return componentomegaarray->GetPostProbArray(*siteomegasuffstatarray,sitepostprobarray);
	}

    //
    // Priors
    //

	double GetLogPrior() {
		double total = 0;
		total += LengthLogProb();
		total += LambdaLogProb();
        total += OmegaHyperLogProb();
        total += NucHyperLogProb();
		return total;
	}

    double NucHyperLogProb()    {
        return 0;
    }

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

    // Exponential prior for lambda
	double LambdaLogProb()	{
		return -lambda / 10;
	}

    //
    // Hyper priors for omega mixture
    //
    
    double OmegaHyperLogProb()  {
        double total = 0;
        total += PurifMeanLogProb();
        total += PurifInvConcLogProb();
        total += PosOmegaLogProb();
        total += PosWeightLogProb();
        total += PurifWeightLogProb();
        return total;
    }

    // Beta prior for purifmean
    double PurifMeanLogProb()   {
        double alpha = purifmeanhypermean / purifmeanhyperinvconc;
        double beta = (1-purifmeanhypermean) / purifmeanhyperinvconc;
        return Random::logGamma(alpha+beta) - Random::logGamma(alpha) - Random::logGamma(beta) + (alpha-1)*log(purifmean) + (beta-1)*log(1.0-purifmean);
	}

    // Gamma prior for purifinvconc
    double PurifInvConcLogProb() {
        double alpha = 1.0 / purifinvconchyperinvshape;
        double beta = alpha / purifinvconchypermean;
        return alpha*log(beta) - Random::logGamma(alpha) + (alpha-1)*log(purifinvconc) - beta*purifinvconc;
	}

    // Gamma prior for dposom
	double PosOmegaLogProb()	{
        double alpha = 1.0 / dposomhyperinvshape;
        double beta = alpha / dposomhypermean;
        return alpha*log(beta) - Random::logGamma(alpha) + (alpha-1)*log(dposom) - beta*dposom;
	}

    // mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1 - pi) for posw
	double PosWeightLogProb()   {
        if (posw)   {
            if (! pi)   {
                cerr << "in PosWeightLogProb: pi == 0 and posw > 0\n";
                exit(1);
            }

            double alpha = poswhypermean / poswhyperinvconc;
            double beta = (1 - poswhypermean) / poswhyperinvconc;
            return log(pi) + Random::logGamma(alpha+beta) - Random::logGamma(alpha) - Random::logGamma(beta) + (alpha-1)*log(posw) + (beta-1)*log(1.0-posw);
        }
        else    {
            return log(1-pi);
        }
	}

    // Dirichlet prior for the weights associated with omega = 0, 0<omega<1 and omega=1
    double PurifWeightLogProb() {
        double tot = Random::logGamma(1.0 / purifweighthyperinvconc);
        for (unsigned int k=0; k<purifweight.size(); k++)    {
            tot += -Random::logGamma(purifweighthypercenter[k]/purifweighthyperinvconc) + (purifweighthypercenter[k]/purifweighthyperinvconc-1)*log(purifweight[k]);
        }
        return tot;
    }

    // Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogProb()   {
        if (posw)   {
            return log(pi);
        }
        return log(1-pi);
    }

    //
    //  Moves 
    //

	void Move()	{

        Chrono mappingtime, reptime, lengthtime, collecttime, omegatime, nuctime;

        mappingtime.Start();
		ResampleSub(0.2);
        mappingtime.Stop();

		int nrep = 30;

        reptime.Start();
		for (int rep=0; rep<nrep; rep++)	{

            if (! fixglobal)    {
                lengthtime.Start();
                ResampleBranchLengths();
                MoveLambda();
                lengthtime.Stop();
            }

            collecttime.Start();
			CollectPathSuffStat();
            collecttime.Stop();

            omegatime.Start();
			MoveOmega();
            omegatime.Stop();

            if (! fixglobal)    {
                nuctime.Start();
                UpdateMatrices();
                MoveNuc();
                nuctime.Stop();
            }
		}
        reptime.Stop();

        // cerr << mappingtime.GetTime() << '\t' << reptime.GetTime() << '\t' << lengthtime.GetTime() + collecttime.GetTime() + omegatime.GetTime() + nuctime.GetTime() << '\t' << lengthtime.GetTime() << '\t' << collecttime.GetTime() << '\t' << omegatime.GetTime() << '\t' << nuctime.GetTime() << '\n';
	}

    void ResampleSub(double frac)  {
        UpdateMatrices();
		phyloprocess->Move(frac);
    }

    //
    // Branch Lengths and hyperparam lambda
    //

	void ResampleBranchLengths()	{

        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthsuffstatarray);
	}

    void CollectLengthSuffStat()    {

		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
    }

	void MoveLambda()	{

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

    //
    // Omega mixture 
    //

	void CollectPathSuffStat()	{

		sitepathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
	}

	void CollectComponentPathSuffStat()	{

		componentpathsuffstatarray->Clear();
		sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
	}

	void MoveOmega() 	{

		CollectOmegaSuffStat();

        OmegaHyperSlidingMove(purifmean,0.1,10,0,1);
        if (ncat > 1)   {
            OmegaHyperScalingMove(purifinvconc,0.1,10);
        }
        if (pi != 0)    {
            OmegaHyperScalingMove(dposom,1,10);
            OmegaHyperSlidingMove(posw,1,10,0,1);
        }
        if ((pi != 0) && (pi != 1))    {
            SwitchPosWeight(10);
        }
        OmegaHyperProfileMove(purifweight,0.3,1,10);
        OmegaHyperProfileMove(purifweight,0.1,1,10);

		ResampleAlloc();
	}

	void CollectOmegaSuffStat()	{

		siteomegasuffstatarray->Clear();
		siteomegasuffstatarray->AddSuffStat(*sitecodonmatrixarray,*sitepathsuffstatarray);
	}

	void ResampleAlloc()	{
		OmegaSuffStatLogProb();
		sitealloc->GibbsResample(sitepostprobarray);
	}

	double OmegaHyperScalingMove(double& x, double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - OmegaHyperLogProb() - OmegaSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
		    x *= e;
			deltalogprob += OmegaHyperLogProb() + OmegaSuffStatLogProb();
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

	double OmegaHyperSlidingMove(double& x, double tuning, int nrep, double min = 0, double max = 0)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - OmegaHyperLogProb() - OmegaSuffStatLogProb();
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
			deltalogprob += OmegaHyperLogProb() + OmegaSuffStatLogProb();
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

	double OmegaHyperProfileMove(vector<double>& x, double tuning, int n, int nrep)	{

		double nacc = 0;
		double ntot = 0;
        vector<double> bk(x.size(),0);
		for (int rep=0; rep<nrep; rep++)	{
            bk = x;
			double deltalogprob = - OmegaHyperLogProb() - OmegaSuffStatLogProb();
            double loghastings = Random::ProfileProposeMove(x,x.size(),tuning,n);
			deltalogprob += OmegaHyperLogProb() + OmegaSuffStatLogProb();
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

    double DrawBetaPosWeight()    {
        double alpha = poswhypermean / poswhyperinvconc;
        double beta = (1-poswhypermean) / poswhyperinvconc;
        double a = Random::sGamma(alpha);
        double b = Random::sGamma(beta);
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
			deltalogprob += PosSwitchLogProb() + OmegaSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				posw = bkposw;
			}
			ntot++;
		}
		return nacc/ntot;
	}

    //
    // nucleotide parameters
    //

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


	// summary statistics

	double GetTotalLength()	{
		double tot = 0;
		for (int j=1; j<Nbranch; j++)	{
			tot += branchlength->GetVal(j);
		}
		return tot;
	}

	double GetMeanOmega()   {
		return posw*(1 + dposom) + (1-posw)*(purifweight[1]*purifmean + purifweight[2]);
	}

	double GetEntropy(const std::vector<double>& profile, int dim) const {
		double tot = 0;
		for (int i=0; i<dim; i++)	{
			tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
		}
		return tot;
	}

    void FixGlobalParameters()  {
        fixglobal = 1;
    }

    void GetGlobalParametersFromFile(istream& is)   {
        for (int j=0; j<Nbranch; j++)   {
            is >> (*branchlength)[j];
        }
        for (int j=0; j<Nrr; j++)   {
            is >> nucrelrate[j];
        }
        for (int j=0; j<Nnuc; j++)  {
            is >> nucstat[j];
        }
        UpdateMatrices();
    }

    void GetHyperParametersFromFile(istream& is)    {
        is >> purifmeanhypermean;
        is >> purifmeanhyperinvconc;
        is >> purifinvconchypermean;
        is >> purifinvconchyperinvshape;
        is >> pi;
        is >> poswhypermean;
        is >> poswhyperinvconc;
        is >> dposomhypermean;
        is >> dposomhyperinvshape;
        is >> purifweighthypercenter[0];
        is >> purifweighthypercenter[1];
        is >> purifweighthypercenter[2];
        is >> purifweighthyperinvconc;
    }

	void TraceHeader(std::ostream& os)  {
		os << "#logprior\tlnL\tlength\t";
		os << "meanomega\tposom\tposw\t";
        os << "w0\tw1\t";
		os << "purifmean\tpurifinvconc\t";
		os << "statent\t";
		os << "rrent\t";
        os << "diagerr\n";
	}

	void Trace(ostream& os) {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
		os << GetMeanOmega() << '\t';
		os << dposom+1 << '\t' << posw << '\t';
        os << purifweight[0] << '\t' << purifweight[2] << '\t';
		os << purifmean << '\t' << purifinvconc << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\t';
        os << SubMatrix::diagerr << '\n';
        SubMatrix::diagerr = 0;
	}

	void TracePostProb(ostream& os) {
		for (int i=0; i<GetNsite(); i++)    {
			os << sitepostprobarray[i][ncat+2] << '\t';
		}
		os << '\n';
	}

	void Monitor(ostream& os) {}

	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

};

