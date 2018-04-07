
#include "SequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "IIDMultiGamma.hpp"
#include "IIDMultiBernoulli.hpp"
#include "DiffSelSparseFitnessArray.hpp"
#include "IIDProfileMask.hpp"
#include "AASubSelSubMatrixArray.hpp"
#include "PathSuffStat.hpp"
#include "MultiGammaSuffStat.hpp"
#include "GTRSuffStat.hpp"

class SparseAASubSelModel : public ProbModel    {

    // -----
    // external parameters
    // -----

    Tree* tree;
    SequenceAlignment* data;
    const TaxonSet* taxonset;

    // number of sites
    int Nsite;
    int Ntaxa;
    int Nbranch;

    // -----
    //  model structure
    // -----

    // branch lengths 
	double lambda;
	BranchIIDGamma* blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise* branchlength;
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
	GammaSuffStat hyperlengthsuffstat;

    // rates across sites
    double alpha;
    IIDGamma* siterate;
    PoissonSuffStatArray* siteratepathsuffstatarray;
    GammaSuffStat siteratehypersuffstat;

    // aa exchange rates
    vector<double> relratehypercenter;
    double relratehyperinvconc;
	vector<double> relrate;

    double profileshape;
    vector<double> profilecenter;
    IIDMultiGamma* preprofile;

    double maskepsilon;
    double pi;
    IIDProfileMask* sitemaskarray;

    // across conditions and across sites
    MutSelSparseFitnessArray* profile;

    // an array of site-specific matrices
	AASubSelSubMatrixArray* sitematrixarray;

    // phyloprocess
    PhyloProcess* phyloprocess;

    // suff stats
	PathSuffStatArray* sitepathsuffstatarray;
    RelRateSuffStat rrsuffstat;

    int blmode;
    int rrmode;
    int maskepsilonmode;
    int maskmode;
    int ratemode;

    Chrono totchrono;
    Chrono subchrono;
    Chrono paramchrono;
    Chrono profilechrono;
    Chrono maskchrono;

    public:

    SparseAASubSelModel(const string& datafile, const string& treefile, int inratemode, double inepsilon) : rrsuffstat(Naa) {

        blmode = 0;
        rrmode = 0;
        ratemode = inratemode;

        if (inepsilon == 1)   {
            maskepsilon = 1;
            maskmode = 3;
            maskepsilonmode = 3;
        }
        else if (inepsilon >= 0) {
            maskepsilon = inepsilon;
            maskepsilonmode = 3;
            maskmode = 0;
        }
        else    {
            maskepsilonmode = 0;
            maskmode = 0;
            maskepsilon = 0.01;
        }

        ReadFiles(datafile, treefile);
    }

    SparseAASubSelModel(const SparseAASubSelModel& from) = delete;

    ~SparseAASubSelModel() {}

    //! read data and tree files
    void ReadFiles(string datafile, string treefile) {

        // sequence alignment
        data = new FileSequenceAlignment(datafile);
        Nsite = data->GetNsite();  // # columns
        Ntaxa = data->GetNtaxa();
        taxonset = data->GetTaxonSet();
        if (data->GetNstate() != Naa)   {
            cerr << "error: model works only for amino-acids\n";
            exit(1);
        }

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        // traversal of the tree, so as to number links, branches and nodes
        // convention is: branches start at 1 (branch number 0 is the null branch behind the root)
        // nodes start at 0 (for the root), and nodes 1..Ntaxa are tip nodes (corresponding to taxa
        // in sequence alignment)
        tree->SetIndices();
        Nbranch = tree->GetNbranch();
    }

    //! allocate the model (data structures)
    void Allocate() {

        // ----------
        // construction of the model
        // ----------

        // allocating data structures and sampling initial configuration

        // branch lengths
		lambda = 10;
        blhypermean = new BranchIIDGamma(*tree,1.0,lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree,*blhypermean,1.0/blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        alpha = 1.0;
        siterate = 0;
        siteratepathsuffstatarray = 0;
        if (ratemode == 1)  {
            siterate = new IIDGamma(Nsite,alpha,alpha);
            for (int i=0; i<Nsite; i++) {
                (*siterate)[i] = 1.0;
            }
            siteratepathsuffstatarray = new PoissonSuffStatArray(Nsite);
        }

        // relative exchange rates
        relratehypercenter.assign(GetNrr(),1.0/GetNrr());
        relratehyperinvconc = 1.0 / GetNrr();
		relrate.assign(GetNrr(),0);
        SampleRelRate();

        profileshape = 20.0;
        profilecenter.assign(Naa,1.0/Naa);
        preprofile = new IIDMultiGamma(Nsite,Naa,profileshape,profilecenter);

        pi = 0.1;
        sitemaskarray = new IIDProfileMask(Nsite,Naa,pi);

        profile = new MutSelSparseFitnessArray(*preprofile,*sitemaskarray,maskepsilon);
        
        // mut sel matrices (based on the profiles of the mixture)
        sitematrixarray = new AASubSelSubMatrixArray(relrate, profile,false);

		phyloprocess = new PhyloProcess(tree,data,branchlength,siterate,sitematrixarray);
		phyloprocess->Unfold();

        // create suffstat arrays
		sitepathsuffstatarray = new PathSuffStatArray(Nsite);
    }

    void SetBLMode(int in)   {
        blmode = in;
    }

    void SetRelRateMode(int in) {
        rrmode = in;
    }

    void SetMaskMode(int in)    {
        maskmode = in;
    }

    void SetMaskEpsilonMode(int in)    {
        maskepsilonmode = in;
    }

    int GetNrr() const {
        return Naa * (Naa-1) / 2;
    }

    void SampleRelRate()    {
        for (int i=0; i<GetNrr(); i++)  {
            // relrate[i] = Random::Gamma(10,1.0);
            relrate[i] = Random::Gamma(relratehypercenter[i]/relratehyperinvconc,1.0);
        }
    }

    // ------------------
    // Update system
    // ------------------

    //! \brief set branch lengths to a new value
    //! 
    //! Used in a multigene context.
    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    //! get a copy of branch lengths into array given as argument
    void GetBranchLengths(BranchArray<double>& inbranchlength) const    {
        inbranchlength.Copy(*branchlength);
    }

    //! set branch lengths hyperparameters to a new value (multi-gene analyses)
    void SetBranchLengthsHyperParameters(const BranchSelector<double>& inblmean, double inblinvshape)   {
        blhypermean->Copy(inblmean);
        blhyperinvshape = inblinvshape;
        branchlength->SetShape(1.0 / blhyperinvshape);
    }

    //! set relative exchange rate hyperparameters to a new value (multi-gene analyses)
    void SetRelRatesHyperParameters(const vector<double>& inrelratehypercenter, double inrelratehyperinvconc)  {
        relratehypercenter = inrelratehypercenter;
        relratehyperinvconc = inrelratehyperinvconc;
    }

    //! set relative exchange rates to a new value (multi-gene analyses)
    void SetRelRates(const vector<double>& inrelrate)  {
        relrate = inrelrate;
        CorruptMatrices();
    }

    //! copy relative exchange rates into vectors given as arguments (multi-gene analyses)
    void GetRelRates(vector<double>& inrelrate) const  {
        inrelrate = relrate;
    }

    //! \brief set value of background fitness of low-fitness amino-acids
    void SetMaskEpsilon(double in)  {
        maskepsilon = in;
        profile->SetEpsilon(maskepsilon);
    }

    //! const ref access to masks across sites
    const vector<vector<int> >& GetMaskArray() const    {
        return sitemaskarray->GetArray();
    }

    //! const access to low-fitness background value (mask epsilon)
    double GetMaskEpsilon() const   {
        return maskepsilon;
    }

    void Update() override {
        if (blmode == 0)    {
            blhypermean->SetAllBranches(1.0/lambda);
        }
        if (ratemode == 1)  {
            siterate->SetShape(alpha);
            siterate->SetScale(alpha);
        }
        UpdateMask();
		preprofile->SetShape(profileshape);
        UpdateAll();
        ResampleSub(1.0);
    }

    //! \brief dummy function that does not do anything.
    //! 
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove (defined in ProbModel),
    //! all of which require a void (*f)(void) function pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief tell the gtr matrices that their parameters have changed and that it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is notified,
    //! such that the matrices know that they will have to recalculate whichever component is requested later on upon demand.
    void CorruptMatrices()  {
        sitematrixarray->UpdateMatrices();
    }

    //! update profiles and matrices across all sites and conditions
    void UpdateAll() {
        profile->SetEpsilon(maskepsilon);
        profile->Update();
        CorruptMatrices();
    }

    //! update profiles and matrices across all conditions for site i
    void UpdateSite(int i) {
        profile->Update(i);
        sitematrixarray->UpdateMatrix(i);
    }

    //! update parameters of sitemaskarray (probability parameter pi)
    void UpdateMask()   {
        sitemaskarray->SetPi(pi);
    }

    // ---------------
    // log priors
    // ---------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) {
            total += BranchLengthsLogPrior();
        }
        if (ratemode == 1)  {
            total += SiteRateHyperLogPrior();
            total += SiteRateLogPrior();
        }
        if (rrmode < 2)    {
            total += RelRatesLogPrior();
        }
        if (maskmode < 2)   {
            total += MaskHyperLogPrior();
            total += MaskLogPrior();
        }
        return total;
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here, lambda ~ exponential of rate 10)
	double BranchLengthsHyperLogPrior() const {
		return -log(10.0) - lambda / 10;
	}

    //! log prior over branch lengths (iid exponential of rate lambda)
	double BranchLengthsLogPrior() const {
		double ret = branchlength->GetLogProb();
        if (blmode == 0)    {
            ret += BranchLengthsHyperLogPrior();
        }
        return ret;
	}

    double SiteRateHyperLogPrior() const {
        return -alpha;
    }

    double SiteRateLogPrior() const {
        return siterate->GetLogProb();
    }

    //! log prior over nuc rates rho and pi (uniform)
    double RelRatesLogPrior() const {
        double ret = 0;
        for (int i=0; i<GetNrr(); i++)  {
            ret += Random::logGammaDensity(relrate[i],relratehypercenter[i]/relratehyperinvconc,1.0);
        }
        return ret;
    }

    //! log prior over mask array hyperparameters
    double MaskHyperLogPrior() const  {
        double ret = 0;
        if (maskepsilonmode < 2)    {
            ret -= 10*maskepsilon;
        }
        return ret;
    }

    //! log prior over mask array
    double MaskLogPrior() const   {
        return sitemaskarray->GetLogProb();
    }

    //! log prior over mask array (only for site i)
    double MaskLogPrior(int i) const   {
        return sitemaskarray->GetLogProb(i);
    }

    //! return log likelihood
    double GetLogLikelihood() const { 
        return phyloprocess->GetLogLikelihood();
    }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const {
        return GetLogPrior() + GetLogLikelihood();
    }

    // ---------------
    // collecting suff stats
    // ---------------

    //! \brief const access to array of length-pathsuffstats across branches
    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! collect sufficient statistics if substitution mappings across sites
	void CollectSitePathSuffStat()	{
		sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
	}

    //! collect sufficient statistics for moving branch lengths (directly from the substitution mappings)
    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! collect sufficient statistics for moving site rates (directly from the substitution mappings)
    void CollectSiteRateSuffStat()    {
		siteratepathsuffstatarray->Clear();
        siteratepathsuffstatarray->AddRatePathSuffStat(*phyloprocess);
    }

    //! return log prob of the current substitution mapping, as a function of the current substitution process
	double SuffStatLogProb() const {
        return sitepathsuffstatarray->GetLogProb(*sitematrixarray);
	}

    //! return log prob of the substitution mappings for site i
    double SiteSuffStatLogProb(int i) const {
        return sitepathsuffstatarray->GetVal(i).GetLogProb(sitematrixarray->GetVal(i));
    }

    //! \brief return log prob of current branch lengths, as a function of branch lengths hyperparameter lambda
	double BranchLengthsHyperSuffStatLogProb() const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    //! \brief return log prob of current site rates, as a function of hyperparameter alpha
    double SiteRateHyperSuffStatLogProb() const {
        return siteratehypersuffstat.GetLogProb(alpha,alpha);
    }

    // ---------------
    // log probs for MH moves
    // ---------------

    //! \brief log prob factor to be recomputed when moving branch lengths hyperparameters (here, lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving mask hyperparameter pi
    double MaskLogProb() const  {
        return MaskHyperLogPrior() + MaskLogPrior();
    }

    //! \brief log prob factor to be recomputed when moving maskepsilon
    double MaskEpsilonLogProb() const  {
        return MaskHyperLogPrior() + SuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving alpha
    double SiteRateHyperLogProb() const  {
        return SiteRateHyperLogPrior() + SiteRateHyperSuffStatLogProb();
    }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
	double Move() override {
        totchrono.Reset();
        subchrono.Reset();
        paramchrono.Reset();
        profilechrono.Reset();
        maskchrono.Reset();

        totchrono.Start();
        subchrono.Start();
        ResampleSub(1.0);
        subchrono.Stop();
        paramchrono.Start();
        MoveParameters(3,1);
        paramchrono.Stop();
        totchrono.Stop();
        return 1.0;
	}

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep0, int nrep) {

        for (int rep0 = 0; rep0 < nrep0; rep0++) {
            if (blmode < 2)    {
                MoveBranchLengths();
            }
            if (ratemode == 1)  {
                MoveSiteRates();
            }
            CollectSitePathSuffStat();
            UpdateAll();
            for (int rep = 0; rep < nrep; rep++) {

                profilechrono.Start();
                MovePreProfile(10);
                CompMovePreProfile(3);
                profilechrono.Stop();

                maskchrono.Start();
                if (maskmode < 3)   {
                    MoveMasks(20);
                }
                if (maskmode < 2)   {
                    MoveMaskHyperParameters();
                }
                if (maskepsilonmode < 2)    {
                    MoveMaskEpsilon(10);
                }
                maskchrono.Stop();
            }

            if (rrmode < 2)    {
                MoveRelRates();
                RelRateLengthCompMoveSchedule(10);
                CorruptMatrices();
            }
        }

        UpdateAll();
    }

    //! Gibbs resampling of substitution histories conditional on current parameter configuration
    void ResampleSub(double frac)   {
        CorruptMatrices();
		phyloprocess->Move(frac);
    }

    //! Gibbs resampling of branch lengths (based on sufficient statistics and current value of lambda)
	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}

    //! MCMC move schedule on branch lengths 
    void MoveBranchLengths()    {
        ResampleBranchLengths();
        if (blmode == 0)    {
            MoveLambda();
        }
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda, based on suffstats for branch lengths)
	void MoveLambda()	{
		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda,1.0,10,&SparseAASubSelModel::BranchLengthsHyperLogProb,&SparseAASubSelModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&SparseAASubSelModel::BranchLengthsHyperLogProb,&SparseAASubSelModel::NoUpdate,this);
        blhypermean->SetAllBranches(1.0/lambda);
	}

    void MoveSiteRates()    {
        ResampleSiteRates();
        MoveAlpha();
    }

    void ResampleSiteRates()    {
        CollectSiteRateSuffStat();
        siterate->GibbsResample(*siteratepathsuffstatarray);
    }

    void MoveAlpha()    {
        siteratehypersuffstat.Clear();
        siteratehypersuffstat.AddSuffStat(*siterate);
        ScalingMove(alpha,1.0,10,&SparseAASubSelModel::SiteRateHyperLogProb,&SparseAASubSelModel::NoUpdate,this);
        ScalingMove(alpha,0.3,10,&SparseAASubSelModel::SiteRateHyperLogProb,&SparseAASubSelModel::NoUpdate,this);
        siterate->SetShape(alpha);
        siterate->SetScale(alpha);
	}

    void CollectRelRateSuffStat()   {
        rrsuffstat.Clear();
        rrsuffstat.AddSuffStat(*sitematrixarray,*sitepathsuffstatarray);
    }

    //! MH moves on relative exchange rates
	void MoveRelRates()	{
        CollectRelRateSuffStat();
        for (int i=0; i<GetNrr(); i++)  {
            relrate[i] = Random::Gamma(relratehypercenter[i]/relratehyperinvconc + rrsuffstat.GetCount(i), 1.0 + rrsuffstat.GetBeta(i));
        }
        CorruptMatrices();
	}

    double RelRateLengthCompMoveSchedule(int nrep)  {
        for (int rep=0; rep<nrep; rep++)    {
            RelRateLengthCompMove(1);
            RelRateLengthCompMove(0.1);
        }
        return 1.0;
    }

    double RelRateLengthCompMove(double tuning)   {

        double m = tuning*(Random::Uniform() - 0.5);
        double e = exp(m);
        double logratio = -BranchLengthsLogPrior() - RelRatesLogPrior();
        for (int j=0; j<Nbranch; j++)   {
            (*branchlength)[j] *= e;
        }
        for (int i=0; i<GetNrr(); i++)  {
            relrate[i] /= e;
        }
        double loghastings = m*(Nbranch-GetNrr());
        logratio += loghastings;
        logratio += BranchLengthsLogPrior() + RelRatesLogPrior();
        int accept = (log(Random::Uniform()) < logratio);
        if (! accept)   {
            for (int j=0; j<Nbranch; j++)   {
                (*branchlength)[j] /= e;
            }
            for (int i=0; i<GetNrr(); i++)  {
                relrate[i] *= e;
            }
        }
        return accept;
    }

    //! MH compensatory move schedule on parameters and hyper-parameters
    void CompMovePreProfile(int nrep)  {
        CompMovePreProfile(1.0,nrep);
    }

    //! \brief MH compensatory move on parameters and hyper-parameters
    //!
    //! for a given amino-acid.
    //! shift by a constant factor e all *active* gamma parameters across all amino-acids and conditions
    double CompMovePreProfile(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*preprofile)[i];
                const vector<int>& mask = (*sitemaskarray)[i];

                double deltalogprob = - SiteSuffStatLogProb(i);

                for (int a=0; a<Naa; a++)   {
                    if (mask[a])	{
                        double alpha = profileshape*profilecenter[a];
                        deltalogprob -= - Random::logGamma(alpha) + (alpha-1)*log(x[a]) - x[a];
                    }
                }

                double m = tuning*(Random::Uniform() - 0.5);
                double e = exp(m);

                int n = 0;
                for (int a=0; a<Naa; a++)   {
                    if (mask[a])	{
                        x[a] *= e;
                        n++;
                    }
                }

                double loghastings = n * m;

                for (int a=0; a<Naa; a++)   {
                    if (mask[a])	{
                        double alpha = profileshape*profilecenter[a];
                        deltalogprob += - Random::logGamma(alpha) + (alpha-1)*log(x[a]) - x[a];
                    }
                }

                UpdateSite(i);
                deltalogprob += SiteSuffStatLogProb(i);

                deltalogprob += loghastings;

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    for (int a=0; a<Naa; a++)   {
                        if (mask[a])  {
                            x[a] /= e;
                        }
                    }
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on baseline gamma parameters (for condition k=0)
    void MovePreProfile(int nrep) {
        // if masks are not activated (all entries equal to 1), move a random subset of entries over the 20 amino-acids (2d parameter of call)
        if (maskmode == 3)  {
            MovePreProfileAll(1.0, 1, nrep);
            MovePreProfileAll(1.0, 3, nrep);
            MovePreProfileAll(1.0, 20, nrep);
            MovePreProfileAll(0.3, 20, nrep);
        }
        // if masks are activated, move all active entries
        else    {
            MovePreProfile(1.0, nrep);
            MovePreProfile(0.3, nrep);
        }
    }

    //! elementary MH move on baseline gamma parameters (for condition k=0)
    double MovePreProfile(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*preprofile)[i];
                const vector<int>& s = (*sitemaskarray)[i];

                bk = x;

                double deltalogprob = -preprofile->GetLogProb(i,s) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, s);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += preprofile->GetLogProb(i,s) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! elementary MH move on baseline fitness parameters (for condition k=0): version used when masks are activated
    double MovePreProfileAll(double tuning, int n, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*preprofile)[i];

                bk = x;

                double deltalogprob = -preprofile->GetLogProb(i) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, n);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += preprofile->GetLogProb(i) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on hyperparameter of mask across sites (maskprob, or pi)
    void MoveMaskHyperParameters()  {
        SlidingMove(pi,1.0,10,0.05,0.975,&SparseAASubSelModel::MaskLogProb,&SparseAASubSelModel::UpdateMask,this);
        SlidingMove(pi,0.1,10,0.05,0.975,&SparseAASubSelModel::MaskLogProb,&SparseAASubSelModel::UpdateMask,this);
    }

    //! MH move schedule on background fitness (maskepsilon)
    void MoveMaskEpsilon(int nrep)  {
        SlidingMove(maskepsilon,1.0,nrep,0,1.0,&SparseAASubSelModel::MaskEpsilonLogProb,&SparseAASubSelModel::UpdateAll,this);
        SlidingMove(maskepsilon,0.1,nrep,0,1.0,&SparseAASubSelModel::MaskEpsilonLogProb,&SparseAASubSelModel::UpdateAll,this);
    }

    //! MH move on masks across sites
    double MoveMasks(int nrep)    {
		double nacc = 0;
		double ntot = 0;
        for (int i=0; i<Nsite; i++) {
            vector<int>& mask = (*sitemaskarray)[i];
            int naa = 0;
            for (int k=0; k<Naa; k++)   {
                naa += mask[k];
            }
            for (int rep=0; rep<nrep; rep++)    {
                int k = (int) (Naa * Random::Uniform());
                if (nrep == 20) {
                    k = rep;
                }
                // for (int k=0; k<Naa; k++)   {
                if ((!mask[k]) || (naa > 1))    {
                    double deltalogprob = -MaskLogPrior(i) - SiteSuffStatLogProb(i);
                    naa -= mask[k];
                    mask[k] = 1-mask[k];
                    naa += mask[k];
                    if (mask[k])    {
                        (*preprofile)[i][k] = Random::sGamma(profileshape * profilecenter[k]);
                        if (! (*preprofile)[i][k]) {
                            (*preprofile)[i][k] = 1e-8;
                        }
                    }
                    UpdateSite(i);
                    deltalogprob += MaskLogPrior(i) + SiteSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted)	{
                        nacc ++;
                    }
                    else	{
                        naa -= mask[k];
                        mask[k] = 1-mask[k];
                        naa += mask[k];
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
		return nacc/ntot;
	}

    //-------------------
    // Accessors
    // ------------------

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //-------------------
    // Traces and monitors
    // ------------------

    //! write header of tracefile
    void TraceHeader(ostream& os) const override {
        os << "#logprior\tlnL\tlength\t";
        if (ratemode == 1)  {
            os << "alpha\tmeanrate\tvarrate\t";
        }
        os << "pi\t";
        os << "width\t";
        os << "epsilon\t";
        os << "statent\t";
        os << "meanrr\t";
        os << "rrent\n";
    }

    //! write trace (one line summarizing current state) into trace file
    void Trace(ostream& os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        if (ratemode == 1)  {
            os << alpha << '\t';
            os << GetMeanSiteRate() << '\t' << GetVarSiteRate() << '\t';
        }
        os << pi << '\t';
        os << sitemaskarray->GetMeanWidth() << '\t';
        os << maskepsilon << '\t';
        os << profile->GetMeanEntropy() << '\t';
        os << GetMeanRelRate() << '\t' << GetRelRateEntropy() << '\n';
    }

    double GetMeanSiteRate() const  {
        double tot = 0;
        for (int i=0; i<Nsite; i++) {
            tot += (*siterate)[i];
        }
        tot /= Nsite;
        return tot;
    }

    double GetVarSiteRate() const   {
        double mean = 0;
        double var = 0;
        for (int i=0; i<Nsite; i++) {
            mean += (*siterate)[i];
            var += (*siterate)[i] * (*siterate)[i];
        }
        mean /= Nsite;
        var /= Nsite;
        var -= mean*mean;
        return var;
    }

    double GetRelRateEntropy() const    {

        double tot = 0;
        for (unsigned int i=0; i<relrate.size(); i++)    {
            tot += relrate[i];
        }
        double ret = 0;
        for (unsigned int i=0; i<relrate.size(); i++)    {
            double tmp = relrate[i] / tot;
            ret -= tmp*log(tmp);
        }
        return ret;
    }

    double GetMeanRelRate() const   {
        double tot = 0;
        for (unsigned int i=0; i<relrate.size(); i++)    {
            tot += relrate[i];
        }
        tot /= relrate.size();
        return tot;
    }

    //! monitoring MCMC statistics
    void Monitor(ostream& os) const override {
        os << "tot   : " << totchrono.GetTime() << '\n';
        os << "sub   : " << 100 * subchrono.GetTime() / totchrono.GetTime() << '\n';
        os << "param : " << 100 * paramchrono.GetTime() / totchrono.GetTime() << '\n';
        os << "prof  : " << 100 * profilechrono.GetTime() / totchrono.GetTime() << '\n';
        os << "mask  : " << 100 * maskchrono.GetTime() / totchrono.GetTime() << '\n';
    }

    //! get complete parameter configuration from stream
    void FromStream(istream& is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (rrmode < 2)    {
            is >> relrate;
        }
        is >> *preprofile;
        if (maskmode < 2)   {
            is >> pi;
        }
        if (maskmode < 3)   {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            is >> maskepsilon;
        }
    }

    //! write complete current parameter configuration to stream
    void ToStream(ostream& os) const override {
        if (blmode < 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        }
        if (rrmode < 2)    {
            os << relrate << '\t';
        }
        os << *preprofile << '\t';
        if (maskmode < 2)   {
            os << pi << '\t';
        }
        if (maskmode < 3)   {
            os << *sitemaskarray << '\t';
        }
        if (maskepsilonmode < 2)    {
            os << maskepsilon << '\t';
        }
    }

    //! return size of model, when put into an MPI buffer (in multigene context)
    unsigned int GetMPISize() const {
        int size = 0;
        if (blmode < 2) {
            size++;
            size += branchlength->GetMPISize();
        }
        if (rrmode < 2)    {
            size += relrate.size();
        }
        size += preprofile->GetMPISize();
        if (maskmode < 2)   {
            size++;
        }
        if (maskmode < 3)   {
            size += sitemaskarray->GetMPISize();
        }
        if (maskepsilonmode < 2)    {
            size++;
        }
        return size;
    }

    //! get complete parameter configuration from MPI buffer
    void MPIGet(const MPIBuffer& is)    {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (rrmode < 2)    {
            is >> relrate;
        }
        is >> *preprofile;
        if (maskmode < 2)   {
            is >> pi;
        }
        if (maskmode < 3)   {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            is >> maskepsilon;
        }
    }

    //! write complete current parameter configuration into MPI buffer
    void MPIPut(MPIBuffer& os) const {
        if (blmode < 2) {
            os << lambda;
            os << *branchlength;
        }
        if (rrmode < 2)    {
            os << relrate;
        }
        os << *preprofile;
        if (maskmode < 2)   {
            os << pi;
        }
        if (maskmode < 3)   {
            os << *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            os << maskepsilon;
        }
    }
};

