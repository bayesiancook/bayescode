
#include "SequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "T92SubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "IIDMultiGamma.hpp"
#include "IIDMultiBernoulli.hpp"
#include "DiffSelSparseFitnessArray.hpp"
#include "IIDProfileMask.hpp"
#include "AAMutSelSubMatrixArray.hpp"
#include "PathSuffStat.hpp"
#include "MultiGammaSuffStat.hpp"

class SparseAAMutSelModel : public ProbModel    {

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

    // mutation process
    double kappa;
    double gc;
    double xi;
    T92SubMatrix* nucmatrix;
    CodonStateSpace* codonstatespace;

    double profileshape;
    vector<double> profilecenter;
    IIDMultiGamma* preprofile;

    double maskepsilon;
    // double pi;
    // IIDProfileMask* sitemaskarray;
    vector<double> pi;
    ProfileMask* sitemaskarray;

    // across conditions and across sites
    MutSelSparseFitnessArray* profile;

    // an array of site-specific matrices
	AAMutSelSubMatrixArray* sitematrixarray;

    // phyloprocess
    PhyloProcess* phyloprocess;

    // suff stats
	PathSuffStatArray* sitepathsuffstatarray;

    int blmode;
    int maskepsilonmode;
    int maskmode;
    int ratemode;
    int profilemode;
    int ximode;
    int gcmode;

    public:

    SparseAAMutSelModel(const string& datafile, const string& treefile, int inratemode, int inprofilemode, double inepsilon, double ingc, double inxi)   {

        blmode = 0;
        ratemode = inratemode;
        profilemode = inprofilemode;

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

        xi = inxi;
        if (xi == -1)   {
            ximode = 0;
            xi = 0.01;
        }
        else    {
            ximode = 3;
        }

        gc = ingc;
        if (gc == -1)   {
            gcmode = 0;
            gc = 0.5;
        }
        else    {
            gcmode = 3;
        }

        ReadFiles(datafile, treefile);
    }

    SparseAAMutSelModel(const SparseAAMutSelModel& from) = delete;

    ~SparseAAMutSelModel() {}

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

        // mutation rate parameters
        kappa = 2.0;
        // xi = 0.01;
        // gc = 0.4;
        nucmatrix = new T92SubMatrix(kappa,gc,true);
        codonstatespace = new CodonStateSpace(Universal);

        profileshape = 20.0;
        profilecenter.assign(Naa,1.0/Naa);
        preprofile = new IIDMultiGamma(Nsite,Naa,profileshape,profilecenter);
        if (profilemode == 3)   {
            preprofile->SetUniform();
        }

        // pi = 0.1;
        // sitemaskarray = new IIDProfileMask(Nsite,Naa,pi);
        pi.assign(Naa,0.1);
        sitemaskarray = new ProfileMask(Nsite,pi);
        FitSiteMaskArray();

        profile = new MutSelSparseFitnessArray(*preprofile,*sitemaskarray,maskepsilon);
        
        // mut sel matrices (based on the profiles of the mixture)
        sitematrixarray = new AAMutSelSubMatrixArray(codonstatespace, nucmatrix, xi, profile, false);

		phyloprocess = new PhyloProcess(tree,data,branchlength,siterate,sitematrixarray);
		phyloprocess->Unfold();

        // create suffstat arrays
		sitepathsuffstatarray = new PathSuffStatArray(Nsite);
    }

    void SetBLMode(int in)   {
        blmode = in;
    }

    void SetProfileMode(int in)   {
        profilemode = in;
    }

    void SetRateMode(int in)   {
        ratemode = in;
    }

    void SetMaskMode(int in)    {
        maskmode = in;
    }

    void SetMaskEpsilonMode(int in)    {
        maskepsilonmode = in;
    }

    void FitSiteMaskArray() {
        for (int i=0; i<Nsite; i++)    {
            vector<int>& x = (*sitemaskarray)[i];
            for (int k=0; k<Naa; k++)   {
                x[k] = 0;
            }
            for (int j=0; j<Ntaxa; j++) {
                if (data->GetState(j,i) != unknown) {
                    x[data->GetState(j,i)] = 1;
                }
            }
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
        nucmatrix->SetKappa(kappa);
        nucmatrix->SetGC(gc);
        nucmatrix->CorruptMatrix();
        sitematrixarray->SetXi(xi);
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
        // sitemaskarray->SetPi(pi);
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
        total += NucRateLogPrior();
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

    double NucRateLogPrior() const {
        double ret = -kappa/10;
        if (ximode < 2) {
           ret -= 10*xi;
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

    //! \brief log prob factor to be recomputed when moving kappa and gc
    double NucRateLogProb() const   {
        return NucRateLogPrior() + SuffStatLogProb();
    }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
	double Move() override {
        ResampleSub(1.0);
        MoveParameters(3,1);
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
                if (profilemode < 2)    {
                    MovePreProfile(10);
                    CompMovePreProfile(3);
                }
                if (maskmode < 3)   {
                    MoveMasks(20);
                }
                if (maskmode < 2)   {
                    MoveMaskHyperParameters();
                }
                if (maskepsilonmode < 2)    {
                    MoveMaskEpsilon(10);
                }
            }
            MoveNucRate();
            CorruptMatrices();
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
        ScalingMove(lambda,1.0,10,&SparseAAMutSelModel::BranchLengthsHyperLogProb,&SparseAAMutSelModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&SparseAAMutSelModel::BranchLengthsHyperLogProb,&SparseAAMutSelModel::NoUpdate,this);
        blhypermean->SetAllBranches(1.0/lambda);
	}

    //! MH move on nuc rate parameters: kappa and gc
	void MoveNucRate()	{
        ScalingMove(kappa,1.0,3,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
        ScalingMove(kappa,0.3,3,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
        ScalingMove(kappa,0.1,3,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
        if (gcmode < 2) {
            SlidingMove(gc,1.0,3,0,1.0,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
            SlidingMove(gc,0.3,3,0,1.0,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
            SlidingMove(gc,0.1,3,0,1.0,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
        }
        if (ximode < 2) {
            ScalingMove(xi,1.0,3,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
            ScalingMove(xi,0.3,3,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
            ScalingMove(xi,0.1,3,&SparseAAMutSelModel::NucRateLogProb,&SparseAAMutSelModel::CorruptMatrices,this);
        }
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
        ScalingMove(alpha,1.0,10,&SparseAAMutSelModel::SiteRateHyperLogProb,&SparseAAMutSelModel::NoUpdate,this);
        ScalingMove(alpha,0.3,10,&SparseAAMutSelModel::SiteRateHyperLogProb,&SparseAAMutSelModel::NoUpdate,this);
        siterate->SetShape(alpha);
        siterate->SetScale(alpha);
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

                double deltalogprob = 0;

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
            MovePreProfileAll(1.0, 1, 10);
            MovePreProfileAll(1.0, 3, 10);
            MovePreProfileAll(1.0, 20, 10);
            MovePreProfileAll(0.3, 20, 10);
        }
        // if masks are activated, move all active entries
        else    {
            MovePreProfile(1.0, 10);
            MovePreProfile(0.3, 10);
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
        for (int k=0; k<Naa; k++)   {
            SlidingMove(pi[k],1.0,10,0.05,0.975,&SparseAAMutSelModel::MaskLogProb,&SparseAAMutSelModel::UpdateMask,this);
            SlidingMove(pi[k],0.1,10,0.05,0.975,&SparseAAMutSelModel::MaskLogProb,&SparseAAMutSelModel::UpdateMask,this);
        }
    }

    //! MH move schedule on background fitness (maskepsilon)
    void MoveMaskEpsilon(int nrep)  {
        SlidingMove(maskepsilon,1.0,nrep,0,1.0,&SparseAAMutSelModel::MaskEpsilonLogProb,&SparseAAMutSelModel::UpdateAll,this);
        SlidingMove(maskepsilon,0.1,nrep,0,1.0,&SparseAAMutSelModel::MaskEpsilonLogProb,&SparseAAMutSelModel::UpdateAll,this);
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
                    if ((profilemode < 2) && mask[k])    {
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
        // os << "pi\t";
        os << "meanpi\t";
        os << "varpi\t";
        os << "width\t";
        os << "epsilon\t";
        os << "kappa\t";
        if (ximode < 2) {
            os << "xi\t";
        }
        if (gcmode < 2) {
            os << "gc\t";
        }
        os << '\n';
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
        // os << pi << '\t';
        os << GetMeanPi() << '\t';
        os << GetVarPi() << '\t';
        os << sitemaskarray->GetMeanWidth() << '\t';
        os << maskepsilon << '\t';
        os << kappa << '\t';
        if (ximode < 2) {
            os << xi << '\t';
        }
        if (gcmode < 2) {
            os << gc << '\t';
        }
        os << '\n';
    }

    double GetMeanPi() const    {
        double mean = 0;
        for (int k=0; k<Naa; k++)   {
            mean += pi[k];
        }
        mean /= Naa;
        return mean;
    }

    double GetVarPi() const {
        double mean = 0;
        double var = 0;
        for (int k=0; k<Naa; k++)   {
            mean += pi[k];
            var += pi[k]*pi[k];
        }
        mean /= Naa;
        var /= Naa;
        var -= mean*mean;
        return var;
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

    //! monitoring MCMC statistics
    void Monitor(ostream&) const override {}

    //! get complete parameter configuration from stream
    void FromStream(istream& is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        is >> kappa;
        if (ximode < 2) {
            is >> xi;
        }
        if (gcmode < 2) {
            is >> gc;
        }
        if (profilemode < 2)    {
            is >> *preprofile;
        }
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
        os << kappa << '\t';
        if (ximode < 2) {
            os << xi << '\t';
        }
        if (gcmode < 2) {
            os << gc << '\t';
        }
        if (profilemode < 2)    {
            os << *preprofile << '\t';
        }
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
        // mut rates
        size += 3;
        if (profilemode < 2)    {
            size += preprofile->GetMPISize();
        }
        if (maskmode < 2)   {
            // size++;
            size += Naa;
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
        is >> kappa;
        if (ximode < 2) {
            is >> xi;
        }
        if (gcmode < 2) {
            is >> gc;
        }
        if (profilemode < 2)    {
            is >> *preprofile;
        }
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
        os << kappa;
        if (ximode < 2) {
            os << xi;
        }
        if (gcmode < 2) {
            os << gc;
        }
        if (profilemode < 2)    {
            os << *preprofile;
        }
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

