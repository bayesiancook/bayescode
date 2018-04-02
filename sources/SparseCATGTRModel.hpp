
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
#include "GTRSubMatrixArray.hpp"
#include "PathSuffStat.hpp"
#include "MultiGammaSuffStat.hpp"

class SparseCATGTRModel : public ProbModel    {

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

    // aa exchange rates
    vector<double> relratehypercenter;
    double relratehyperinvconc;
	std::vector<double> relrate;

    double shape;
    vector<double> center;
    IIDMultiGamma* preprofile;

    double maskepsilon;
    double pi;
    IIDProfileMask* sitemaskarray;

    // across conditions and across sites
    MutSelSparseFitnessArray* profile;

    // an array of site-specific matrices
	GTRSubMatrixArray* sitematrixarray;

    // phyloprocess
    PhyloProcess* phyloprocess;

    // suff stats
	PathSuffStatArray* sitepathsuffstatarray;
    MultiGammaSuffStat hyperprofilesuffstat;

    int blmode;
    int rrmode;
    int maskepsilonmode;
    int maskmode;

    public:

    SparseCATGTRModel(const string& datafile, const string& treefile, double inepsilon) : hyperprofilesuffstat(Naa) {

        blmode = 0;
        rrmode = 0;

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

    SparseCATGTRModel(const SparseCATGTRModel& from) = delete;

    ~SparseCATGTRModel() {}

    //! read data and tree files
    void ReadFiles(string datafile, string treefile) {

        // sequence alignment
        data = new FileSequenceAlignment(datafile);
        Nsite = data->GetNsite();  // # columns
        Ntaxa = data->GetNtaxa();
        taxonset = data->GetTaxonSet();

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

        // relative exchange rates
        int Nrr = Naa*(Naa-1)/2; 
        relratehypercenter.assign(Nrr,1.0/Nrr);
        relratehyperinvconc = 1.0 / Nrr;
		relrate.assign(Nrr,0);
        Random::DirichletSample(relrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

        shape = 20.0;
        center.assign(Naa,1.0/Naa);
        preprofile = new IIDMultiGamma(Nsite,Naa,shape,center);

        pi = 0.1;
        sitemaskarray = new IIDProfileMask(Nsite,Naa,pi);

        profile = new MutSelSparseFitnessArray(*preprofile,*sitemaskarray,maskepsilon);
        
        // mut sel matrices (based on the profiles of the mixture)
        sitematrixarray = new GTRSubMatrixArray(relrate, profile,false);

		phyloprocess = new PhyloProcess(tree,data,branchlength,0,sitematrixarray);
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
    void SetRelRatesHyperParameters(const std::vector<double>& inrelratehypercenter, double inrelratehyperinvconc)  {
        relratehypercenter = inrelratehypercenter;
        relratehyperinvconc = inrelratehyperinvconc;
    }

    //! set relative exchange rates to a new value (multi-gene analyses)
    void SetRelRates(const std::vector<double>& inrelrate)  {
        relrate = inrelrate;
        CorruptMatrices();
    }

    //! copy relative exchange rates into vectors given as arguments (multi-gene analyses)
    void GetRelRates(std::vector<double>& inrelrate) const  {
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
        UpdateMask();
		preprofile->SetShape(shape);
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
        if (rrmode < 2)    {
            total += RelRateLogPrior();
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

    //! log prior over nuc rates rho and pi (uniform)
    double RelRateLogPrior() const {
        return 0;
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

    // ---------------
    // log probs for MH moves
    // ---------------

    //! \brief log prob factor to be recomputed when moving branch lengths hyperparameters (here, lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation rate parameters (nucrelrate and nucstat)
    double RelRateLogProb() const {
        return RelRateLogPrior() + SuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving mask hyperparameter pi
    double MaskLogProb() const  {
        return MaskHyperLogPrior() + MaskLogPrior();
    }

    //! \brief log prob factor to be recomputed when moving maskepsilon
    double MaskEpsilonLogProb() const  {
        return MaskHyperLogPrior() + SuffStatLogProb();
    }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
	double Move() override {
        ResampleSub(1.0);
        MoveParameters(3,20);
        return 1.0;
	}

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep0, int nrep) {

        for (int rep0 = 0; rep0 < nrep0; rep0++) {
            if (blmode < 2)    {
                MoveBranchLengths();
            }
            CollectSitePathSuffStat();
            UpdateAll();
            for (int rep = 0; rep < nrep; rep++) {
                MovePreProfile();
                CompMovePreProfile();
                if (maskmode < 3)   {
                    MoveMasks();
                }
                if (maskmode < 2)   {
                    MoveMaskHyperParameters();
                }
                if (maskepsilonmode < 2)    {
                    MoveMaskEpsilon();
                }
            }
            if (rrmode < 2)    {
                MoveRelRates();
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
        ScalingMove(lambda,1.0,10,&SparseCATGTRModel::BranchLengthsHyperLogProb,&SparseCATGTRModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&SparseCATGTRModel::BranchLengthsHyperLogProb,&SparseCATGTRModel::NoUpdate,this);
        blhypermean->SetAllBranches(1.0/lambda);
	}

    //! MH moves on relative exchange rates
	void MoveRelRates()	{

        CorruptMatrices();

        ProfileMove(relrate,0.1,1,10,&SparseCATGTRModel::RelRateLogProb,&SparseCATGTRModel::CorruptMatrices,this);
        ProfileMove(relrate,0.03,3,10,&SparseCATGTRModel::RelRateLogProb,&SparseCATGTRModel::CorruptMatrices,this);
        ProfileMove(relrate,0.01,3,10,&SparseCATGTRModel::RelRateLogProb,&SparseCATGTRModel::CorruptMatrices,this);
        CorruptMatrices();
	}

    //! MH compensatory move schedule on parameters and hyper-parameters
    void CompMovePreProfile()  {
        CompMovePreProfile(1.0,10);
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
                        double alpha = shape*center[a];
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
                        double alpha = shape*center[a];
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
    void MovePreProfile() {
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
        SlidingMove(pi,1.0,10,0.05,0.975,&SparseCATGTRModel::MaskLogProb,&SparseCATGTRModel::UpdateMask,this);
        SlidingMove(pi,0.1,10,0.05,0.975,&SparseCATGTRModel::MaskLogProb,&SparseCATGTRModel::UpdateMask,this);
    }

    //! MH move schedule on background fitness (maskepsilon)
    void MoveMaskEpsilon()  {
        SlidingMove(maskepsilon,1.0,10,0,1.0,&SparseCATGTRModel::MaskEpsilonLogProb,&SparseCATGTRModel::UpdateAll,this);
        SlidingMove(maskepsilon,0.1,10,0,1.0,&SparseCATGTRModel::MaskEpsilonLogProb,&SparseCATGTRModel::UpdateAll,this);
    }

    //! MH move on masks across sites
    double MoveMasks()    {
		double nacc = 0;
		double ntot = 0;
        for (int i=0; i<Nsite; i++) {
            vector<int>& mask = (*sitemaskarray)[i];
            int naa = 0;
            for (int k=0; k<Naa; k++)   {
                naa += mask[k];
            }
            for (int k=0; k<Naa; k++)   {
                if ((!mask[k]) || (naa > 1))    {
                    double deltalogprob = -MaskLogPrior(i) - SiteSuffStatLogProb(i);
                    naa -= mask[k];
                    mask[k] = 1-mask[k];
                    naa += mask[k];
                    if (mask[k])    {
                        (*preprofile)[i][k] = Random::sGamma(shape * center[k]);
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
        os << "pi\t";
        os << "width\t";
        os << "epsilon\t";
        os << "rrent\n";
    }

    //! write trace (one line summarizing current state) into trace file
    void Trace(ostream& os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << pi << '\t';
        os << sitemaskarray->GetMeanWidth() << '\t';
        os << maskepsilon << '\t';
        os << Random::GetEntropy(relrate) << '\n';
    }

    //! monitoring MCMC statistics
    void Monitor(ostream&) const override {}

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

