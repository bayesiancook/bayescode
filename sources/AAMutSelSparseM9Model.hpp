/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent evolution using Bayesian
phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/


#include "CodonSequenceAlignment.hpp"
#include "GTRSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "IIDMultiGamma.hpp"
#include "IIDMultiBernoulli.hpp"
#include "MultiGammaSuffStat.hpp"
#include "DiffSelSparseFitnessArray.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "CodonSuffStat.hpp"
#include "SubMatrixSelector.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "PathSuffStat.hpp"
#include "IIDProfileMask.hpp"
#include "IIDMechM9.hpp"

class AAMutSelSparseM9Model : public ProbModel {

    // -----
    // model selectors
    // -----

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    int blmode;
    int nucmode;
    int maskepsilonmode;
    int maskmode;
    int fitnesshypermode;
    int fixaa;

    // -----
    // external parameters
    // -----

    const Tree* tree;
    FileSequenceAlignment* data;
    const TaxonSet* taxonset;
    const CodonSequenceAlignment* codondata;

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

    // nucleotide exchange rates and equilibrium frequencies (stationary probabilities)
    // hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    // parameters
	std::vector<double> nucrelrate;
	std::vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    // Omega

    // weight of positive selection component
    double pospi;
    double poswhypermean;
    double poswhyperinvconc;
	double posw;

    double posmeanhypermean;
    double posmeanhyperinvshape;
    double posmean;

    double posinvshapehypermean;
    double posinvshapehyperinvshape;
    double posinvshape;

    IIDMechM9 *omegaarray;

    OmegaPathSuffStatArray *omegapathsuffstatarray;
    MechM9SuffStat omegahypersuffstat;

    double fitnessshape;
    vector<double> fitnesscenter;
    IIDMultiGamma* fitness;

    double maskepsilon;
    double pi;
    IIDProfileMask* sitemaskarray;

    // across conditions and across sites
    MutSelSparseFitnessArray* fitnessprofile;

    // an array of site-specific codon matrices
	AAMutSelOmegaCodonSubMatrixArray* sitecodonmatrixarray;

    // phyloprocess
    PhyloProcess* phyloprocess;

    // suff stats
	PathSuffStatArray* sitepathsuffstatarray;
    MultiGammaSuffStat hyperfitnesssuffstat;

    int chainsize;
    int burnin;

  public:

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions, such as specified by branch names)
    AAMutSelSparseM9Model(const std::string& datafile, const std::string& treefile, int infitnesshypermode, int infixaa, double inepsilon, double inpi, double inpospi) : hyperfitnesssuffstat(Naa) {

        blmode = 0;
        nucmode = 0;
        fitnesshypermode = infitnesshypermode;
        fixaa = infixaa;

        chainsize = 0;
        burnin = 20;

        //
        // maskmode = 3: no masks: aa fitness ~ iid uniform
        // maskmode = 2: fixed pi
        // maskmode = 0: pi estimated
        //
        // maskepsilonmode = 3: free epsilon
        // maskepsilonmode = 0: fixed epsilon
        //
        // fitnesshypermode = 3: fixed hyperparameters for fitness profiles
        // fitnesshypermode = 0:: hyperparameters estimated
        //
        // fixaa: flat fitness profiles
        //

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
        }

        pi = 0.1;
        if (maskmode < 2)   {
            if (inpi != -1.0) {
                maskmode = 2;
                pi = inpi;
            }
        }

        pospi = inpospi;

        poswhypermean = 0.5;
        poswhyperinvconc = 0.1;
        posmeanhypermean = 1.0;
        posmeanhyperinvshape = 1.0;
        posinvshapehypermean = 1.0;
        posinvshapehyperinvshape = 1.0;

        ReadFiles(datafile, treefile);
    }

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions, such as specified by branch names)
    AAMutSelSparseM9Model(const CodonSequenceAlignment* incodondata, const Tree* intree, int infitnesshypermode, int infixaa, double inepsilon, double inpi, double inpospi) : hyperfitnesssuffstat(Naa) {

        blmode = 0;
        nucmode = 0;
        fitnesshypermode = infitnesshypermode;
        fixaa = infixaa;

        chainsize = 0;
        burnin = 20;

        //
        // maskmode = 3: no masks: aa fitness ~ iid uniform
        // maskmode = 2: fixed pi
        // maskmode = 0: pi estimated
        //
        // maskepsilonmode = 3: free epsilon
        // maskepsilonmode = 0: fixed epsilon
        //
        // fitnesshypermode = 3: fixed hyperparameters for fitness profiles
        // fitnesshypermode = 0:: hyperparameters estimated
        //
        // fixaa: flat fitness profiles
        //

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
        }

        pi = 0.1;
        if (maskmode < 2)   {
            if (inpi != -1.0) {
                maskmode = 2;
                pi = inpi;
            }
        }
        codondata = incodondata;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();

        pospi = inpospi;

        poswhypermean = 0.5;
        poswhyperinvconc = 0.1;
        posmeanhypermean = 1.0;
        posmeanhyperinvshape = 1.0;
        posinvshapehypermean = 1.0;
        posinvshapehyperinvshape = 1.0;
    }

    void SetChainSize(int insize)    {
        chainsize = insize;
    }

    void SetBurnin(int inburnin)    {
        burnin = inburnin;
    }

    // AAMutSelSparseM9Model(const AAMutSelSparseM9Model&) = delete;

    ~AAMutSelSparseM9Model() {}

    //! read files (and read out the distribution of conditions across branches, based on the tree read from treefile)
    void ReadFiles(string datafile, string treefile) {
        // nucleotide sequence alignment
        data = new FileSequenceAlignment(datafile);

        // translated into codon sequence alignment
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        std::cerr << "-- Number of sites: " << Nsite << std::endl;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        Tree* tmptree = new Tree(treefile);
        // check whether tree and data fits together
        tmptree->RegisterWith(taxonset);
        tmptree->SetIndices();
        tree = tmptree;

        // traversal of the tree, so as to number links, branches and nodes
        // convention is: branches start at 1 (branch number 0 is the null branch behind the root)
        // nodes start at 0 (for the root), and nodes 1..Ntaxa are tip nodes (corresponding to taxa
        // in sequence alignment)
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

        nucrelratehypercenter.assign(Nrr,1.0/Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc,1.0/Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix
		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));
		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));
		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        // omega
        posw = poswhypermean;
        posmean = posmeanhypermean;
        posinvshape = posinvshapehypermean;

        omegaarray = new IIDMechM9(Nsite, posw, posmean, posinvshape);
        omegapathsuffstatarray = new OmegaPathSuffStatArray(Nsite);

        fitnessshape = 20.0;
        fitnesscenter.assign(Naa,1.0/Naa);
        fitness = new IIDMultiGamma(Nsite,Naa,fitnessshape,fitnesscenter);

        // pi = 0.1;
        sitemaskarray = new IIDProfileMask(Nsite,Naa,pi);

        fitnessprofile = new MutSelSparseFitnessArray(*fitness,*sitemaskarray,maskepsilon);
        
        // mut sel codon matrices (based on the fitness profiles of the mixture)
        sitecodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, fitnessprofile, omegaarray);

		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitecodonmatrixarray);
		phyloprocess->Unfold();

        // create suffstat arrays
		sitepathsuffstatarray = new PathSuffStatArray(Nsite);
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
	const CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //! \brief set estimation method for branch lengths
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetBLMode(int in)   {
        blmode = in;
    }

    //! \brief set estimation method for nuc rates
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetNucMode(int in) {
        nucmode = in;
    }

    //! \brief set estimation method for fitness hyperparameters
    //!
    //! thus far, mask model gives reasonable and interesting results only with fixed hyper params
    void SetFitnessHyperMode(int in)    {
        fitnesshypermode = in;
    }

    void SetFixAA(int in)   {
        fixaa = in;
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! \brief set estimation method for site profile masks
    //!
    //! Used in a multigene context.
    //! - mode == 3: no mask
    //! - mode == 2: parameter (maskprob) shared across genes
    //! - mode == 1: gene-specific parameter (maskprob), hyperparameters estimated across genes
    //! - mode == 0: gene-specific parameter (maskprob) with fixed hyperparameters
    void SetMaskMode(int in)    {
        maskmode = in;
    }

    //! \brief set estimation method for background fitness (maskepsilon)
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

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc, const std::vector<double>& innucstathypercenter, double innucstathyperinvconc) {
        nucrelratehypercenter = innucrelratehypercenter;
        nucrelratehyperinvconc = innucrelratehyperinvconc;
        nucstathypercenter = innucstathypercenter;
        nucstathyperinvconc = innucstathyperinvconc;
    }

    //! set nucleotide rates to a new value (multi-gene analyses)
    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        CorruptMatrices();
    }

    //! copy nucleotide rates into vectors given as arguments (multi-gene analyses)
    void GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! \brief set value of background fitness of low-fitness amino-acids
    void SetEpsilon(double in)  {
        maskepsilon = in;
        fitnessprofile->SetEpsilon(maskepsilon);
    }

    double GetEpsilon() const   {
        return maskepsilon;
    }

    void SetPi(double inpi) {
        pi = inpi;
        sitemaskarray->SetPi(pi);
    }

    double GetPi() const    {
        return pi;
    }

    void GetMixtureParameters(double& inposw, double& inposmean, double& inposinvshape) const   {
        inposw = posw;
        inposmean = posmean;
        inposinvshape = posinvshape;
    }

    void SetMixtureParameters(double inposw, double inposmean, double inposinvshape)    {
        posw = inposw;
        posmean = inposmean;
        posinvshape = inposinvshape;
        UpdateOmega();
    }

    void SetMixtureHyperParameters(
            double inpospi, 
            double inposmeanhypermean, double inposmeanhyperinvshape,
            double inposinvshapehypermean, double inposinvshapehyperinvshape,
            double inposwhypermean, double inposwhyperinvconc)  {

        pospi = inpospi;

        posmeanhypermean = inposmeanhypermean;
        posmeanhyperinvshape = inposmeanhyperinvshape;
        posinvshapehypermean = inposinvshapehypermean;
        posinvshapehyperinvshape = inposinvshapehyperinvshape;

        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;

        if (!pospi) {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }

        if (! posmeanhyperinvshape)  {
            posmean = posmeanhypermean;
        }

        if (! posinvshapehyperinvshape) {
            posinvshape = posinvshapehypermean;
        }

        if (! poswhyperinvconc) {
            posw = poswhypermean;
        }
    }

    void Update() override {
        if (blmode == 0)    {
            blhypermean->SetAllBranches(1.0/lambda);
        }
        UpdateMask();
		fitness->SetShape(fitnessshape);
        UpdateOmega();
        UpdateAll();
        ResampleSub(1.0);
    }

    void PostPred(string name) override {
        if (blmode == 0)    {
            blhypermean->SetAllBranches(1.0/lambda);
        }
        UpdateMask();
		fitness->SetShape(fitnessshape);
        UpdateOmega();
        UpdateAll();
        phyloprocess->PostPredSample(name);
    }

    void UpdateOmega()  {
        omegaarray->SetParameters(posw, posmean, posinvshape);
    }

    //! \brief dummy function that does not do anything.
    //! 
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove (defined in ProbModel),
    //! all of which require a void (*f)(void) function pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief tell the nucleotide and the codon matrices that their parameters have changed and that it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is notified,
    //! such that the matrices know that they will have to recalculate whichever component is requested later on upon demand.
    void CorruptMatrices()  {
        CorruptNucMatrix();
        CorruptCodonMatrices();
    }

    //! \brief tell the codon matrices that their parameters have changed and that it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is notified,
    //! such that the matrices know that they will have to recalculate whichever component is requested later on upon demand.
    void CorruptCodonMatrices() {
        sitecodonmatrixarray->UpdateCodonMatrices();
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is notified,
    //! such that the matrix knows that it will have to recalculate whichever component is requested later on upon demand.
    void CorruptNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! update fitness profiles and matrices across all sites and conditions
    void UpdateAll() {
        fitnessprofile->SetEpsilon(maskepsilon);
        fitnessprofile->Update();
        CorruptMatrices();
    }

    //! update fitness profiles and matrices across all conditions for site i
    void UpdateSite(int i) {
        fitnessprofile->Update(i);
        (*sitecodonmatrixarray)[i].CorruptMatrix();
    }

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
        if (nucmode < 2)    {
            total += NucRatesLogPrior();
        }
        if (fitnesshypermode < 2)   {
            total += FitnessHyperLogPrior();
        }
        if (! fixaa)    {
            total += FitnessLogPrior();
        }
        if (maskmode < 3)   {
            total += MaskLogPrior();
        }
        if (maskmode < 2)   {
            total += MaskHyperLogPrior();
        }
        total += OmegaHyperLogPrior();
        total += OmegaLogPrior();
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
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate,nucrelratehypercenter,1.0/nucrelratehyperinvconc);
        total += Random::logDirichletDensity(nucstat,nucstathypercenter,1.0/nucstathyperinvconc);
        return total;
    }

    //! log prior over fitness hyperparameters
    double FitnessHyperLogPrior() const {
        // uniform on center
        // exponential on shape
        return -fitnessshape;
    }

    //! log prior over input fitness parameters
    double FitnessLogPrior() const  {
        return fitness->GetLogProb();
    }

    double FitnessLogPrior(int i) const {
        return fitness->GetLogProb(i);
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

    //! log prior over mask array
    double MaskLogPrior(int i) const   {
        return sitemaskarray->GetLogProb(i);
    }

    // Omega
    //! log prior over omega mixture
    double OmegaHyperLogPrior() const    {
        double total = 0;
        total += PosOmegaLogPrior();
        total += PosWeightLogPrior();
        if (std::isinf(total))  {
            cerr << "omega hyper log prior: inf\n";
            exit(1);
        }
        return total;
    }

    //! gamma prior for dposom
    double PosOmegaLogPrior() const {
        double total = 0;
        if (posmeanhyperinvshape)   {
            double shape = 1.0 / posmeanhyperinvshape;
            double scale = shape / posmeanhypermean;
            total += Random::logGammaDensity(posmean, shape, scale);
        }
        if (posinvshapehyperinvshape)   {
            double shape = 1.0 / posinvshapehyperinvshape;
            double scale = shape / posinvshapehypermean;
            total += Random::logGammaDensity(posinvshape, shape, scale);
        }
        return total;
    }

    //! mixture of point mass at 0 (with prob pi) and Beta distribution (with prob
    //! 1 - pi) for posw
    double PosWeightLogPrior() const    {
        if (posw) {
            if (!pi) {
                cerr << "in PosWeightLogProb: pi == 0 and posw > 0\n";
                exit(1);
            }

            if (! poswhyperinvconc)   {
                return 0;
            }
            double alpha = poswhypermean / poswhyperinvconc;
            double beta = (1 - poswhypermean) / poswhyperinvconc;
            return log(pi) + Random::logBetaDensity(posw, alpha, beta);
        } else {
            return log(1 - pi);
        }
    }

    double OmegaLogPrior() const { return omegaarray->GetLogProb(); }

    //! Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogPrior() const    {
        if (posw) {
            return log(pi);
        }
        return log(1 - pi);
    }

    //! return log likelihood
    double GetLogLikelihood() const { 
        return phyloprocess->GetLogLikelihood();
    }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override {
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

    //! return log prob of the current substitution mapping, as a function of the current codon substitution process
	double SuffStatLogProb() const {
        return sitepathsuffstatarray->GetLogProb(*sitecodonmatrixarray);
	}

    //! return log prob of the substitution mappings for site i
    double SiteSuffStatLogProb(int i) const {
        return sitepathsuffstatarray->GetVal(i).GetLogProb(sitecodonmatrixarray->GetVal(i));
    }

    //! \brief return log prob of current branch lengths, as a function of branch lengths hyperparameter lambda
	double BranchLengthsHyperSuffStatLogProb() const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    //! return log prob of current fitness parameters, conditional on their hyperparameters
	double FitnessHyperSuffStatLogProb() const {
		double ret = hyperfitnesssuffstat.GetLogProb(fitnessshape,fitnesscenter);
        if (std::isinf(ret)) {
            cerr << "fitness hypersuffstat log prob is inf\n";
            exit(1);
        }
        return ret;
	}

    // Omega Hyper
    double OmegaHyperSuffStatLogProb() const {
        return omegahypersuffstat.GetLogProb(posw, posmean, posinvshape);
    }

    // ---------------
    // log probs for MH moves
    // ---------------

    //! \brief log prob factor to be recomputed when moving branch lengths hyperparameters (here, lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + SuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving fitness hyperparameters
    double FitnessHyperLogProb() const  {
        return FitnessHyperLogPrior() + FitnessHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving mask hyperparameter pi
    double MaskLogProb() const  {
        return MaskHyperLogPrior() + MaskLogPrior();
    }

    //! \brief log prob factor to be recomputed when moving maskepsilon
    double MaskEpsilonLogProb() const  {
        return MaskHyperLogPrior() + SuffStatLogProb();
    }

    // for moving omegamean and omegainvshape
    double OmegaHyperLogProb() const { return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb(); }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
	double Move() override {
        ResampleSub(1.0);
        MoveParameters(30,2);
        chainsize++;
        return 1.0;
	}

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep0, int nrep) {

        for (int rep0 = 0; rep0 < nrep0; rep0++) {
            if (blmode < 2)    {
                MoveBranchLengths();
            }
            CollectSitePathSuffStat();
            MoveAA(nrep);
            if (nucmode < 2)    {
                MoveNucRates();
            }
            MoveOmega();
            MoveOmegaHyperParameters();
        }

        UpdateAll();
    }

    void MoveAA(int nrep)   {
        UpdateAll();
        for (int rep = 0; rep < nrep; rep++) {
            if (! fixaa)    {
                MoveFitness();
                CompMoveFitness();
            }
            if (maskmode < 3)   {
                MoveMasks();
            }
            if (maskmode < 2)   {
                MoveMaskHyperParameters();
            }
            // works best when not used
            if (fitnesshypermode < 2)   {
                MoveFitnessHyperParameters();
            }
            if (maskepsilonmode < 2)    {
                MoveMaskEpsilon();
            }
        }
    }


    //! Gibbs resample substitution mappings conditional on current parameter configuration
    void ResampleSub(double frac)   {
        CorruptMatrices();
		phyloprocess->Move(frac);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current value of lambda)
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
        ScalingMove(lambda,1.0,10,&AAMutSelSparseM9Model::BranchLengthsHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelSparseM9Model::BranchLengthsHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);
        blhypermean->SetAllBranches(1.0/lambda);
	}

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using ProfileMove)
	void MoveNucRates()	{

        CorruptMatrices();

        ProfileMove(nucrelrate,0.1,1,10,&AAMutSelSparseM9Model::NucRatesLogProb,&AAMutSelSparseM9Model::CorruptMatrices,this);
        ProfileMove(nucrelrate,0.03,3,10,&AAMutSelSparseM9Model::NucRatesLogProb,&AAMutSelSparseM9Model::CorruptMatrices,this);
        ProfileMove(nucrelrate,0.01,3,10,&AAMutSelSparseM9Model::NucRatesLogProb,&AAMutSelSparseM9Model::CorruptMatrices,this);

        ProfileMove(nucstat,0.1,1,10,&AAMutSelSparseM9Model::NucRatesLogProb,&AAMutSelSparseM9Model::CorruptMatrices,this);
        ProfileMove(nucstat,0.01,1,10,&AAMutSelSparseM9Model::NucRatesLogProb,&AAMutSelSparseM9Model::CorruptMatrices,this);

        CorruptMatrices();
	}

    //! MH compensatory move on fitness parameters and hyper-parameters
    void CompMoveFitness()  {
        CompMoveFitness(1.0,10);
    }

    //! MH compensatory move on fitness parameters and hyper-parameters
    double CompMoveFitness(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*fitness)[i];
                const vector<int>& mask = (*sitemaskarray)[i];
                // vector<int> mask(Naa,1);
                

                double deltalogprob = 0;

                for (int a=0; a<Naa; a++)   {
                    if (mask[a])	{
                        double alpha = fitnessshape*fitnesscenter[a];
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
                        double alpha = fitnessshape*fitnesscenter[a];
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

    //! MH moves on baseline fitness parameters (for condition k=0)
    void MoveFitness() {
        // if masks are not activated (all entries equal to 1), move a random subset of entries over the 20 amino-acids (2d parameter of call)
        if (maskmode == 3)  {
            MoveFitnessAll(1.0, 1, 10);
            MoveFitnessAll(1.0, 3, 10);
            MoveFitnessAll(1.0, 20, 10);
            MoveFitnessAll(0.3, 20, 10);
        }
        // if masks are activated, move all active entries
        else    {
            MoveFitness(1.0, 10);
            MoveFitness(0.3, 10);
        }
    }

    //! MH moves on baseline fitness parameters (for condition k=0)
    double MoveFitness(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*fitness)[i];
                const vector<int>& s = (*sitemaskarray)[i];

                bk = x;

                double deltalogprob = -fitness->GetLogProb(i,s) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, s);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(i,s) + SiteSuffStatLogProb(i);

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

    //! MH moves on baseline fitness parameters (for condition k=0)
    double MoveFitnessAll(double tuning, int n, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*fitness)[i];

                bk = x;

                double deltalogprob = -fitness->GetLogProb(i) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, n);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(i) + SiteSuffStatLogProb(i);

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

    void MoveMaskHyperParameters()  {
        SlidingMove(pi,1.0,10,0.05,0.975,&AAMutSelSparseM9Model::MaskLogProb,&AAMutSelSparseM9Model::UpdateMask,this);
        SlidingMove(pi,0.1,10,0.05,0.975,&AAMutSelSparseM9Model::MaskLogProb,&AAMutSelSparseM9Model::UpdateMask,this);
    }

    void MoveMaskEpsilon()  {
        SlidingMove(maskepsilon,1.0,10,0,1.0,&AAMutSelSparseM9Model::MaskEpsilonLogProb,&AAMutSelSparseM9Model::UpdateAll,this);
        SlidingMove(maskepsilon,0.1,10,0,1.0,&AAMutSelSparseM9Model::MaskEpsilonLogProb,&AAMutSelSparseM9Model::UpdateAll,this);
    }

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
                        (*fitness)[i][k] = Random::sGamma(fitnessshape * fitnesscenter[k]);
                        if (! (*fitness)[i][k]) {
                            (*fitness)[i][k] = 1e-8;
                            // cerr << "null fitness : " << fitnessshape << '\t' << fitnesscenter[k] << '\n';
                            // exit(1);
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

    //! MH moves on hyperparameters of distribution of fitness factors
    void MoveFitnessHyperParameters() {
        // collect suff stats across all active fitness parameters
        hyperfitnesssuffstat.Clear();
        // hyperfitnesssuffstat.AddSuffStat(*fitness);
        hyperfitnesssuffstat.AddSuffStat(*fitness,*sitemaskarray);

        ScalingMove(fitnessshape,1.0,100,&AAMutSelSparseM9Model::FitnessHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);
        ScalingMove(fitnessshape,0.3,100,&AAMutSelSparseM9Model::FitnessHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);
        ScalingMove(fitnessshape,0.1,100,&AAMutSelSparseM9Model::FitnessHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);

        ProfileMove(fitnesscenter,0.3,1,100,&AAMutSelSparseM9Model::FitnessHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);
        ProfileMove(fitnesscenter,0.1,1,100,&AAMutSelSparseM9Model::FitnessHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);
        ProfileMove(fitnesscenter,0.1,3,100,&AAMutSelSparseM9Model::FitnessHyperLogProb,&AAMutSelSparseM9Model::NoUpdate,this);

		fitness->SetShape(fitnessshape);
        // fitness->PriorResample(*sitemaskarray,1e-8);
    }

    // Omega

    void MoveOmega() {
        omegapathsuffstatarray->Clear();
        omegapathsuffstatarray->AddSuffStat(*sitecodonmatrixarray, *sitepathsuffstatarray);
        omegaarray->MultipleTryMove(100, *omegapathsuffstatarray);
        CorruptCodonMatrices();
    }

    void MoveOmegaHyperParameters()    {

        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegaarray);

        if (pi != 0) {
            if (posmeanhyperinvshape)   {
                ScalingMove(posmean, 0.1, 10, &AAMutSelSparseM9Model::OmegaHyperLogProb, &AAMutSelSparseM9Model::NoUpdate, this);
            }
            if (posinvshapehyperinvshape)   {
                ScalingMove(posinvshape, 0.1, 10, &AAMutSelSparseM9Model::OmegaHyperLogProb, &AAMutSelSparseM9Model::NoUpdate, this);
            }
            if (poswhyperinvconc)   {
                SlidingMove(posw, 0.1, 10, 0, 1, &AAMutSelSparseM9Model::OmegaHyperLogProb, &AAMutSelSparseM9Model::NoUpdate, this);
            }
        }
        if ((pi != 0) && (pi != 1)) {
            SwitchPosWeight(10);
        }
        UpdateOmega();
    }

    //! reversible jump move on posw
    double SwitchPosWeight(int nrep)    {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double bkposw = posw;
            double deltalogprob = -PosSwitchLogPrior() - OmegaHyperLogProb();
            if (posw) {
                posw = 0;
            } else {
                double alpha = poswhypermean / poswhyperinvconc;
                double beta = (1 - poswhypermean) / poswhyperinvconc;
                posw = Random::BetaSample(alpha, beta);
            }
            deltalogprob += PosSwitchLogPrior() + OmegaHyperLogProb();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                posw = bkposw;
            }
            ntot++;
        }
        return nacc / ntot;
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    double GetEmpiricalPosFrac() const {
        double tot = 0;
        for (int i = 0; i < Nsite; i++) {
            if ((*omegaarray)[i] > 1.0) {
                tot++;
            }
        }
        return tot / Nsite;
    }

    void TraceOmega(ostream &os) const {
        for (int i = 0; i < GetNsite(); i++) {
            os << omegaarray->GetVal(i) << '\t';
        }
        os << '\n';
    }

    void GetSiteOmega(double *array) const {
        for (int i = 0; i < GetNsite(); i++) {
            array[i] = omegaarray->GetVal(i);
        }
    }

    double GetMeanOmega() const {
        return omegaarray->GetMean();
    }

    double GetMeanActiveFitness() const    {
        double mean = 0;
        int n = 0;
        for (int i=0; i<GetNsite(); i++)    {
            const vector<double>& f = fitness->GetVal(i);
            const vector<int>& m = sitemaskarray->GetVal(i);
            for (int a=0; a<Naa; a++)   {
                mean += m[a]*f[a];
                n += m[a];
            }
        }
        mean /= n;
        return mean;
    }

    double GetMeanAAEntropy() const {
        return fitnessprofile->GetMeanEntropy();
    }

    double GetMeanWidth() const {
        return sitemaskarray->GetMeanWidth();
    }

    void GetSitePredictedDNDS(double* array)  const   {
        for (int i=0; i<GetNsite(); i++)    {
            array[i] = sitecodonmatrixarray->GetVal(i).GetPredictedDNDS();
        }
    }

    double GetPredictedDNDS() const  {
        double mean = 0;
        for (int i=0; i<GetNsite(); i++) {
            mean += sitecodonmatrixarray->GetVal(i).GetPredictedDNDS();
        }
        mean /= Nsite;
        return mean;
    }

    void TraceHeader(ostream& os) const override {
        os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "dnds\t";
        os << "posfrac\t";
        os << "posw\t";
        os << "posmean\tinvshape\t";
        if (maskmode == 0) {
            os << "pi\t";
        }
        if (maskmode < 3)  {
            os << "width\t";
        }
        // os << "meanfitness\t";
        if (maskepsilonmode < 2)    {
            os << "epsilon\t";
        }
        os << "aaent\t";
        if (fitnesshypermode < 2)   {
            os << "shape\t";
            os << "center\t";
        }
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream& os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << 3*branchlength->GetTotalLength() << '\t';
		os << GetMeanOmega() << '\t';
        os << GetPredictedDNDS() << '\t';
        os << GetEmpiricalPosFrac() << '\t';
        os << posw << '\t';
        os << posmean << '\t' << posinvshape << '\t';
        if (maskmode == 0) {
            os << pi << '\t';
        }
        if (maskmode < 3)  {
            os << sitemaskarray->GetMeanWidth() << '\t';
        }
        // os << GetMeanActiveFitness() << '\t';
        if (maskepsilonmode < 2)    {
            os << maskepsilon << '\t';
        }
        os << GetMeanAAEntropy() << '\t';
        if (fitnesshypermode < 2)   {
            os << fitnessshape << '\t';
            os << Random::GetEntropy(fitnesscenter) << '\t';
        }
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream&) const override {}

    void FromStream(istream& is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2)    {
            is >> nucrelrate;
            is >> nucstat;
        }
        is >> posw;
        is >> posmean >> posinvshape;
        is >> *omegaarray;
        if (fitnesshypermode < 2)   {
            is >> fitnessshape;
            is >> fitnesscenter;
        }
        if (! fixaa)    {
            is >> *fitness;
        }
        if (maskmode == 0)   {
            is >> pi;
        }
        if (maskmode < 3)   {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            is >> maskepsilon;
        }
    }

    void ToStream(ostream& os) const override {
        if (blmode < 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        }
        if (nucmode < 2)    {
            os << nucrelrate << '\t';
            os << nucstat << '\t';
        }
        os << posw << '\t';
        os << posmean << '\t' << posinvshape << '\t';
        os << *omegaarray << '\t';
        if (fitnesshypermode < 2)   {
            os << fitnessshape << '\t';
            os << fitnesscenter << '\t';
        }
        if (! fixaa)    {
            os << *fitness << '\t';
        }
        if (maskmode == 0)   {
            os << pi << '\t';
        }
        if (maskmode < 3)   {
            os << *sitemaskarray << '\t';
        }
        if (maskepsilonmode < 2)    {
            os << maskepsilon << '\t';
        }
    }

    //! return size of model, when put into an MPI buffer (in multigene context -- only omegatree)
    unsigned int GetMPISize() const {
        int size = 0;
        if (blmode < 2) {
            size++;
            size += branchlength->GetMPISize();
        }
        if (nucmode < 2)    {
            size += nucrelrate.size();
            size += nucstat.size();
        }
        size += 3;
        size += GetNsite();
        if (fitnesshypermode < 2)   {
            size ++;
            size += fitnesscenter.size();
        }
        if (! fixaa)    {
            size += fitness->GetMPISize();
        }
        // pi and epsilon
        if (maskmode  == 0)   {
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

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& is)    {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2)    {
            is >> nucrelrate;
            is >> nucstat;
        }
        is >> posw;
        is >> posmean >> posinvshape;
        is >> *omegaarray;
        if (fitnesshypermode < 2)   {
            is >> fitnessshape;
            is >> fitnesscenter;
        }
        if (! fixaa)    {
            is >> *fitness;
        }
        if (maskmode  == 0)   {
            is >> pi;
        }
        if (maskmode < 3)   {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2)    {
            is >> maskepsilon;
        }
    }

    //! write array into MPI buffer
    void MPIPut(MPIBuffer& os) const {
        if (blmode < 2) {
            os << lambda;
            os << *branchlength;
        }
        if (nucmode < 2)    {
            os << nucrelrate;
            os << nucstat;
        }
        os << posw;
        os << posmean << posinvshape;
        os << *omegaarray;
        if (fitnesshypermode < 2)   {
            os << fitnessshape;
            os << fitnesscenter;
        }
        if (! fixaa)   {
            os << *fitness;
        }
        if (maskmode == 0)   {
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
