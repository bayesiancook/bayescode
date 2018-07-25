/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS)
(2017-06-14). Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent
evolution using Bayesian phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding
by the rules of distribution of free software. You can use, modify and/ or
redistribute the software under the terms of the CeCILL-C license as circulated
by CEA, CNRS and INRIA at the following URL "http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and
redistribute granted by the license, users are provided only with a limited
warranty and the software's author, the holder of the economic rights, and the
successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with
loading, using, modifying and/or developing or reproducing the software by the
user in light of its specific status of free software, that may mean that it is
complicated to manipulate, and that also therefore means that it is reserved for
developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as
regards their requirements in conditions enabling the security of their systems
and/or data to be ensured and, more generally, to use and operate it in the same
conditions as regards security.

The fact that you are presently reading this means that you have had knowledge
of the CeCILL-C license and that you accept its terms.*/

#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSuffStat.hpp"
#include "DiffSelSparseFitnessArray.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "IIDMultiBernoulli.hpp"
#include "IIDMultiGamma.hpp"
#include "IIDProfileMask.hpp"
#include "MultiGammaSuffStat.hpp"
#include "PathSuffStat.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "SubMatrixSelector.hpp"
#include "Tree.hpp"

/**
 * \brief The mutation-selection model with constant fitness landscape over the
 * tree -- sparse multi-gamma version.
 *
 * The model is parameterized by
 * - an fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for
 * j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process
 * (assumed homogeneous across sites and lineages)
 * - G_ia: an array of site-specific amino-acid (gamma-distributed) pre-fitness
 * vectors G_ia, for site i and amino-acid a
 * - d_ia: an array of site-specific masks over the 20 amino-acids (specifying
 * which are the low- and high-fitness amino-acids for that site)
 * - 0 < epsilon < 1: the background fitness (for those amino-acids that have
 * been masked)
 *
 * Based on this parameterization, the fitness of amino-acids is then
 * essentially determined by the gamma variates G_ia, except for those
 * amino-acids that have been masked (i.e. for which d_ia == 0), which are given
 * a low background fitness (equal to epsilon). Mathematically, the fitness of
 * amino-acid a at site i is defined by:
 * - F_ia = G_ia * d_ia + epsilon (1-d_ia)
 *
 * The model also implements an omega multiplier, capturing deviations of the
 * non-syn rate from the model (see Rodrigue and Lartillot, 2107); this
 * parameter is fixed to 1 by default.
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 * - omega*: fixed to 1 or exponential of rate 1
 * - (d_ia)_a=1..20: iid Bernoulli of probability pi, conditional on at least
 * one entry over the 20 amino-acids being equal to 1
 * - (G_ia)_a=1..20: multi-gamma: i.e., G_ia ~ Gamma(shape = fitnessshape *
 * fitnesscenter[a],scale = 1.0)
 * - pi: uniform
 * - epsilon: uniform
 * - fitnesscenter: uniform
 * - fitnessshape: exponential of mean 20
 */

class AAMutSelSparseOmegaModel : public ProbModel {
    // -----
    // model selectors
    // -----

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    int blmode;
    int nucmode;
    int omegamode;
    int maskepsilonmode;
    int maskmode;
    int fitnessshapemode;
    int fitnesscentermode;

    // -----
    // external parameters
    // -----

    Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    // number of sites
    int Nsite;
    int Ntaxa;
    int Nbranch;

    // -----
    //  model structure
    // -----

    // branch lengths
    double lambda;
    BranchIIDGamma *blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise *branchlength;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStat hyperlengthsuffstat;

    // nucleotide exchange rates and equilibrium frequencies (stationary
    // probabilities) hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    // parameters
    std::vector<double> nucrelrate;
    std::vector<double> nucstat;
    GTRSubMatrix *nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
    double omega;
    OmegaPathSuffStat omegapathsuffstat;

    double fitnessshape;
    vector<double> fitnesscenter;
    IIDMultiGamma *fitness;

    double maskepsilon;
    double pi;
    IIDProfileMask *sitemaskarray;

    // across conditions and across sites
    MutSelSparseFitnessArray *fitnessprofile;

    // an array of site-specific codon matrices
    AAMutSelOmegaCodonSubMatrixArray *sitecodonmatrixarray;

    // phyloprocess
    PhyloProcess *phyloprocess;

    // suff stats
    PathSuffStatArray *sitepathsuffstatarray;
    MultiGammaSuffStat hyperfitnesssuffstat;

  public:
    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions,
    //! such as specified by branch names)
    //! - inomegamode: omega fixed (3), shared across genes (2) or estimated with
    //! shrinkage across genes (1) or without shrinkage (0)
    //! - inepsilon: background fitness for low-fitness amino-acids: if
    //! 0<inepsilon<1, then epsilon is fixed, if epsilon == 1, then this is the
    //! model without masks, if epsilon == -1, then epsilon is estimated from the
    //! data
    //! - inshape: shape parameter of the Gamma distribution of pre-fitness
    //! parameters. If inshape>0, shape parameter is fixed, if inshape == -1,
    //! shape parameter is estimated
    AAMutSelSparseOmegaModel(const std::string &datafile, const std::string &treefile,
                             int inomegamode, double inepsilon, double inshape)
        : hyperfitnesssuffstat(Naa) {
        blmode = 0;
        nucmode = 0;
        omegamode = inomegamode;
        fitnesscentermode = 3;
        fitnessshapemode = 3;
        fitnessshape = 20.0;

        if (inshape > 0) {
            fitnessshapemode = 3;
            fitnessshape = inshape;
        } else {
            fitnessshapemode = 0;
            fitnessshape = 20.0;
        }

        if (inepsilon == 1) {
            maskepsilon = 1;
            maskmode = 3;
            maskepsilonmode = 3;
        } else if (inepsilon >= 0) {
            maskepsilon = inepsilon;
            maskepsilonmode = 3;
            maskmode = 0;
        } else {
            maskepsilonmode = 0;
            maskmode = 0;
            maskepsilon = 0.01;
        }

        ReadFiles(datafile, treefile);
    }

    AAMutSelSparseOmegaModel(const AAMutSelSparseOmegaModel &) = delete;

    ~AAMutSelSparseOmegaModel() {}

    //! read data and tree files
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
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        // traversal of the tree, so as to number links, branches and nodes
        // convention is: branches start at 1 (branch number 0 is the null branch
        // behind the root) nodes start at 0 (for the root), and nodes 1..Ntaxa are
        // tip nodes (corresponding to taxa in sequence alignment)
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
        blhypermean = new BranchIIDGamma(*tree, 1.0, lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree, *blhypermean, 1.0 / blhyperinvshape);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));
        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        // global omega (fixed to 1 by default)
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        omega = 1.0;

        fitnessshape = 20.0;
        fitnesscenter.assign(Naa, 1.0 / Naa);
        fitness = new IIDMultiGamma(Nsite, Naa, fitnessshape, fitnesscenter);

        pi = 0.1;
        sitemaskarray = new IIDProfileMask(Nsite, Naa, pi);

        fitnessprofile = new MutSelSparseFitnessArray(*fitness, *sitemaskarray, maskepsilon);

        // mut sel codon matrices (based on the fitness profiles of the mixture)
        sitecodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix,
                                                                    fitnessprofile, omega);

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, sitecodonmatrixarray);
        phyloprocess->Unfold();

        // create suffstat arrays
        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
    }

    //! \brief set estimation method for branch lengths
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetBLMode(int in) { blmode = in; }

    //! \brief set estimation method for nuc rates
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetNucMode(int in) { nucmode = in; }

    //! \brief set estimation method for fitness hyperparameter (center of
    //! multi-gamma distribution)
    //!
    //! - mode == 3: fixed (uniform)
    //! - mode == 2: shared across genes, estimated
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetFitnessCenterMode(int in) { fitnesscentermode = in; }

    //! \brief set estimation method for fitness hyperparameter (shape of
    //! multi-gamma distribution)
    //!
    //! - mode == 3: fixed
    //! - mode == 2: shared across genes, estimated
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetFitnessShapeMode(int in, double inshape = 0) {
        fitnessshapemode = in;
        if (inshape > 0) {
            fitnessshape = inshape;
        }
    }

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 3: fixed to 1
    //! - mode == 2: shared and estimated across genes: currently not implemented
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, either mode 3 and mode 0 can be used -- default
    //! mode is 3.
    void SetOmegaMode(int mode) { omegamode = mode; }

    //! \brief set estimation method for site profile masks
    //!
    //! Used in a multigene context.
    //! - mode == 3: no mask
    //! - mode == 2: parameter (maskprob) shared across genes
    //! - mode == 1: gene-specific parameter (maskprob), hyperparameters estimated
    //! across genes
    //! - mode == 0: gene-specific parameter (maskprob) with fixed hyperparameters
    void SetMaskMode(int in) { maskmode = in; }

    //! \brief set estimation method for background fitness (maskepsilon)
    void SetMaskEpsilonMode(int in) { maskepsilonmode = in; }

    // ------------------
    // Update system
    // ------------------

    //! \brief set branch lengths to a new value
    //!
    //! Used in a multigene context.
    void SetBranchLengths(const BranchSelector<double> &inbranchlength) {
        branchlength->Copy(inbranchlength);
    }

    //! get a copy of branch lengths into array given as argument
    void GetBranchLengths(BranchArray<double> &inbranchlength) const {
        inbranchlength.Copy(*branchlength);
    }

    //! set branch lengths hyperparameters to a new value (multi-gene analyses)
    void SetBranchLengthsHyperParameters(const BranchSelector<double> &inblmean,
                                         double inblinvshape) {
        blhypermean->Copy(inblmean);
        blhyperinvshape = inblinvshape;
        branchlength->SetShape(1.0 / blhyperinvshape);
    }

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double> &innucrelratehypercenter,
                                    double innucrelratehyperinvconc,
                                    const std::vector<double> &innucstathypercenter,
                                    double innucstathyperinvconc) {
        nucrelratehypercenter = innucrelratehypercenter;
        nucrelratehyperinvconc = innucrelratehyperinvconc;
        nucstathypercenter = innucstathypercenter;
        nucstathyperinvconc = innucstathyperinvconc;
    }

    //! set nucleotide rates to a new value (multi-gene analyses)
    void SetNucRates(const std::vector<double> &innucrelrate,
                     const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        CorruptMatrices();
    }

    //! copy nucleotide rates into vectors given as arguments (multi-gene
    //! analyses)
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! return current omega value
    double GetOmega() const { return omega; }

    //! set omega to new value (multi-gene analyses)
    void SetOmega(double inomega) {
        omega = inomega;
        CorruptCodonMatrices();
    }

    //! set omega hyperparams to new value (multi-gene analyses)
    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape) {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    //! \brief set value of background fitness of low-fitness amino-acids
    void SetMaskEpsilon(double in) {
        maskepsilon = in;
        fitnessprofile->SetEpsilon(maskepsilon);
    }

    //! const ref access to masks across sites
    const vector<vector<int>> &GetMaskArray() const { return sitemaskarray->GetArray(); }

    //! const access to low-fitness background value (mask epsilon)
    double GetMaskEpsilon() const { return maskepsilon; }

    void Update() override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        UpdateMask();
        fitness->SetShape(fitnessshape);
        UpdateAll();
        ResampleSub(1.0);
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is
    //! notified, such that the matrices know that they will have to recalculate
    //! whichever component is requested later on upon demand.
    void CorruptMatrices() {
        CorruptNucMatrix();
        CorruptCodonMatrices();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is
    //! notified, such that the matrices know that they will have to recalculate
    //! whichever component is requested later on upon demand.
    void CorruptCodonMatrices() {
        sitecodonmatrixarray->SetOmega(omega);
        sitecodonmatrixarray->UpdateCodonMatrices();
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
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

    //! update parameters of sitemaskarray (probability parameter pi)
    void UpdateMask() { sitemaskarray->SetPi(pi); }

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
        if (nucmode < 2) {
            total += NucRatesLogPrior();
        }
        if ((fitnessshapemode < 2) || (fitnesscentermode < 2)) {
            total += FitnessHyperLogPrior();
        }
        // total += FitnessLogPrior();
        if (maskmode < 2) {
            total += MaskHyperLogPrior();
            total += MaskLogPrior();
        }
        if (omegamode < 2) {
            total += OmegaLogPrior();
        }
        return total;
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here,
    //! lambda ~ exponential of rate 10)
    double BranchLengthsHyperLogPrior() const { return -log(10.0) - lambda / 10; }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const {
        double ret = branchlength->GetLogProb();
        if (blmode == 0) {
            ret += BranchLengthsHyperLogPrior();
        }
        return ret;
    }

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate, nucrelratehypercenter,
                                             1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //! log prior over omega (gamma of mean omegahypermean and inverse shape
    //! omegahyperinvshape)
    double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(omega) -
               beta * omega;
    }

    //! log prior over fitness hyperparameters
    double FitnessHyperLogPrior() const {
        double ret = 0;
        if (fitnessshapemode < 2) {
            ret += FitnessShapeLogPrior();
        }
        if (fitnesscentermode < 2) {
            ret += FitnessCenterLogPrior();
        }
        return ret;
    }

    //! log prior over fitness center hyperparameter
    double FitnessCenterLogPrior() const { return 0; }

    //! log prior over fitness shape hyperparameter
    double FitnessShapeLogPrior() const {
        // exponential on shape
        return -fitnessshape;
    }

    //! log prior over input fitness parameters
    double FitnessLogPrior() const { return fitness->GetLogProb(); }

    //! log prior over input fitness parameters
    double FitnessLogPrior(int i) const { return fitness->GetLogProb(i); }

    //! log prior over mask array hyperparameters
    double MaskHyperLogPrior() const {
        double ret = 0;
        if (maskepsilonmode < 2) {
            ret -= 10 * maskepsilon;
        }
        return ret;
    }

    //! log prior over mask array
    double MaskLogPrior() const { return sitemaskarray->GetLogProb(); }

    //! log prior over mask array (only for site i)
    double MaskLogPrior(int i) const { return sitemaskarray->GetLogProb(i); }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    // ---------------
    // collecting suff stats
    // ---------------

    //! \brief const access to array of length-pathsuffstats across branches
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! collect sufficient statistics if substitution mappings across sites
    void CollectSitePathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double SuffStatLogProb() const {
        return sitepathsuffstatarray->GetLogProb(*sitecodonmatrixarray);
    }

    //! return log prob of the substitution mappings for site i
    double SiteSuffStatLogProb(int i) const {
        return sitepathsuffstatarray->GetVal(i).GetLogProb(sitecodonmatrixarray->GetVal(i));
    }

    //! \brief return log prob of current branch lengths, as a function of branch
    //! lengths hyperparameter lambda
    double BranchLengthsHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    //! return log prob of current fitness parameters, conditional on their
    //! hyperparameters
    double FitnessHyperSuffStatLogProb() const {
        double ret = hyperfitnesssuffstat.GetLogProb(fitnessshape, fitnesscenter);
        if (std::isinf(ret)) {
            cerr << "fitness hypersuffstat log prob is inf\n";
            exit(1);
        }
        return ret;
    }

    // ---------------
    // log probs for MH moves
    // ---------------

    //! \brief log prob factor to be recomputed when moving branch lengths
    //! hyperparameters (here, lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + SuffStatLogProb(); }

    //! \brief log prob factor to be recomputed when moving fitness
    //! hyperparameters
    double FitnessHyperLogProb() const {
        return FitnessHyperLogPrior() + FitnessHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving mask hyperparameter pi
    double MaskLogProb() const { return MaskHyperLogPrior() + MaskLogPrior(); }

    //! \brief log prob factor to be recomputed when moving maskepsilon
    double MaskEpsilonLogProb() const { return MaskHyperLogPrior() + SuffStatLogProb(); }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
    double Move() override {
        ResampleSub(1.0);
        MoveParameters(3, 20);
        return 1.0;
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep0, int nrep) {
        for (int rep0 = 0; rep0 < nrep0; rep0++) {
            if (blmode < 2) {
                MoveBranchLengths();
            }
            CollectSitePathSuffStat();
            UpdateAll();
            for (int rep = 0; rep < nrep; rep++) {
                MoveFitness();
                CompMoveFitness();
                if (maskmode < 3) {
                    MoveMasks();
                }
                if (maskmode < 2) {
                    MoveMaskHyperParameters();
                }
                if ((fitnessshapemode < 2) || (fitnesscentermode < 2)) {
                    MoveFitnessHyperParameters();
                }
                if (maskepsilonmode < 2) {
                    MoveMaskEpsilon();
                }
            }
            if (nucmode < 2) {
                MoveNucRates();
            }
            if (omegamode < 2) {
                MoveOmega();
            }
        }

        UpdateAll();
    }

    //! Gibbs resampling of substitution histories conditional on current
    //! parameter configuration
    void ResampleSub(double frac) {
        CorruptMatrices();
        phyloprocess->Move(frac);
    }

    //! Gibbs resampling of branch lengths (based on sufficient statistics and
    //! current value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MCMC move schedule on branch lengths
    void MoveBranchLengths() {
        ResampleBranchLengths();
        if (blmode == 0) {
            MoveLambda();
        }
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &AAMutSelSparseOmegaModel::BranchLengthsHyperLogProb,
                    &AAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &AAMutSelSparseOmegaModel::BranchLengthsHyperLogProb,
                    &AAMutSelSparseOmegaModel::NoUpdate, this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CorruptMatrices();

        ProfileMove(nucrelrate, 0.1, 1, 10, &AAMutSelSparseOmegaModel::NucRatesLogProb,
                    &AAMutSelSparseOmegaModel::CorruptMatrices, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &AAMutSelSparseOmegaModel::NucRatesLogProb,
                    &AAMutSelSparseOmegaModel::CorruptMatrices, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &AAMutSelSparseOmegaModel::NucRatesLogProb,
                    &AAMutSelSparseOmegaModel::CorruptMatrices, this);

        ProfileMove(nucstat, 0.1, 1, 10, &AAMutSelSparseOmegaModel::NucRatesLogProb,
                    &AAMutSelSparseOmegaModel::CorruptMatrices, this);
        ProfileMove(nucstat, 0.01, 1, 10, &AAMutSelSparseOmegaModel::NucRatesLogProb,
                    &AAMutSelSparseOmegaModel::CorruptMatrices, this);

        CorruptMatrices();
    }

    //! MH move on omega
    void MoveOmega() {
        omegapathsuffstat.Clear();
        omegapathsuffstat.AddSuffStat(*sitecodonmatrixarray, *sitepathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omega = Random::GammaSample(alpha + omegapathsuffstat.GetCount(),
                                    beta + omegapathsuffstat.GetBeta());
        CorruptCodonMatrices();
    }

    //! MH compensatory move schedule on fitness parameters and hyper-parameters
    void CompMoveFitness() { CompMoveFitness(1.0, 10); }

    //! \brief MH compensatory move on fitness parameters and hyper-parameters
    //!
    //! for a given amino-acid.
    //! shift by a constant factor e all *active* gamma fitness parameters across
    //! all amino-acids and conditions
    double CompMoveFitness(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                vector<double> &x = (*fitness)[i];
                const vector<int> &mask = (*sitemaskarray)[i];
                // vector<int> mask(Naa,1);

                double deltalogprob = 0;

                for (int a = 0; a < Naa; a++) {
                    if (mask[a]) {
                        double alpha = fitnessshape * fitnesscenter[a];
                        deltalogprob -= -Random::logGamma(alpha) + (alpha - 1) * log(x[a]) - x[a];
                    }
                }

                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);

                int n = 0;
                for (int a = 0; a < Naa; a++) {
                    if (mask[a]) {
                        x[a] *= e;
                        n++;
                    }
                }

                double loghastings = n * m;

                for (int a = 0; a < Naa; a++) {
                    if (mask[a]) {
                        double alpha = fitnessshape * fitnesscenter[a];
                        deltalogprob += -Random::logGamma(alpha) + (alpha - 1) * log(x[a]) - x[a];
                    }
                }

                deltalogprob += loghastings;

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    for (int a = 0; a < Naa; a++) {
                        if (mask[a]) {
                            x[a] /= e;
                        }
                    }
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on baseline gamma fitness parameters (for condition k=0)
    void MoveFitness() {
        // if masks are not activated (all entries equal to 1), move a random subset
        // of entries over the 20 amino-acids (2d parameter of call)
        if (maskmode == 3) {
            MoveFitnessAll(1.0, 1, 10);
            MoveFitnessAll(1.0, 3, 10);
            MoveFitnessAll(1.0, 20, 10);
            MoveFitnessAll(0.3, 20, 10);
        }
        // if masks are activated, move all active entries
        else {
            MoveFitness(1.0, 10);
            MoveFitness(0.3, 10);
        }
    }

    //! elementary MH move on baseline gamma fitness parameters (for condition
    //! k=0)
    double MoveFitness(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                vector<double> &x = (*fitness)[i];
                const vector<int> &s = (*sitemaskarray)[i];

                bk = x;

                double deltalogprob = -fitness->GetLogProb(i, s) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, s);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(i, s) + SiteSuffStatLogProb(i);

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

    //! elementary MH move on baseline fitness parameters (for condition k=0):
    //! version used when masks are activated
    double MoveFitnessAll(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                vector<double> &x = (*fitness)[i];

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

    //! MH move schedule on hyperparameter of fitness mask across sites (maskprob,
    //! or pi)
    void MoveMaskHyperParameters() {
        SlidingMove(pi, 1.0, 10, 0.05, 0.975, &AAMutSelSparseOmegaModel::MaskLogProb,
                    &AAMutSelSparseOmegaModel::UpdateMask, this);
        SlidingMove(pi, 0.1, 10, 0.05, 0.975, &AAMutSelSparseOmegaModel::MaskLogProb,
                    &AAMutSelSparseOmegaModel::UpdateMask, this);
    }

    //! MH move schedule on background fitness (maskepsilon)
    void MoveMaskEpsilon() {
        SlidingMove(maskepsilon, 1.0, 10, 0, 1.0, &AAMutSelSparseOmegaModel::MaskEpsilonLogProb,
                    &AAMutSelSparseOmegaModel::UpdateAll, this);
        SlidingMove(maskepsilon, 0.1, 10, 0, 1.0, &AAMutSelSparseOmegaModel::MaskEpsilonLogProb,
                    &AAMutSelSparseOmegaModel::UpdateAll, this);
    }

    //! MH move on fitness masks across sites
    double MoveMasks() {
        double nacc = 0;
        double ntot = 0;
        for (int i = 0; i < Nsite; i++) {
            vector<int> &mask = (*sitemaskarray)[i];
            int naa = 0;
            for (int k = 0; k < Naa; k++) {
                naa += mask[k];
            }
            for (int k = 0; k < Naa; k++) {
                if ((!mask[k]) || (naa > 1)) {
                    double deltalogprob = -MaskLogPrior(i) - SiteSuffStatLogProb(i);
                    naa -= mask[k];
                    mask[k] = 1 - mask[k];
                    naa += mask[k];
                    if (mask[k]) {
                        (*fitness)[i][k] = Random::sGamma(fitnessshape * fitnesscenter[k]);
                        if (!(*fitness)[i][k]) {
                            (*fitness)[i][k] = 1e-8;
                            // cerr << "null fitness : " << fitnessshape << '\t' <<
                            // fitnesscenter[k] << '\n'; exit(1);
                        }
                    }
                    UpdateSite(i);
                    deltalogprob += MaskLogPrior(i) + SiteSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        naa -= mask[k];
                        mask[k] = 1 - mask[k];
                        naa += mask[k];
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on hyperparameters of distribution of fitness factors
    void MoveFitnessHyperParameters() {
        // collect suff stats across all active fitness parameters
        hyperfitnesssuffstat.Clear();
        hyperfitnesssuffstat.AddSuffStat(*fitness, *sitemaskarray);

        if (fitnessshapemode < 2) {
            ScalingMove(fitnessshape, 1.0, 100, &AAMutSelSparseOmegaModel::FitnessHyperLogProb,
                        &AAMutSelSparseOmegaModel::NoUpdate, this);
            ScalingMove(fitnessshape, 0.3, 100, &AAMutSelSparseOmegaModel::FitnessHyperLogProb,
                        &AAMutSelSparseOmegaModel::NoUpdate, this);
            ScalingMove(fitnessshape, 0.1, 100, &AAMutSelSparseOmegaModel::FitnessHyperLogProb,
                        &AAMutSelSparseOmegaModel::NoUpdate, this);
        }
        fitness->SetShape(fitnessshape);

        if (fitnesscentermode < 2) {
            ProfileMove(fitnesscenter, 0.3, 1, 100, &AAMutSelSparseOmegaModel::FitnessHyperLogProb,
                        &AAMutSelSparseOmegaModel::NoUpdate, this);
            ProfileMove(fitnesscenter, 0.1, 1, 100, &AAMutSelSparseOmegaModel::FitnessHyperLogProb,
                        &AAMutSelSparseOmegaModel::NoUpdate, this);
            ProfileMove(fitnesscenter, 0.1, 3, 100, &AAMutSelSparseOmegaModel::FitnessHyperLogProb,
                        &AAMutSelSparseOmegaModel::NoUpdate, this);
        }
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    const CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //-------------------
    // Traces and monitors
    // ------------------

    //! mean fitness over active (high-fitness) amino-acids
    double GetMeanActiveFitness() const {
        double mean = 0;
        int n = 0;
        for (int i = 0; i < GetNsite(); i++) {
            const vector<double> &f = fitness->GetVal(i);
            const vector<int> &m = sitemaskarray->GetVal(i);
            for (int a = 0; a < Naa; a++) {
                mean += m[a] * f[a];
                n += m[a];
            }
        }
        return mean / n;
    }

    //! write header of tracefile
    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "omega\t";
        os << "pi\t";
        os << "width\t";
        os << "meanfitness\t";
        os << "epsilon\t";
        os << "shape\t";
        os << "center\t";
        os << "statent\t";
        os << "rrent\n";
    }

    //! write trace (one line summarizing current state) into trace file
    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << 3 * branchlength->GetTotalLength() << '\t';
        os << omega << '\t';
        os << pi << '\t';
        os << sitemaskarray->GetMeanWidth() << '\t';
        os << GetMeanActiveFitness() << '\t';
        os << maskepsilon << '\t';
        os << fitnessshape << '\t';
        os << Random::GetEntropy(fitnesscenter) << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    //! monitoring MCMC statistics
    void Monitor(ostream &) const override {}

    //! get complete parameter configuration from stream
    void FromStream(istream &is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2) {
            is >> nucrelrate;
            is >> nucstat;
        }
        if (omegamode < 2) {
            is >> omega;
        }
        if (fitnessshapemode < 2) {
            is >> fitnessshape;
        }
        if (fitnesscentermode < 2) {
            is >> fitnesscenter;
        }
        is >> *fitness;
        if (maskmode < 2) {
            is >> pi;
        }
        if (maskmode < 3) {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2) {
            is >> maskepsilon;
        }
    }

    //! write complete current parameter configuration to stream
    void ToStream(ostream &os) const override {
        if (blmode < 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        }
        if (nucmode < 2) {
            os << nucrelrate << '\t';
            os << nucstat << '\t';
        }
        if (omegamode < 2) {
            os << omega << '\t';
        }
        if (fitnessshapemode < 2) {
            os << fitnessshape << '\t';
        }
        if (fitnesscentermode < 2) {
            os << fitnesscenter << '\t';
        }
        os << *fitness << '\t';
        if (maskmode < 2) {
            os << pi << '\t';
        }
        if (maskmode < 3) {
            os << *sitemaskarray << '\t';
        }
        if (maskepsilonmode < 2) {
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
        if (nucmode < 2) {
            size += nucrelrate.size();
            size += nucstat.size();
        }
        if (omegamode < 2) {
            size++;
        }
        if (fitnessshapemode < 2) {
            size++;
        }
        if (fitnesscentermode < 2) {
            size += fitnesscenter.size();
        }
        size += fitness->GetMPISize();
        // pi and epsilon
        if (maskmode < 2) {
            size++;
        }
        if (maskmode < 3) {
            size += sitemaskarray->GetMPISize();
        }
        if (maskepsilonmode < 2) {
            size++;
        }
        return size;
    }

    //! get complete parameter configuration from MPI buffer
    void MPIGet(const MPIBuffer &is) {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2) {
            is >> nucrelrate;
            is >> nucstat;
        }
        if (omegamode < 2) {
            is >> omega;
        }
        if (fitnessshapemode < 2) {
            is >> fitnessshape;
        }
        if (fitnesscentermode < 2) {
            is >> fitnesscenter;
        }
        is >> *fitness;
        if (maskmode < 2) {
            is >> pi;
        }
        if (maskmode < 3) {
            is >> *sitemaskarray;
        }
        if (maskepsilonmode < 2) {
            is >> maskepsilon;
        }
    }

    //! write complete current parameter configuration into MPI buffer
    void MPIPut(MPIBuffer &os) const {
        if (blmode < 2) {
            os << lambda;
            os << *branchlength;
        }
        if (nucmode < 2) {
            os << nucrelrate;
            os << nucstat;
        }
        if (omegamode < 2) {
            os << omega;
        }
        if (fitnessshapemode < 2) {
            os << fitnessshape;
        }
        if (fitnesscentermode < 2) {
            os << fitnesscenter;
        }
        os << *fitness;
        if (maskmode < 2) {
            os << pi;
        }
        if (maskmode < 3) {
            os << *sitemaskarray;
        }
        if (maskepsilonmode < 2) {
            os << maskepsilon;
        }
    }
};
