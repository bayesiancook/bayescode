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

/**
 * \brief The differential selection (diffsel) model, such as implemented in
 * Parto and Lartillot, 2017, BMC Evolutionary Biology.
 *
 * This is a codon model, based on the mutation-selection formalism,
 * such that the fitness landscape (defined by site-specific amimo-acid fitness
 * profiles) is modulated across a set of alternative 'conditions', or
 * environments, specified over the tree. The model is meant to identify
 * positions displaying significant patterns of directional convergent selection
 * (i.e. sites that tend to substitute their amino-acid state in a consistent
 * manner, upon repeated transitions into a specific ecological condition).
 *
 * Technically, the model defines K conditions;
 * condition 0 defines the background, and then condition k=1..K-1 represent the
 * alternative conditions. Branches are a priori allocated to any one of the K
 * conditions. The parameterization of the model is as follows: An array of
 * site-specific basal fitness profiles G_ia, for site i and amino-acid a, is
 * defined for the background condition. Then, for each alternative condition, a
 * differential (log-fitness) contrast is defined: D_kia, for condition k, site
 * i and amino-acid a. A positive (resp negative) D_kia means that the fitness
 * of amino-acid a at site i under condition k is increased (resp decreased)
 * compared to the baseline. Quantitatively, the fitness vector at site i under
 * condition k is defined by F_kia = G_kia * exp(D_kia).
 *
 * The constrasts D_kia are iid normal of mean 0 and variance that is either
 * fixed (equal to 1) or estimated across all amino-acids and sites (separately
 * for each condition) -- see BidimIIDMVNormal.
 *
 * Statistical support for a differential effect between conditions, for a given
 * site i and a given amino acid a, is quantified by the tail area of the
 * posterior distribution of D_kia (i.e. to what extent this posterior
 * distribution is mostly under either the right-hand or the left-hand site of
 * 0).
 */

#include "AADiffSelCodonMatrixBidimArray.hpp"
#include "BranchAllocationSystem.hpp"
#include "CodonSequenceAlignment.hpp"
#include "DiffSelFitnessArray.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "IIDMVNormal.hpp"
#include "PathSuffStat.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "SubMatrixSelector.hpp"
#include "Tree.hpp"

class DiffSelModel : public ProbModel {
    // -----
    // model selectors
    // -----

    int codonmodel;
    int blmode;
    int nucmode;
    int fixglob;
    int fixvar;

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

    // number of diff sel categories
    int Ncond;

    // number of levels of the model
    // with 2 levels, structure of the model is as follows:
    // baseline (condition 0)
    // baseline  * exp(*delta1) (condition 1)
    // baseline * exp(delta1 + deltak) (for condition k=2..Ncond)
    int Nlevel;

    // which branch is under which condition
    BranchAllocationSystem *branchalloc;

    // auxiliary array
    // condalloc[k][l] is 1 iff condition l should recompute its fitness profiles
    // whenever condition k has changed has changed
    vector<vector<int>> condalloc;

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

    // baseline (global) fitness profiles across sites
    IIDDirichlet *baseline;

    // variance parameters (across conditions k=1..Ncond)
    IIDGamma *varsel;

    // differential selection factors across conditions k=1..Ncond and across
    // sites
    BidimIIDMVNormal *delta;

    // fitness profiles (combinations of baseline and delta)
    // across conditions and across sites
    DiffSelFitnessArray *fitnessprofile;

    // codon substitution matrices
    // across conditions and sites
    AADiffSelCodonMatrixBidimArray *condsubmatrixarray;

    // branch- and site-substitution matrices (for phyloprocess)
    SubMatrixSelector *submatrixarray;
    // and for root (condition 0)
    RootSubMatrixSelector *rootsubmatrixarray;

    // phyloprocess
    PhyloProcess *phyloprocess;

    // suff stats

    // path suff stats across conditions and sites
    PathSuffStatBidimArray *suffstatarray;

  public:
    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions,
    //! such as specified by branch names)
    //! - Ncond: number of conditions (K)
    //! - Nlevel: number of levels (if Nlevel == 1: each condition is defined
    //! w.r.t. condition 0; if Nlevel == 2, condition 1 is defined w.r.t.
    //! condition 0, and condition 2..K-1 all defined w.r.t. condition 1)
    //! - fixglob: whether the hyperparameters for the global selection profiles
    //! (G_i) are fixed (i.e. fitness profiles iid uniform Dirichlet) or estimated
    //! across sites
    //! - fixvar: whether the variance of the differential effects D_ki is fixed
    //! (equal to 1) or estimated across sites (separately for each condition)
    //! - codonmodel: type of codon substitution model (1: canonical
    //! mutation-selection model, 0: square-root model, see Parto and Lartillot,
    //! 2017)
    DiffSelModel(const std::string &datafile, const std::string &treefile, int inNcond,
                 int inNlevel, int infixglob, int infixvar, int incodonmodel) {
        blmode = 0;
        nucmode = 0;

        fixglob = infixglob;
        if (!fixglob) {
            cerr << "error: free hyperparameters for baseline (global profile) not "
                    "yet "
                    "implemented\n";
            exit(1);
        }
        fixvar = infixvar;
        codonmodel = incodonmodel;
        Ncond = inNcond;

        Nlevel = inNlevel;
        if (Ncond <= 2) {
            Nlevel = 1;
        }

        ReadFiles(datafile, treefile);
        MakePattern();

        // specifies which condition for which branch
        branchalloc = new BranchAllocationSystem(*tree, Ncond);
        std::cerr << "-- conditions over branches ok\n";
    }

    DiffSelModel(const DiffSelModel &) = delete;

    ~DiffSelModel() {}

    //! read files (and read out the distribution of conditions across branches,
    //! based on the tree read from treefile)
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

        std::cerr << "-- Number of taxa : " << Ntaxa << '\n';
        std::cerr << "-- Number of branches : " << Nbranch << '\n';

        std::cerr << "-- Tree and data fit together\n";
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

        // baseline (global) profile
        // uniform Dirichlet distributed
        vector<double> center(20, 1.0 / 20);
        double concentration = 20.0;
        baseline = new IIDDirichlet(Nsite, center, concentration);
        // baseline->SetUniform();

        // variance parameters (one for each condition, 1..Ncond)
        double shape = 1.0;
        double scale = 1.0;
        varsel = new IIDGamma(Ncond - 1, shape, scale);
        for (int k = 0; k < Ncond - 1; k++) {
            (*varsel)[k] = 1.0;
        }

        // differential selection effects
        // normally distributed
        delta = new BidimIIDMVNormal(Nsite, 20, *varsel);
        // delta->SetToZero();

        // fitnessprofiles...
        fitnessprofile = new DiffSelFitnessArray(*baseline, *delta, Nlevel);
        // fitnessprofile = new DiffSelFitnessArray(*baseline,*delta,condalloc);

        // codon matrices
        // per condition and per site
        condsubmatrixarray =
            new AADiffSelCodonMatrixBidimArray(*fitnessprofile, *GetCodonStateSpace(), *nucmatrix);

        // sub matrices per branch and per site
        submatrixarray = new SubMatrixSelector(*condsubmatrixarray, *branchalloc);
        // sub matrices for root, across sites
        rootsubmatrixarray = new RootSubMatrixSelector(*condsubmatrixarray);

        // create phyloprocess
        phyloprocess =
            new PhyloProcess(tree, codondata, branchlength, 0, submatrixarray, rootsubmatrixarray);
        phyloprocess->Unfold();

        // create suffstat arrays
        suffstatarray = new PathSuffStatBidimArray(Ncond, Nsite);
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
        condsubmatrixarray->Corrupt();
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

    void Update() override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        UpdateAll();
        ResampleSub(1.0);
    }

    //! update fitness profiles and matrices across all sites and conditions
    void UpdateAll() {
        fitnessprofile->Update();
        CorruptMatrices();
    }

    //! update fitness profiles and matrices across all conditions for site i
    void UpdateSite(int i) {
        fitnessprofile->UpdateColumn(i);
        condsubmatrixarray->CorruptColumn(i);
        CorruptMatrices();
    }

    //! update fitness profile and matrix for site i and condition k
    void UpdateSiteCond(int i, int k) {
        fitnessprofile->UpdateColumn(i, condalloc[k + 1]);
        condsubmatrixarray->CorruptColumn(i, condalloc[k + 1]);
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
        if (nucmode < 2) {
            total += NucRatesLogPrior();
        }

        // branchlengths
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();

        // nuc rates
        total += NucRatesLogPrior();

        // uniform on baseline
        total += BaselineLogPrior();

        // variance parameters
        total += VarSelLogPrior();

        // differential selection effects
        total += DeltaLogPrior();

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

    //! log prior over baseline (G)
    double BaselineLogPrior() const { return Nsite * Random::logGamma((double)Naa); }

    //! log prior over differential effect (D)
    double DeltaLogPrior() const { return delta->GetLogProb(); }

    //! log hyper prior over variance parameters of D
    double VarSelLogPrior() const { return varsel->GetLogProb(); }

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

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        suffstatarray->Clear();
        suffstatarray->AddSuffStat(*phyloprocess, *branchalloc);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    double SuffStatLogProb() const { return suffstatarray->GetLogProb(*condsubmatrixarray); }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process, at site i
    double SiteSuffStatLogProb(int site) const {
        return suffstatarray->GetLogProb(site, *condsubmatrixarray);
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process, at site i and for branches
    //! under condition k
    double SiteCondSuffStatLogProb(int site, int k) {
        return suffstatarray->GetLogProb(site, condalloc[k + 1], *condsubmatrixarray);
    }

    //! \brief return log prob of current branch lengths, as a function of branch
    //! lengths hyperparameter lambda
    double BranchLengthsHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
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

    // ---------------
    // Moves
    // ---------------

    /*
    // old move schedule
    double Move() override {
        int nrep0 = 3;
        int nrep = 20;

        for (int rep0 = 0; rep0 < nrep0; rep0++) {
            CollectLengthSuffStat();

            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();

            CollectPathSuffStat();

            UpdateAll();

            for (int rep = 0; rep < nrep; rep++) {
                MoveBaseline();
                UpdateAll();
                MoveDelta();
                if (!fixvar) {
                    MoveVarSel();
                }
            }
            MoveNucRates();
        }

        UpdateAll();

        ResampleSub(1.0);
        return 1.0;
    }
    */

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

            CollectPathSuffStat();
            UpdateAll();

            for (int rep = 0; rep < nrep; rep++) {
                MoveBaseline();
                UpdateAll();
                MoveDelta();
                if (!fixvar) {
                    MoveVarSel();
                }
            }

            if (nucmode < 2) {
                MoveNucRates();
            }
        }
        UpdateAll();
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) { phyloprocess->Move(frac); }

    //! MCMC move schedule on branch lengths
    void MoveBranchLengths() {
        ResampleBranchLengths();
        if (blmode == 0) {
            MoveLambda();
        }
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &DiffSelModel::BranchLengthsHyperLogProb,
                    &DiffSelModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &DiffSelModel::BranchLengthsHyperLogProb,
                    &DiffSelModel::NoUpdate, this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CorruptMatrices();

        ProfileMove(nucrelrate, 0.1, 1, 10, &DiffSelModel::NucRatesLogProb,
                    &DiffSelModel::CorruptMatrices, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &DiffSelModel::NucRatesLogProb,
                    &DiffSelModel::CorruptMatrices, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &DiffSelModel::NucRatesLogProb,
                    &DiffSelModel::CorruptMatrices, this);

        ProfileMove(nucstat, 0.1, 1, 10, &DiffSelModel::NucRatesLogProb,
                    &DiffSelModel::CorruptMatrices, this);
        ProfileMove(nucstat, 0.01, 1, 10, &DiffSelModel::NucRatesLogProb,
                    &DiffSelModel::CorruptMatrices, this);

        CorruptMatrices();
    }

    //! MH moves on baseline fitness profiles (G)
    void MoveBaseline() { MoveBaseline(0.15, 10, 1); }

    //! MH moves on baseline fitness profiles (G)
    double MoveBaseline(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                bk = (*baseline)[i];

                double deltalogprob = -baseline->GetLogProb(i) - SiteSuffStatLogProb(i);
                double loghastings = Random::ProfileProposeMove((*baseline)[i], Naa, tuning, n);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += baseline->GetLogProb(i) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    (*baseline)[i] = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH moves on differential fitness effects (D)
    void MoveDelta() {
        for (int k = 0; k < Ncond - 1; k++) {
            // for (int k = 1; k < Ncond; k++) {
            MoveDelta(k, 5, 1, 10);
            MoveDelta(k, 3, 5, 10);
            MoveDelta(k, 1, 10, 10);
            MoveDelta(k, 1, 20, 10);
        }
    }

    //! MH moves on differential fitness effects (D)
    double MoveDelta(int k, double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);
        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                bk = delta->GetVal(k, i);
                double deltalogprob = -delta->GetLogProb(k, i) - SiteCondSuffStatLogProb(i, k);
                double loghastings = Random::RealVectorProposeMove((*delta)(k, i), Naa, tuning, n);
                deltalogprob += loghastings;
                UpdateSiteCond(i, k);
                deltalogprob += delta->GetLogProb(k, i) + SiteCondSuffStatLogProb(i, k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    (*delta)(k, i) = bk;
                    UpdateSiteCond(i, k);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH moves on variance parameters of differential fitness effects
    void MoveVarSel() {
        MoveVarSel(1.0, 10);
        MoveVarSel(0.3, 10);
    }

    //! MH moves on variance parameters of differential fitness effects
    double MoveVarSel(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int k = 1; k < Ncond; k++) {
                double deltalogprob = -varsel->GetLogProb(k) - delta->GetRowLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*varsel)[k] *= e;
                deltalogprob += varsel->GetLogProb(k) + delta->GetRowLogProb(k);
                deltalogprob += m;
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    (*varsel)[k] /= e;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    const CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return number of sites
    int GetNsite() const { return Nsite; }
    //! return number of conditions
    int GetNcond() const { return Ncond; }

    const vector<vector<double>> &GetCondDeltaArray(int k) const {
        return delta->GetSubArray(k - 1);
    }

    //! get a copy of baseline array
    void GetBaselineArray(double *array) const {
        int j = 0;
        for (int i = 0; i < GetNsite(); i++) {
            for (int a = 0; a < Naa; a++) {
                array[j++] = baseline->GetVal(i)[a];
            }
        }
    }

    //! get a copy of delta array for condition k
    void GetDeltaArray(int k, double *array) const {
        int j = 0;
        for (int i = 0; i < GetNsite(); i++) {
            for (int a = 0; a < Naa; a++) {
                array[j++] = delta->GetVal(k - 1, i)[a];
            }
        }
    }

    //-------------------
    // Traces and monitors
    // ------------------

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "globent\t";
        for (int k = 1; k < Ncond; k++) {
            os << "var" << k << '\t';
        }
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << baseline->GetMeanEntropy() << '\t';
        for (int k = 0; k < Ncond - 1; k++) {
            os << delta->GetMeanVar(k) << '\t';
        }
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    //! trace the current value of baselines, across all sites and all amino-acids
    void TraceBaseline(ostream &os) const {
        for (int i = 0; i < GetNsite(); i++) {
            for (int a = 0; a < Naa; a++) {
                os << baseline->GetVal(i)[a] << '\t';
            }
        }
        os << '\n';
    }

    //! trace the current value of deltas, across all sites and all amino-acids,
    //! under condition k (one single line in output stream)
    void TraceDelta(int k, ostream &os) const {
        for (int i = 0; i < GetNsite(); i++) {
            for (int a = 0; a < Naa; a++) {
                os << delta->GetVal(k - 1, i)[a] << '\t';
            }
        }
        os << '\n';
    }

    void Monitor(ostream &) const override {}

    void FromStream(istream &is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2) {
            is >> nucrelrate;
            is >> nucstat;
        }
        is >> *baseline;
        is >> *varsel;
        is >> *delta;
    }

    void ToStream(ostream &os) const override {
        if (blmode < 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        }
        if (nucmode < 2) {
            os << nucrelrate << '\t';
            os << nucstat << '\t';
        }
        os << *baseline << '\n';
        os << *varsel << '\n';
        os << *delta << '\n';
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
        size += varsel->GetMPISize();
        size += baseline->GetMPISize();
        size += delta->GetMPISize();
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
        is >> *varsel;
        is >> *baseline;
        is >> *delta;
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
        os << *varsel;
        os << *baseline;
        os << *delta;
    }

  private:
    void MakePattern() {
        condalloc.assign(Ncond, vector<int>(Ncond, 0));
        for (int l = 0; l < Ncond; l++) {
            condalloc[0][l] = 1;
        }
        if (Nlevel == 2) {
            for (int l = 1; l < Ncond; l++) {
                condalloc[1][l] = 1;
            }
        }
        for (int l = Nlevel; l < Ncond; l++) {
            condalloc[l][l] = 1;
        }
    }
};
