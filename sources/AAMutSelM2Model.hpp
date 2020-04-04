
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "Chrono.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "MultinomialAllocationVector.hpp"
#include "Permutation.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "StickBreakingProcess.hpp"
#include "Tree.hpp"
#include "MechM2aMix.hpp"
#include "AAMutSelOmegaCodonSubMatrixBidimArray.hpp"

/**
 * \brief The mutation-selection model with constant fitness landscape over the
 * tree -- double Dirichlet process version.
 *
 * The model is parameterized by
 * - an fixed unrooted phylogenetic tree tau, with branch lengths l = (l_j), for
 * j running over branches
 * - a GTR nucleotide matrix Q = rho * pi, specifying the mutation process
 * (assumed homogeneous across sites and lineages)
 * - an array of site-specific amino-acid fitness profiles F_ia, for site i and
 * amino-acid a
 * - an omega multiplier, capturing deviations of the non-syn rate from the
 * model (see Rodrigue and Lartillot, 2107); this parameter is fixed to 1 by
 * default.
 *
 * Site-specific amino-acid fitness profiles are drawn from a Dirichlet process,
 * implemented using a truncated stick-breaking process, of concentration
 * parameter kappa, and truncated at Ncat.
 *
 * The base distribution of this Dirichlet process, G_0, is itself a (truncated
 * stick breaking) mixture of Dirichlet distributions, parameterized by
 * basekappa, and truncated at baseNcat. Each component of this mixture
 * Dirichlet is parameterized by a center (a 20-freqeuncy vector) and a
 * concentration.
 *
 * Concerning the first-level stick-breaking process, by default, Ncat == 100
 * (can be changed using the  -ncat option). As for baseNcat, it is equal to 1,
 * in which case the base distribution G_0 is just a simple Dirichlet (as in
 * Rodrigue et al, 2010). This simple model is probably the one that should be
 * used by default for single-gene analyses. The more complex model with
 * baseNcat > 1 is meant to be used in a multi-gene context (although, even in
 * that case, mixing is still challenging, and not sure whether this model
 * brings important improvement in the end). Thus, baseNcat = 1 appears to be
 * the most reasonable model settings for now.
 *
 * Priors (in a single-gene context):
 * - branch lengths iid exponential, of rate lambda
 * - lambda exponential of rate 10
 * - rho and pi uniform Dirichlet
 * - omega: fixed to 1 or exponential of rate 1
 * - kappa: exponential of rate 0.1
 * - center of base distribution: uniform Dirichlet
 * - concentration of base distribution: exponential of mean 20.
 *
 * In a multi-gene context, shrinkage across genes can be applied to branch
 * lengths, omega, nucleotide rate parameters (rho and pi), and to the
 * parameters of the base distribution (center and concentration) -- see
 * MultiGeneAAMutSelDSBDPModel.
 *
 */

class AAMutSelM2Model : public ProbModel {
    const Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    const CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    double lambda;
    BranchIIDGamma *blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise *branchlength;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStat hyperlengthsuffstat;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucstat;
    std::vector<double> nucrelrate;
    GTRSubMatrix *nucmatrix;

    // prior probability for the gene to be under positive selection (i.e. prior
    // prob that posw > 0)
    double pi;

    // Beta prior for posw (assuming posw>0)
    double poswhypermean;
    double poswhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape
    // parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    double posw;
    double dposom;

    MechM2aMix* componentomegaarray;
    MultinomialAllocationVector *omegasitealloc;
    mutable vector<vector<double>> omegasitepostprobarray;

    OmegaPathSuffStatArray* siteomegapathsuffstatarray;

    // base distribution G0 is itself a stick-breaking mixture of Dirichlet
    // distributions

    int baseNcat;
    int basemin;
    double basekappa;
    StickBreakingProcess *baseweight;
    OccupancySuffStat *baseoccupancy;

    vector<double> basecenterhypercenter;
    double basecenterhyperinvconc;
    IIDDirichlet *basecenterarray;

    double baseconchypermean;
    double baseconchyperinvshape;
    IIDGamma *baseconcentrationarray;

    MultinomialAllocationVector *componentalloc;
    MixtureSelector<vector<double>> *componentcenterarray;
    MixtureSelector<double> *componentconcentrationarray;

    // aa fitness arrays across sites are a SBDP process of base G0 defined above
    int Ncat;
    double kappa;
    StickBreakingProcess *weight;
    OccupancySuffStat *occupancy;

    MultiDirichlet *componentaafitnessarray;
    DirichletSuffStatArray *basesuffstatarray;

    MultinomialAllocationVector *sitealloc;

    // an array of codon matrices (one for each distinct aa fitness profile)
    AAMutSelOmegaCodonSubMatrixBidimArray *componentcodonmatrixbidimarray;

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
    DoubleMixtureSelector<SubMatrix> *sitesubmatrixarray;
    DoubleMixtureSelector<AAMutSelOmegaCodonSubMatrix> *sitecodonmatrixarray;

    PhyloProcess *phyloprocess;

    PathSuffStatArray *sitepathsuffstatarray;
    PathSuffStatBidimArray *componentpathsuffstatbidimarray;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    // currently: shared across genes
    int blmode;
    // currently, free without shrinkage: shared across genes
    int nucmode;
    // currently, shared across genes.
    // free without shrinkage, only with baseNcat = 1
    // free with shrinkage: not really interesting
    int basemode;
    // currently: fixed or free with shrinkage

    Chrono aachrono;
    Chrono basechrono;
    Chrono totchrono;

    double acca1, acca2, acca3, acca4;
    double tota1, tota2, tota3, tota4;
    double accb1, accb2, accb3, accb4;
    double totb1, totb2, totb3, totb4;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor
    //!
    //! parameters:
    //! - datafile: name of file containing codon sequence alignment
    //! - treefile: name of file containing tree topology (and branch conditions,
    //! such as specified by branch names)
    //! - Ncat: truncation of the first-level stick-breaking process (by default:
    //! 100)
    //! - baseNcat: truncation of the second-level stick-breaking process (by
    //! default: 1)

    AAMutSelM2Model(string datafile, string treefile, int inNcat, int inbaseNcat) {
        blmode = 0;
        nucmode = 0;
        basemode = 0;

        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        Ncat = inNcat;
        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100) {
                Ncat = 100;
            }
        }

        basemin = 0;
        if (inbaseNcat < 0) {
            basemin = 1;
            baseNcat = -inbaseNcat;
            if (baseNcat != 2) {
                cerr << "error in basencat\n";
                exit(1);
            }
        } else {
            baseNcat = inbaseNcat;
        }

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        cerr << "ncat : " << Ncat << '\n';
        cerr << "basencat : " << baseNcat << '\n';

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        Tree* tmptree = new Tree(treefile);
        // check whether tree and data fits together
        tmptree->RegisterWith(taxonset);
        tmptree->SetIndices();
        tree = tmptree;

        Nbranch = tree->GetNbranch();

        acca1 = acca2 = acca3 = acca4 = 0;
        tota1 = tota2 = tota3 = tota4 = 0;
        accb1 = accb2 = accb3 = accb4 = 0;
        totb1 = totb2 = totb3 = totb4 = 0;

        // Allocate();
    }

    AAMutSelM2Model(const CodonSequenceAlignment* incodondata, const Tree* intree, int inNcat, int inbaseNcat) {
        blmode = 0;
        nucmode = 0;
        basemode = 0;

        codondata = incodondata;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        Ncat = inNcat;
        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100) {
                Ncat = 100;
            }
        }

        basemin = 0;
        if (inbaseNcat < 0) {
            basemin = 1;
            baseNcat = -inbaseNcat;
            if (baseNcat != 2) {
                cerr << "error in basencat\n";
                exit(1);
            }
        } else {
            baseNcat = inbaseNcat;
        }

        std::cerr << "-- Number of sites: " << Nsite << std::endl;
        cerr << "ncat : " << Ncat << '\n';
        cerr << "basencat : " << baseNcat << '\n';

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();

        acca1 = acca2 = acca3 = acca4 = 0;
        tota1 = tota2 = tota3 = tota4 = 0;
        accb1 = accb2 = accb3 = accb4 = 0;
        totb1 = totb2 = totb3 = totb4 = 0;

        // Allocate();
    }

    //! \brief set estimation method for branch lengths
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //! \brief set estimation method for branch lengths
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetBLMode(int mode) { blmode = mode; }

    //! \brief set estimation method for nuc rates
    //!
    //! - mode == 2: shared and estimated across genes
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! (with shrinkage)
    //! - mode == 0: gene-specific, with fixed hyperparameters (without shrinkage)
    //!
    //! for single-gene analyses, only mode 0 can be used.
    void SetNucMode(int mode) { nucmode = mode; }

    void SetBaseMode(int mode) { basemode = mode; }

    //! allocate the model (data structures)
    void Allocate() {
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

        // base distribution (can be skipped)
        basekappa = 1.0;
        baseweight = new StickBreakingProcess(baseNcat, basekappa);
        baseoccupancy = new OccupancySuffStat(baseNcat);

        basecenterhypercenter.assign(Naa, 1.0 / Naa);
        basecenterhyperinvconc = 1.0 / Naa;

        basecenterarray =
            new IIDDirichlet(baseNcat, basecenterhypercenter, 1.0 / basecenterhyperinvconc);
        basecenterarray->SetUniform();

        baseconchypermean = Naa;
        baseconchyperinvshape = 1.0;
        double alpha = 1.0 / baseconchyperinvshape;
        double beta = alpha / baseconchypermean;

        baseconcentrationarray = new IIDGamma(baseNcat, alpha, beta);
        for (int k = 0; k < baseNcat; k++) {
            (*baseconcentrationarray)[k] = 20.0;
        }
        if (basemin == 1) {
            (*baseconcentrationarray)[0] = 1.0;
        }
        // suff stats for component aa fitness arrays
        basesuffstatarray = new DirichletSuffStatArray(baseNcat, Naa);
        componentalloc = new MultinomialAllocationVector(Ncat, baseweight->GetArray());
        componentcenterarray = new MixtureSelector<vector<double>>(basecenterarray, componentalloc);
        componentconcentrationarray =
            new MixtureSelector<double>(baseconcentrationarray, componentalloc);

        //
        // (truncated) Dirichlet mixture of aa fitness profiles
        //

        // Ncat fitness profiles iid from the base distribution
        componentaafitnessarray =
            new MultiDirichlet(componentcenterarray, componentconcentrationarray);

        // mixture weights (truncated stick breaking process)
        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat, kappa);

        // site allocations to the mixture (multinomial allocation)
        sitealloc = new MultinomialAllocationVector(Nsite, weight->GetArray());

        // occupancy suff stats of site allocations (for resampling weights)
        occupancy = new OccupancySuffStat(Ncat);

        // omega
        dposom = dposomhypermean;
        posw = poswhypermean;
        componentomegaarray = new MechM2aMix(dposom+1, posw);
        omegasitealloc = new MultinomialAllocationVector(GetNsite(), componentomegaarray->GetWeights());
        omegasitepostprobarray.assign(GetNsite(), vector<double>(2, 0));
        siteomegapathsuffstatarray = new OmegaPathSuffStatArray(GetNsite());

        // Ncat mut sel codon matrices (based on the Ncat fitness profiles of the
        // mixture)
        componentcodonmatrixbidimarray = new AAMutSelOmegaCodonSubMatrixBidimArray(
            GetCodonStateSpace(), nucmatrix, componentaafitnessarray, componentomegaarray);

        // selector, specifying which codon matrix should be used for each site
        sitesubmatrixarray = new DoubleMixtureSelector<SubMatrix>(componentcodonmatrixbidimarray, sitealloc, omegasitealloc);
        sitecodonmatrixarray = new DoubleMixtureSelector<AAMutSelOmegaCodonSubMatrix>(componentcodonmatrixbidimarray, sitealloc, omegasitealloc);

        // create polyprocess

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, sitesubmatrixarray);
        phyloprocess->Unfold();

        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatbidimarray = new PathSuffStatBidimArray(Ncat, 2);
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //! \brief const access to array of length-pathsuffstats across branches
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //! \brief const access to array of suff stats for base mixture
    const DirichletSuffStatArray *GetBaseSuffStatArray() const { return basesuffstatarray; }

    //! \brief const access to array of occupancy suff stats for base mixture
    const OccupancySuffStat *GetBaseOccupancies() const { return baseoccupancy; }



    //-------------------
    // Setting and updating
    // ------------------

    //! set branch lengths to a new value (multi-gene analyses)
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

    void SetOmegaMixtureParameters(double inposw, double indposom)   {
        posw = inposw;
        dposom = indposom;
        componentomegaarray->SetParameters(dposom + 1, posw);
    }

    void GetOmegaMixtureParameters(double &inposw, double &indposom) const   {
        inposw = posw;
        indposom = dposom;
    }

    void SetOmegaMixtureHyperParameters(double inpi, double inposwhypermean, double inposwhyperinvconc,
            double indposomhypermean, double indposomhyperinvshape) {

        pi = inpi;
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
        if (!pi) {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }
        if (! dposomhyperinvshape)  {
            dposom = dposomhypermean;
        }
        if (! poswhyperinvconc) {
            posw = poswhypermean;
        }
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
        UpdateMatrices();
    }

    //! copy nucleotide rates into vectors given as arguments (multi-gene
    //! analyses)
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! set base mixture concentration and center parameters to new value
    //! (multi-gene analyses)
    void SetBaseMixture(const Selector<vector<double>> &inbasecenterarray,
                        const Selector<double> &inbaseconcentrationarray,
                        const Selector<double> &inbaseweight, const Selector<int> &inpermut) {
        basecenterarray->Copy(inbasecenterarray);
        baseconcentrationarray->Copy(inbaseconcentrationarray);
        baseweight->Copy(inbaseweight);
        componentalloc->Permute(inpermut);
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatrices() {
        componentcodonmatrixbidimarray->UpdateCodonMatrices();
    }

    //! \brief tell codon matrix k that its parameters have changed and that it
    //! should be updated
    void UpdateCodonMatrix(int k) {
        (*componentcodonmatrixbidimarray)(k,0).CorruptMatrix();
        (*componentcodonmatrixbidimarray)(k,1).CorruptMatrix();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    void Update() override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        componentomegaarray->SetParameters(dposom + 1, posw);
        baseweight->SetKappa(basekappa);
        weight->SetKappa(kappa);
        UpdateBaseOccupancies();
        UpdateOccupancies();
        UpdateMatrices();
        GetIntegratedLogLikelihood();
        ResampleSub(1.0);
    }

    void PostPred(string name) override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        componentomegaarray->SetParameters(dposom + 1, posw);
        baseweight->SetKappa(basekappa);
        weight->SetKappa(kappa);
        UpdateBaseOccupancies();
        UpdateOccupancies();
        UpdateMatrices();
        GetIntegratedLogLikelihood();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior (up to some constant)
    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) {
            total += BranchLengthsLogPrior();
        }
        if (nucmode < 2) {
            total += NucRatesLogPrior();
        }
        if (basemode < 2) {
            if (baseNcat > 1) {
                total += BaseStickBreakingHyperLogPrior();
                total += BaseStickBreakingLogPrior();
            }
            total += BaseLogPrior();
        }
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += AALogPrior();
        total += OmegaLogPrior();
        return total;
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }


    double GetIntegratedLogLikelihood() const {
        int ncat = 2;

        double total = 0;
        double logp[ncat];
        const vector<double> &w = componentomegaarray->GetWeights();
        double max = 0;
        for (int i = 0; i < GetNsite(); i++) {
            for (int k = 0; k < ncat; k++) {
                (*omegasitealloc)[i] = k;
                logp[k] = phyloprocess->SiteLogLikelihood(i);
                if ((!k) || (max < logp[k])) {
                    max = logp[k];
                }
            }

            double p = 0;
            for (int k = 0; k < ncat; k++) {
                double tmp = w[k] * exp(logp[k] - max);
                p += tmp;
                omegasitepostprobarray[i][k] = tmp;
            }
            double logl = log(p) + max;
            total += logl;
            for (int k = 0; k < ncat; k++) {
                omegasitepostprobarray[i][k] /= p;
            }
        }

        omegasitealloc->GibbsResample(omegasitepostprobarray);
        return total;
    }

	void GetSitesPostProb(double *array) const {
	    for (int i = 0; i < GetNsite(); i++) {
            array[i] = omegasitepostprobarray[i][2];
            if (omegasitepostprobarray[i][2] < 0) {
                cerr << "error in AAMutSelM2Model::GetSitesPostProb: negative prob\n";
                exit(1);
            }
	    }
	}

    //! \brief log prior over hyperparameter of prior over branch lengths (here,
    //! lambda ~ exponential of rate 10)
    double LambdaHyperLogPrior() const { return -lambda / 10; }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const {
        double ret = branchlength->GetLogProb();
        if (blmode == 0) {
            ret += LambdaHyperLogPrior();
        }
        return ret;
    }

    //! log prior over dposom and posw
    double OmegaLogPrior() const {
        double total = 0;
        total += PosOmegaLogPrior();
        total += PosWeightLogPrior();
        if (std::isnan(total) || std::isinf(total))  {
            cerr << "error: nan; \n";
            cerr << posw << '\t' << dposom << '\n';
            cerr << dposomhypermean << '\t' << dposomhyperinvshape << '\n';
            cerr << pi << '\t' << poswhypermean << '\t' << poswhyperinvconc << '\n';
            exit(1);
        }
        return total;
    }

    // Gamma prior for dposom
    double PosOmegaLogPrior() const {
        if (! dposomhyperinvshape)  {
            return 0;
        }
        double alpha = 1.0 / dposomhyperinvshape;
        double beta = alpha / dposomhypermean;
        return Random::logGammaDensity(dposom, alpha, beta);
    }

    // mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1
    // - pi) for posw
    double PosWeightLogPrior() const {
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

    // Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogPrior() const {
        if (posw) {
            return log(pi);
        }
        return log(1 - pi);
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

    //! log prior over concentration parameters basekappa of mixture of base
    //! distribution
    double BaseStickBreakingHyperLogPrior() const { return -basekappa / 10; }

    //! log prior over weights of stick breaking process of base distribution
    double BaseStickBreakingLogPrior() const { return baseweight->GetLogProb(basekappa); }

    //! log prior over concentration parameters kappa of stick-breaking mixture of
    //! amino-acid profiles
    double StickBreakingHyperLogPrior() const { return -kappa / 10; }

    //! log prior over weights of stick breaking process of amino-acid profiles
    double StickBreakingLogPrior() const { return weight->GetLogProb(kappa); }

    //! log prior over base center and concentration parameters
    double BaseLogPrior() const {
        double total = 0;
        total += basecenterarray->GetLogProb();
        total += baseconcentrationarray->GetLogProb();
        if (std::isinf(total)) {
            cerr << "in BaseLogPrior: inf\n";
            exit(1);
        }
        return total;
    }

    //! log prior over base center and concentration parameters of component k of
    //! base distribution
    double BaseLogPrior(int k) const {
        double total = 0;
        total += basecenterarray->GetLogProb(k);
        total += baseconcentrationarray->GetLogProb(k);
        return total;
    }

    //! log prior of amino-acid fitness profiles
    double AALogPrior() const { return componentaafitnessarray->GetLogProb(); }

    //! log prior of amino-acid fitness profile k
    double AALogPrior(int k) const { return componentaafitnessarray->GetLogProb(k); }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return componentpathsuffstatbidimarray->GetLogProb(*componentcodonmatrixbidimarray);
    }

    //! return log prob of the substitution mappings over sites allocated to
    //! component k of the mixture
    double PathSuffStatLogProb(int k) const {
        double l0 = componentpathsuffstatbidimarray->GetVal(k,0).GetLogProb(componentcodonmatrixbidimarray->GetVal(k,0));
        double l1 = componentpathsuffstatbidimarray->GetVal(k,1).GetLogProb(componentcodonmatrixbidimarray->GetVal(k,1));
        return l0 + l1;
    }

    //! return log prob of current branch lengths, as a function of branch lengths
    //! hyperparameter lambda
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    //! return log prob of first-level mixture components (i.e. all amino-acid
    //! profiles drawn from component k of the base distribution), as a function
    //! of the center and concentration parameters of this component
    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(basecenterarray->GetVal(k),
                                                       baseconcentrationarray->GetVal(k));
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! log prob factor to be recomputed when moving branch lengths hyperparameter
    //! lambda
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    //! log prob factor to be recomputed when moving nucleotide mutation rate
    //! parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + PathSuffStatLogProb(); }

    //! log prob factor to be recomputed when moving aa hyper params (center and
    //! concentration) for component k of base distribution
    double BaseLogProb(int k) const { return BaseLogPrior(k) + BaseSuffStatLogProb(k); }

    //! log prob factor to be recomputed when moving basekappa, the concentration
    //! parameter of the second-level strick breaking process (base distribution)
    double BaseStickBreakingHyperLogProb() const {
        return BaseStickBreakingHyperLogPrior() + BaseStickBreakingLogPrior();
    }

    //! log prob factor to be recomputed when moving kappa, the concentration
    //! parameter of the first-level strick breaking process
    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

    //-------------------
    //  Collecting Suff Stats
    //-------------------

    //! collect sufficient statistics if substitution mappings across sites
    void CollectSitePathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    //! gather site-specific sufficient statistics component-wise
    void CollectComponentPathSuffStat() {
        componentpathsuffstatbidimarray->Clear();
        componentpathsuffstatbidimarray->Add(*sitepathsuffstatarray, *sitealloc, *omegasitealloc);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    void CollectOmegaPathSuffStat() {
        siteomegapathsuffstatarray->Clear();
        siteomegapathsuffstatarray->AddSuffStat(*sitecodonmatrixarray, *sitepathsuffstatarray);
    }

    double OmegaPathSuffStatLogProb() const {
        componentomegaarray->SetParameters(dposom + 1, posw);
        return componentomegaarray->GetPostProbArray(*siteomegapathsuffstatarray, omegasitepostprobarray);
    }

    double OmegaLogProb() const { return OmegaLogPrior() + OmegaPathSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    //! complete MCMC move schedule
    double Move() override {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        UpdateMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            totchrono.Start();
            if (blmode < 2) {
                MoveBranchLengths();
            }

            CollectSitePathSuffStat();
            CollectComponentPathSuffStat();

            if (nucmode < 2) {
                MoveNucRates();
            }

            MoveOmega();

            aachrono.Start();
            MoveAAMixture(3);
            aachrono.Stop();

            basechrono.Start();
            if (basemode < 2) {
                MoveBase(3);
            }
            basechrono.Stop();

            totchrono.Stop();
        }
    }

    //! MH move on base mixture
    void MoveBase(int nrep) {
        if (baseNcat > 1) {
            ResampleBaseAlloc();
        }
        MoveBaseMixture(nrep);
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
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
        ScalingMove(lambda, 1.0, 10, &AAMutSelM2Model::LambdaHyperLogProb,
                    &AAMutSelM2Model::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &AAMutSelM2Model::LambdaHyperLogProb,
                    &AAMutSelM2Model::NoUpdate, this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    //! MH move on omega
    void MoveOmega() {
        CollectOmegaPathSuffStat();
        if (pi != 0) {
            if (dposomhyperinvshape)    {
                ScalingMove(dposom, 1.0, 10, &AAMutSelM2Model::OmegaLogProb, &AAMutSelM2Model::NoUpdate, this);
            }
            if (poswhyperinvconc)   {
                SlidingMove(posw, 1.0, 10, 0, 1, &AAMutSelM2Model::OmegaLogProb, &AAMutSelM2Model::NoUpdate,
                        this);
            }
        }
        if ((pi != 0) && (pi != 1)) {
            SwitchPosWeight(10);
        }
        ResampleOmegaAlloc();
        UpdateCodonMatrices();
    }

    void ResampleOmegaAlloc() {
        OmegaPathSuffStatLogProb();
        omegasitealloc->GibbsResample(omegasitepostprobarray);
    }

    double SwitchPosWeight(int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double bkposw = posw;
            double deltalogprob = -PosSwitchLogPrior() - OmegaPathSuffStatLogProb();
            if (posw) {
                posw = 0;
            } else {
                double alpha = poswhypermean / poswhyperinvconc;
                double beta = (1 - poswhypermean) / poswhyperinvconc;
                posw = Random::BetaSample(alpha, beta);
            }
            deltalogprob += PosSwitchLogPrior() + OmegaPathSuffStatLogProb();
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

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        ProfileMove(nucrelrate, 0.1, 1, 3, &AAMutSelM2Model::NucRatesLogProb,
                    &AAMutSelM2Model::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &AAMutSelM2Model::NucRatesLogProb,
                    &AAMutSelM2Model::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &AAMutSelM2Model::NucRatesLogProb,
                    &AAMutSelM2Model::UpdateMatrices, this);

        ProfileMove(nucstat, 0.1, 1, 3, &AAMutSelM2Model::NucRatesLogProb,
                    &AAMutSelM2Model::UpdateMatrices, this);
        ProfileMove(nucstat, 0.01, 1, 3, &AAMutSelM2Model::NucRatesLogProb,
                    &AAMutSelM2Model::UpdateMatrices, this);
    }

    //! MCMC module for the mixture amino-acid fitness profiles
    void MoveAAMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveAAProfiles();
            ResampleEmptyComponents();
            ResampleAlloc();
            LabelSwitchingMove();
            ResampleWeights();
            MoveKappa();
            CollectComponentPathSuffStat();
            UpdateCodonMatrices();
        }
    }

    //! resample empty components of the mixture from prior
    void ResampleEmptyComponents() {
        componentaafitnessarray->PriorResample(*occupancy);
        // componentcodonmatrixbidimarray->UpdateCodonMatrices(*occupancy);
    }

    //! MH move on amino-acid fitness profiles (occupied components only)
    void MoveAAProfiles() {
        CompMoveAAProfiles(3);
        MulMoveAAProfiles(3);
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double CompMoveAAProfiles(int nrep) {
        accb1 += MoveAA(1.0, 1, nrep);
        // accb2 += MoveAA(1.0,3,nrep);
        // accb3 += MoveAA(0.3,3,nrep);
        accb4 += MoveAA(0.1, 3, nrep);
        totb1++;
        totb2++;
        totb3++;
        totb4++;
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MulMoveAAProfiles(int nrep) {
        acca1 += MoveAAGamma(3.0, nrep);
        acca2 += MoveAAGamma(1.0, nrep);
        // acca3 += MoveAAGamma(0.3,nrep);
        // acca4 += MoveAAGamma(0.1,nrep);
        tota1++;
        tota2++;
        tota3++;
        tota4++;
        return 1.0;
    }

    //! MH move on amino-acid fitness profiles: additive compensated move on pairs
    //! of entries of the vector
    double MoveAA(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Naa];
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                vector<double> &aa = (*componentaafitnessarray)[i];
                for (int rep = 0; rep < nrep; rep++) {
                    for (int l = 0; l < Naa; l++) {
                        bk[l] = aa[l];
                    }
                    double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                    double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                    deltalogprob += loghastings;
                    UpdateCodonMatrix(i);
                    deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        for (int l = 0; l < Naa; l++) {
                            aa[l] = bk[l];
                        }
                        UpdateCodonMatrix(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! helper function: log density of 20 gammas
    double GammaAALogPrior(const vector<double> &x, const vector<double> &aacenter, double aaconc) {
        double total = 0;
        for (int l = 0; l < Naa; l++) {
            total += (aaconc * aacenter[l] - 1) * log(x[l]) - x[l] -
                     Random::logGamma(aaconc * aacenter[l]);
        }
        return total;
    }

    //! MH move on amino-acid fitness profiles: multiplicative move (using the
    //! Gamma representation of the Dirichlet)
    double MoveAAGamma(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                double aaconc = componentconcentrationarray->GetVal(i);
                const vector<double> &aacenter = componentcenterarray->GetVal(i);

                vector<double> &aa = (*componentaafitnessarray)[i];
                vector<double> x(Naa, 0);
                double z = Random::sGamma(aaconc);
                for (int l = 0; l < Naa; l++) {
                    x[l] = z * aa[l];
                }

                double bkz = z;
                vector<double> bkx = x;
                vector<double> bkaa = aa;

                for (int rep = 0; rep < nrep; rep++) {
                    double deltalogprob =
                        -GammaAALogPrior(x, aacenter, aaconc) - PathSuffStatLogProb(i);

                    double loghastings = 0;
                    z = 0;
                    for (int l = 0; l < Naa; l++) {
                        double m = tuning * (Random::Uniform() - 0.5);
                        double e = exp(m);
                        x[l] *= e;
                        z += x[l];
                        loghastings += m;
                    }
                    for (int l = 0; l < Naa; l++) {
                        aa[l] = x[l] / z;
                        if (aa[l] < 1e-50) {
                            aa[l] = 1e-50;
                        }
                    }

                    deltalogprob += loghastings;

                    UpdateCodonMatrix(i);

                    deltalogprob += GammaAALogPrior(x, aacenter, aaconc) + PathSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                        bkaa = aa;
                        bkx = x;
                        bkz = z;
                    } else {
                        aa = bkaa;
                        x = bkx;
                        z = bkz;
                        UpdateCodonMatrix(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! Gibbs resample mixture allocations
    void ResampleAlloc() {
        vector<double> postprob(Ncat, 0);
        for (int i = 0; i < Nsite; i++) {
            GetAllocPostProb(i, postprob);
            sitealloc->GibbsResample(i, postprob);
        }
        UpdateOccupancies();
    }

    //! update mixture occupancy suff stats (for resampling mixture weights)
    void UpdateOccupancies() {
        occupancy->Clear();
        occupancy->AddSuffStat(*sitealloc);
    }

    //! get allocation posterior probabilities for a given site
    void GetAllocPostProb(int site, vector<double> &postprob) {
        double max = 0;
        const vector<double> &w = weight->GetArray();
        const PathSuffStat &suffstat = sitepathsuffstatarray->GetVal(site);
        int omalloc = omegasitealloc->GetVal(site);
        for (int i = 0; i < Ncat; i++) {
            double tmp = suffstat.GetLogProb(componentcodonmatrixbidimarray->GetVal(i, omalloc));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) {
                max = tmp;
            }
        }

        double total = 0;
        for (int i = 0; i < Ncat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < Ncat; i++) {
            postprob[i] /= total;
        }
    }

    //! MCMC sequence for label switching moves
    void LabelSwitchingMove() {
        Permutation permut(Ncat);
        weight->LabelSwitchingMove(5, *occupancy, permut);
        sitealloc->Permute(permut);
        componentaafitnessarray->Permute(permut);
    }

    //! Gibbs resample mixture weights (based on occupancy suff stats)
    void ResampleWeights() { weight->GibbsResample(*occupancy); }

    //! MH move on kappa, concentration parameter of the mixture
    void MoveKappa() {
        ScalingMove(kappa, 1.0, 10, &AAMutSelM2Model::StickBreakingHyperLogProb,
                    &AAMutSelM2Model::NoUpdate, this);
        ScalingMove(kappa, 0.3, 10, &AAMutSelM2Model::StickBreakingHyperLogProb,
                    &AAMutSelM2Model::NoUpdate, this);
        weight->SetKappa(kappa);
    }

    //! MCMC module for the base mixture
    void MoveBaseMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseComponents(10);
            ResampleBaseEmptyComponents();
            if (baseNcat > 1) {
                if (!basemin) {
                    BaseLabelSwitchingMove();
                }
                ResampleBaseWeights();
                MoveBaseKappa();
            }
        }
    }

    //! MCMC module for moving the center and concentration parameters of the
    //! components of the the base mixture
    void MoveBaseComponents(int nrep) {
        CollectBaseSuffStat();
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseCenters(1.0, 1);
            MoveBaseCenters(1.0, 3);
            MoveBaseCenters(0.3, 3);
            MoveBaseConcentrations(1.0);
            MoveBaseConcentrations(0.3);
        }
    }

    //! collect suff stats for moving center and concentration parameters of the
    //! base mixture
    void CollectBaseSuffStat() {
        basesuffstatarray->Clear();
        componentaafitnessarray->AddSuffStat(*basesuffstatarray, *componentalloc);
    }

    //! MCMC module for moving the center parameters of the components of the the
    //! base mixture
    double MoveBaseCenters(double tuning, int n) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);
        for (int k = 0; k < baseNcat; k++) {
            if (baseoccupancy->GetVal(k)) {
                vector<double> &aa = (*basecenterarray)[k];
                bk = aa;
                double deltalogprob = -BaseLogProb(k);
                double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                deltalogprob += loghastings;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    aa = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MCMC module for moving the concentration parameters of the components of
    //! the the base mixture
    double MoveBaseConcentrations(double tuning) {
        double nacc = 0;
        double ntot = 0;
        for (int k = basemin; k < baseNcat; k++) {
            if (baseoccupancy->GetVal(k)) {
                double &c = (*baseconcentrationarray)[k];
                double bk = c;
                double deltalogprob = -BaseLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                c *= e;
                deltalogprob += m;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    c = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! resample empty components of the base mixture from the prior
    void ResampleBaseEmptyComponents() {
        basecenterarray->PriorResample(*baseoccupancy);
        baseconcentrationarray->PriorResample(*baseoccupancy);
        if (basemin == 1) {
            (*baseconcentrationarray)[0] = 1.0;
        }
    }

    //! Gibbs resample base mixture allocations
    void ResampleBaseAlloc() {
        vector<double> postprob(baseNcat, 0);
        for (int i = 0; i < Ncat; i++) {
            GetBaseAllocPostProb(i, postprob);
            componentalloc->GibbsResample(i, postprob);
            if ((componentalloc->GetVal(i) < 0) || (componentalloc->GetVal(i) >= baseNcat)) {
                cerr << "error in ResampleBaseAlloc: out of bound\n";
                exit(1);
            }
        }
        UpdateBaseOccupancies();
    }

    //! update base occupancy suff stats (for moving base weights)
    void UpdateBaseOccupancies() {
        baseoccupancy->Clear();
        baseoccupancy->AddSuffStat(*componentalloc);
    }

    //! get allocation posterior probability of a component of the first-level
    //! mixture to the components of the second-level mixture
    void GetBaseAllocPostProb(int cat, vector<double> &postprob) {
        double max = 0;
        const vector<double> &w = baseweight->GetArray();
        for (int i = 0; i < baseNcat; i++) {
            double tmp = Random::logDirichletDensity(componentaafitnessarray->GetVal(cat),
                                                     basecenterarray->GetVal(i),
                                                     baseconcentrationarray->GetVal(i));
            postprob[i] = tmp;
            if ((!i) || (max < tmp)) {
                max = tmp;
            }
        }

        double total = 0;
        for (int i = 0; i < baseNcat; i++) {
            postprob[i] = w[i] * exp(postprob[i] - max);
            total += postprob[i];
        }

        for (int i = 0; i < baseNcat; i++) {
            postprob[i] /= total;
        }
    }

    //! MCMC sequence for label switching moves of the base mixture
    void BaseLabelSwitchingMove() {
        Permutation permut(baseNcat);
        baseweight->LabelSwitchingMove(5, *baseoccupancy, permut);
        componentalloc->Permute(permut);
        basecenterarray->Permute(permut);
        baseconcentrationarray->Permute(permut);
        basesuffstatarray->Permute(permut);
    }

    //! Gibbs resample base mixture weights (based on occupancy suff stats)
    void ResampleBaseWeights() { baseweight->GibbsResample(*baseoccupancy); }

    //! MH move on basekappa, concentration parameter of the base mixture
    void MoveBaseKappa() {
        ScalingMove(basekappa, 1.0, 10, &AAMutSelM2Model::BaseStickBreakingHyperLogProb,
                    &AAMutSelM2Model::NoUpdate, this);
        ScalingMove(basekappa, 0.3, 10, &AAMutSelM2Model::BaseStickBreakingHyperLogProb,
                    &AAMutSelM2Model::NoUpdate, this);
        baseweight->SetKappa(basekappa);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    //! return number of occupied components in first-level mixture (mixture of
    //! amino-acid fitness profiles)
    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                n++;
            }
        }
        return n;
    }

    //! return number of occupied components in base distribution
    int GetBaseNcluster() const {
        int n = 0;
        for (int i = 0; i < baseNcat; i++) {
            if (baseoccupancy->GetVal(i)) {
                n++;
            }
        }
        return n;
    }

    //! return mean entropy of amino-acd fitness profiles
    double GetMeanAAEntropy() const { return componentaafitnessarray->GetMeanEntropy(); }

    //! return mean of concentration parameters of base distribution
    double GetMeanComponentAAConcentration() const {
        double tot = 0;
        double totw = 0;
        for (int i = basemin; i < baseNcat; i++) {
            tot += baseoccupancy->GetVal(i) * baseconcentrationarray->GetVal(i);
            totw += baseoccupancy->GetVal(i);
        }
        return tot / totw;
    }

    //! return mean entropy of centers of base distribution
    double GetMeanComponentAAEntropy() const {
        double tot = 0;
        for (int i = 0; i < baseNcat; i++) {
            tot += baseoccupancy->GetVal(i) * Random::GetEntropy(basecenterarray->GetVal(i));
        }
        return tot / Ncat;
    }

    //! return entropy of vector of nucleotide exchange rates
    double GetNucRREntropy() const { return Random::GetEntropy(nucrelrate); }

    //! return entropy of vector of equilibrium nucleotide composition
    double GetNucStatEntropy() const { return Random::GetEntropy(nucrelrate); }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "posw\t";
        os << "posom\t";
        os << "ncluster\t";
        os << "kappa\t";
        if (baseNcat > 1) {
            if (basemin) {
                os << "basencluster\t";
                os << "baseweight1\t";
            } else {
                os << "basencluster\t";
                os << "basekappa\t";
            }
        }
        os << "aaent\t";
        os << "meanaaconc\t";
        os << "aacenterent\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        // 3x: per coding site (and not per nucleotide site)
        os << 3 * branchlength->GetTotalLength() << '\t';
        os << posw << '\t';
        os << 1.0 + dposom << '\t';
        os << GetNcluster() << '\t';
        os << kappa << '\t';
        if (baseNcat > 1) {
            if (basemin) {
                os << GetBaseNcluster() << '\t';
                os << baseweight->GetVal(1) << '\t';
            } else {
                os << GetBaseNcluster() << '\t';
                os << basekappa << '\t';
            }
        }
        os << GetMeanAAEntropy() << '\t';
        os << GetMeanComponentAAConcentration() << '\t';
        os << GetMeanComponentAAEntropy() << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream &os) const override {
        os << totchrono.GetTime() << '\t' << aachrono.GetTime() << '\t' << basechrono.GetTime()
           << '\n';
        os << "prop time in aa moves  : " << aachrono.GetTime() / totchrono.GetTime() << '\n';
        os << "prop time in base moves: " << basechrono.GetTime() / totchrono.GetTime() << '\n';
    }

    void FromStream(istream &is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2) {
            is >> nucrelrate;
            is >> nucstat;
        }
        is >> posw >> dposom;
        if (basemode < 2) {
            is >> basekappa;
            baseweight->FromStreamSB(is);
            // is >> *baseweight;
            is >> *componentalloc;
            is >> *basecenterarray;
            is >> *baseconcentrationarray;
        }
        is >> kappa;
        weight->FromStreamSB(is);
        // is >> *weight;
        is >> *componentaafitnessarray;
        is >> *sitealloc;
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
        os << posw << '\t' << dposom << '\t';
        if (basemode < 2) {
            os << basekappa << '\t';
            baseweight->ToStreamSB(os);
            // os << *baseweight << '\t';
            os << *componentalloc << '\t';
            os << *basecenterarray << '\t';
            os << *baseconcentrationarray << '\t';
        }
        os << kappa << '\t';
        weight->ToStreamSB(os);
        // os << *weight << '\t';
        os << *componentaafitnessarray << '\t';
        os << *sitealloc << '\t';
    }

    //! return size of model, when put into an MPI buffer (in multigene context --
    //! only omegatree)
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
        size += 2;
        if (basemode < 2) {
            size++;
            size += baseweight->GetMPISizeSB();
            size += componentalloc->GetMPISize();
            size += basecenterarray->GetMPISize();
            size += baseconcentrationarray->GetMPISize();
        }
        size++;
        size += weight->GetMPISizeSB();
        size += componentaafitnessarray->GetMPISize();
        size += sitealloc->GetMPISize();
        return size;
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer &is) {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2) {
            is >> nucrelrate;
            is >> nucstat;
        }
        is >> posw >> dposom;
        if (basemode < 2) {
            is >> basekappa;
            baseweight->MPIGetSB(is);
            is >> *componentalloc;
            is >> *basecenterarray;
            is >> *baseconcentrationarray;
        }
        is >> kappa;
        weight->MPIGetSB(is);
        is >> *componentaafitnessarray;
        is >> *sitealloc;
    }

    //! write array into MPI buffer
    void MPIPut(MPIBuffer &os) const {
        if (blmode < 2) {
            os << lambda;
            os << *branchlength;
        }
        if (nucmode < 2) {
            os << nucrelrate;
            os << nucstat;
        }
        os << posw << dposom;
        if (basemode < 2) {
            os << basekappa;
            baseweight->MPIPutSB(os);
            os << *componentalloc;
            os << *basecenterarray;
            os << *baseconcentrationarray;
        }
        os << kappa;
        weight->MPIPutSB(os);
        os << *componentaafitnessarray;
        os << *sitealloc;
    }
};
