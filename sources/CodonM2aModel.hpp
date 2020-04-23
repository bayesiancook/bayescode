
/**
 * \brief The M2amodel of codeml (Muse and Gaut version)
 *
 * the omega_i's across sites (i=1..Nsite) are a mixture with 3 components
 * - omega0 < 1, with weight w0
 * - omega1 = 1, with weight w1
 * - omega2 > 1, with weight w2
 *
 * This model is used to test for the presence of positive selection (i.e to
 * test whether w2>0) in a gene and then to select those sites that have a dN/dS
 * > 1 (i.e. that are allocated to the third category of the mixture with high
 * post prob). Here, the model is parameterized as follows:
 * - omega0 = purom,
 * - omega1 = 1,
 * - omega2 = 1 + dposom,
 * where 0 < purom < 1 and dposom > 0;
 * purom has a beta prior (hyperparams: puromhypermean and puromhyperinvconc);
 * dposom has a gamma prior (hyperparams: dposomhypermean and
 * dposomhyperinvshape).
 *
 * The weights of the mixture are parameterized as follows:
 * - w0 = purw * (1 - posw)
 * - w1 = (1-purw) * (1-posw)
 * - w2 = posw
 * where 0<purw<1 and 0<=posw<1;
 * purw has a beta prior (hyperparams: purwhypermean and purwhyperinvconc);
 * the prior on posw is a mixture:
 * - with probability 1-pi, posw = 0
 * - with probability pi, 0 < posw < 1, in which case is it from a beta prior
 * (hyperparams: poswhypermean and poswhyperinvconc). Thus, setting pi = 0
 * imposes a model without positive selection.
 *
 * In total, the 9 hyperparameters of the mixture of omegas are:
 * puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape,
 * purwhypermean, purwhyperinvconc, pi, poswhypermean, poswhyperinvconc.
 * In a single-gene context, these hyperparameters are fixed;
 * in a multigene context, they can be either fixed or estimated across genes
 * (see MultiGeneCodonM2aModel).
 */

#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "M2aMix.hpp"
#include "MultinomialAllocationVector.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "WhiteNoise.hpp"

class CodonM2aModel : public ProbModel {
  public:
    //-------------------
    // Constructors
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files, and pi,
    //! proportion of sites under positive selection
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    CodonM2aModel(const CodonSequenceAlignment* incodondata, const Tree* intree, double inpi)   {
        blmode = 0;
        nucmode = 0;
        data = 0;
        codondata = incodondata;
        pi = inpi;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        tree = intree;
        Nbranch = tree->GetNbranch();
    }

    CodonM2aModel(string datapath, string datafile, string treefile, double inpi)   {
        blmode = 0;
        nucmode = 0;
        data = new FileSequenceAlignment(datapath + datafile);
        if (data->GetNsite() % 3) {
            cerr << "error : not a correctly formatted codon-alignment: " << datafile << '\n';
            cerr << "path  : " << datapath << '\n';
            exit(1);
        }
        codondata = new CodonSequenceAlignment(data, true);
        pi = inpi;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        Tree* tmptree = new Tree(datapath + treefile);
        // check whether tree and data fits together
        tmptree->RegisterWith(taxonset);
        tmptree->SetIndices();
        tree = tmptree;

        Nbranch = tree->GetNbranch();
    }


    //! model allocation
    void Allocate() {

        lambda = 10.0;
        blhypermean = new BranchIIDGamma(*tree, 1.0, lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree, *blhypermean, 1.0 / blhyperinvshape);

        purom = puromhypermean;
        dposom = dposomhypermean;
        purw = purwhypermean;
        posw = poswhypermean;

        componentomegaarray = new M2aMix(purom, dposom + 1, purw, posw);
        sitealloc = new MultinomialAllocationVector(GetNsite(), componentomegaarray->GetWeights());
        sitepostprobarray.assign(GetNsite(), vector<double>(3, 0));

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);

        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        componentcodonmatrixarray = new MGOmegaCodonSubMatrixArray(
            (CodonStateSpace *)codondata->GetStateSpace(), nucmatrix, componentomegaarray);

        sitesubmatrixarray = new MixtureSelector<SubMatrix>(componentcodonmatrixarray, sitealloc);
        sitecodonmatrixarray =
            new MixtureSelector<MGOmegaCodonSubMatrix>(componentcodonmatrixarray, sitealloc);

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, sitesubmatrixarray);
        phyloprocess->Unfold();

        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        sitepathsuffstatarray = new PathSuffStatArray(GetNsite());
        componentpathsuffstatarray = new PathSuffStatArray(3);
        siteomegapathsuffstatarray = new OmegaPathSuffStatArray(GetNsite());
    }


    //-------------------
    // Accessors
    // ------------------

    //! number of aligned positions
    int GetNsite() const { return codondata->GetNsite(); }

    //! const access to codon state space
    const CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return value of omega_0 < 1
    double GetPurOm() const { return purom; }

    //! return value of omega_2 > 1
    double GetPosOm() const { return 1.0 + dposom; }

    //! return value of dposom  = omega_2 - 1 > 0
    double GetDPosOm() const { return dposom; }

    //! return proportion of sites under strictly purifying selection (under
    //! omega_0 < 1)
    double GetPurW() const { return purw; }

    //! return proportion of sites under positive selection
    double GetPosW() const { return posw; }

    //! whether branch lengths are fixed externally (e.g. when branch lengths are
    //! shared across genes in a multi-gene context)
    bool FixedBranchLengths() const { return blmode == 2; }

    //! whether nuc rates are fixed externally (e.g. when nuc rates are shared
    //! across genes in a multi-gene context)
    bool FixedNucRates() const { return nucmode == 2; }

    //-------------------
    // Setting and updating
    // ------------------

    // Setting model features and (hyper) parameters

    //! \brief set estimation method for branch lengths and nuc rates
    //!
    //! Used in a multigene context.
    //! - mode == 2: global
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode) {
        blmode = inblmode;
        nucmode = innucmode;
    }

    //! set branch lengths to a new value (multi-gene analyses)
    void SetBranchLengths(const BranchSelector<double> &inbranchlength) {
        branchlength->Copy(inbranchlength);
    }

    //! get a copy of branch lengths into array given as argument
    void GetBranchLengths(BranchArray<double> &inbranchlength) const    {
        inbranchlength.Copy(*branchlength);
    }

    //! set branch lengths hyperparameters to a new value (multi-gene analyses)
    void SetBranchLengthsHyperParameters(const BranchSelector<double> &inblmean,
                                         double inblinvshape)   {
        blhypermean->Copy(inblmean);
        blhyperinvshape = inblinvshape;
        branchlength->SetShape(1.0 / blhyperinvshape);
        // branchlength->ResampleEmptyBranches(*lengthpathsuffstatarray);
    }

    //! resample all branches not conditioned by sequence data from prior (as indicated by lengthpathsuffstats)
    void ResampleEmptyBranches()    {
        branchlength->ResampleEmptyBranches(*lengthpathsuffstatarray);
    }

    //! set nucleotide rates (relative exchangeabilities and eq. frequencies) to a
    //! new value (multi-gene analyses)
    void SetNucRates(const std::vector<double> &innucrelrate, const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

    //! get a copy of nucleotide rates into arrays given as arguments
    void GetNucRates(std::vector<double> &innucrelrate, std::vector<double> &innucstat) const   {
        innucrelrate = nucrelrate;
        innucstat = nucstat;
    }

    //! set nucleotide rates hyperparameters to a new value (multi-gene analyses)
    void SetNucRatesHyperParameters(const std::vector<double> &innucrelratehypercenter,
                                    double innucrelratehyperinvconc,
                                    const std::vector<double> &innucstathypercenter,
                                    double innucstathyperinvconc)   {
        nucrelratehypercenter = innucrelratehypercenter;
        nucrelratehyperinvconc = innucrelratehyperinvconc;
        nucstathypercenter = innucstathypercenter;
        nucstathyperinvconc = innucstathyperinvconc;
    }

    //! set omega mixture parameters to a new value
    void SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw)    {
        purom = inpurom;
        dposom = indposom;
        purw = inpurw;
        posw = inposw;
        componentomegaarray->SetParameters(purom, dposom + 1, purw, posw);
    }

    //! shrink omega mixture parameters to a new value
    void ShrinkMixtureParameters(double shrinkposw, double shrinkdposom)    {
        SetMixtureParameters(purom, dposom*shrinkdposom, purw, posw*shrinkposw);
    }

    //! get omega mixture parameter values
    void GetMixtureParameters(double &inpurom, double &indposom, double &inpurw,
                              double &inposw) const {
        inpurom = purom;
        indposom = dposom;
        inpurw = purw;
        inposw = posw;
    }

    //! set omega mixture hyperparameters to a new value
    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc,
                                   double indposomhypermean, double indposomhyperinvshape,
                                   double inpi, double inpurwhypermean, double inpurwhyperinvconc,
                                   double inposwhypermean, double inposwhyperinvconc)   {
        puromhypermean = inpuromhypermean;
        puromhyperinvconc = inpuromhyperinvconc;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
        pi = inpi;
        purwhypermean = inpurwhypermean;
        purwhyperinvconc = inpurwhyperinvconc;
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;

        if (!pi) {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }

        if (!puromhyperinvconc) {
            purom = puromhypermean;
        }
        if (! dposomhyperinvshape)  {
            dposom = dposomhypermean;
        }
        if (! purwhyperinvconc) {
            purw = purwhypermean;
        }
        if (! poswhyperinvconc) {
            posw = poswhypermean;
        }
    }

    //-------------------
    // Matrices
    //-------------------

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() override  {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        componentomegaarray->SetParameters(purom, dposom + 1, purw, posw);
        UpdateMatrices();
        GetIntegratedLogLikelihood();
        ResampleSub(1.0);
    }


    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name) override {
        if (blmode == 0) {
            blhypermean->SetAllBranches(1.0 / lambda);
        }
        componentomegaarray->SetParameters(purom, dposom + 1, purw, posw);
        UpdateMatrices();
        sitealloc->SampleAlloc();
        phyloprocess->PostPredSample(name + ".ali");
        ofstream os((name + ".truesiteom").c_str());
        TraceSiteOmega(os);
        ofstream pos((name + ".trueparam").c_str());
        ToStreamHeader(pos);
        ToStream(pos);
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void UpdateNucMatrix()  {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! it should be updated
    //!
    //! The matrices are not directly updated at that step. Instead, corruption is
    //! notified, such that the matrices know that they will have to recalculate
    //! whichever component is requested later on upon demand.
    void UpdateCodonMatrices()  {
        componentcodonmatrixarray->UpdateCodonMatrices(); 
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to UpdateNucMatrix() and then UpdateCodonMatrices();
    void UpdateMatrices()   {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    void UpdateOmegaMixture()   {
        componentomegaarray->SetParameters(purom, dposom + 1, purw, posw);
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //-------------------
    // Traces and Monitors
    // ------------------

    //! brief return current mean omega value
    double GetMeanOmega() const {
        return posw * (1 + dposom) + (1 - posw) * (purw * purom + (1 - purw));
    }

    void TraceHeader(ostream &os) const override    {
        os << "#logprior\tlnL\tlength\t";
        os << "purom\tposom\tpurw\tposw\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const override  {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        // os << GetIntegratedLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << purom << '\t' << dposom + 1 << '\t' << purw << '\t' << posw << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    //! \brief write current site post probs (of being under positive selection)
    //! on one line
    void TracePostProb(ostream &os) const   {
        for (int i = 0; i < GetNsite(); i++) {
            os << sitepostprobarray[i][2] << '\t';
        }
        os << '\n';
    }

    //! \brief write current site omega values implied by mixture model
    //! on one line
    void TraceSiteOmega(ostream &os) const  {
        for (int i = 0; i < GetNsite(); i++) {
            os << componentomegaarray->GetVal(sitealloc->GetVal(i)) << '\t';
        }
        os << '\n';
    }

    //! \brief get a copy of current site post probs (of being under positive
    //! selection) into array
    void GetSitesPostProb(double *array) const  {
        for (int i = 0; i < GetNsite(); i++) {
            array[i] = sitepostprobarray[i][2];
            if (sitepostprobarray[i][2] < 0) {
                cerr << "error in CodonM2aModel::GetSitesPostProb: negative prob\n";
                exit(1);
            }
        }
    }

    void FromStream(istream &is) override   {
        is >> purom >> dposom >> purw >> posw;
        is >> nucrelrate;
        is >> nucstat;
        is >> lambda;
        is >> *branchlength;
    }

    void ToStream(ostream &os) const override   {
        os << purom << '\t' << dposom << '\t' << purw << '\t' << posw << '\t';
        os << nucrelrate << '\t';
        os << nucstat << '\t';
        os << lambda << '\t';
        os << *branchlength << '\n';
    }

    void ToStreamHeader(ostream &os) const override {
        os << "purom" << '\t' << "dposom" << '\t' << "purw" << '\t' << "posw" << '\t';
        os << "A2C\tA2G\tA2T\tC2G\tC2T\tG2T\t";
        os << "piA\tpiC\tpiG\tpiT\t";
        os << "lambda" << '\t';
        for (int i=0; i<Nbranch; i++)  {
            os << "bl" << i << '\t';
        }
        os << '\n';
    }

    //-------------------
    // Likelihood
    //-------------------

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override {
        return GetLogPrior() + GetLogLikelihood();
        // return GetLogPrior() + GetIntegratedLogLikelihood();
    }

    //! return current value of likelihood, conditional on omega mixture
    //! allocations
    double GetLogLikelihood() const {
        return phyloprocess->GetLogLikelihood();
    }

    //! return current value of likelihood, averaged over omega mixture
    //! allocations
    double GetIntegratedLogLikelihood() const   {
        int ncat = 3;

        double total = 0;
        double logp[ncat];
        const vector<double> &w = componentomegaarray->GetWeights();
        double max = 0;
        for (int i = 0; i < GetNsite(); i++) {
            // int bkalloc = sitealloc->GetVal(i);

            for (int k = 0; k < ncat; k++) {
                (*sitealloc)[i] = k;
                logp[k] = phyloprocess->SiteLogLikelihood(i);
                if ((!k) || (max < logp[k])) {
                    max = logp[k];
                }
            }

            double p = 0;
            for (int k = 0; k < ncat; k++) {
                double tmp = w[k] * exp(logp[k] - max);
                p += tmp;
                sitepostprobarray[i][k] = tmp;
            }
            double logl = log(p) + max;
            total += logl;
            for (int k = 0; k < ncat; k++) {
                sitepostprobarray[i][k] /= p;
            }

            // (*sitealloc)[i] = bkalloc;
        }

        sitealloc->GibbsResample(sitepostprobarray);
        return total;
    }

    //-------------------
    // Priors
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const  {
        double total = 0;

        if (!FixedBranchLengths()) {
            total += BranchLengthsLogPrior();
        }
        if (!FixedNucRates()) {
            total += NucRatesLogPrior();
        }
        total += OmegaLogPrior();
        return total;
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here,
    //! lambda ~ exponential of rate 10)
    double LambdaHyperLogPrior() const  {
        return -lambda / 10; 
    }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const    {
        double total = 0;
        if (blmode == 0) {
            total += LambdaHyperLogPrior();
        }
        total += branchlength->GetLogProb();
        return total;
    }

    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq.
    //! freqs. (nucstat) -- uniform Dirichlet in both cases
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate, nucrelratehypercenter,
                                             1.0 / nucrelratehyperinvconc);
        total += Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //! log prior over omega mixture
    double OmegaLogPrior() const    {
        double total = 0;
        total += PurOmegaLogPrior();
        total += PosOmegaLogPrior();
        total += PurWeightLogPrior();
        total += PosWeightLogPrior();
        return total;
    }

    //! beta prior for purom
    double PurOmegaLogPrior() const {
        if (! puromhyperinvconc)    {
            return 0;
        }
        double alpha = puromhypermean / puromhyperinvconc;
        double beta = (1 - puromhypermean) / puromhyperinvconc;
        return Random::logBetaDensity(purom, alpha, beta);
    }

    //! gamma prior for dposom
    double PosOmegaLogPrior() const {
        if (! dposomhyperinvshape)  {
            return 0;
        }
        double alpha = 1.0 / dposomhyperinvshape;
        double beta = alpha / dposomhypermean;
        return Random::logGammaDensity(dposom, alpha, beta);
    }

    //! beta prior for purw
    double PurWeightLogPrior() const    {
        if (! purwhyperinvconc) {
            return 0;
        }
        double alpha = purwhypermean / purwhyperinvconc;
        double beta = (1 - purwhypermean) / purwhyperinvconc;
        return Random::logBetaDensity(purw, alpha, beta);
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

    //! Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogPrior() const    {
        if (posw) {
            return log(pi);
        }
        return log(1 - pi);
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current
    //! substitution mapping
    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const    {
        return lengthpathsuffstatarray;
    }

    //! \brief const acess to nuc-pathsuffstat
    //!
    //! Useful for resampling nucleotide relative exchangeabilities (nucrelrate)
    //! and equilibrium frequencies (nucstat) conditional on the current
    //! substitution mapping.
    const NucPathSuffStat &GetNucPathSuffStat() const   {
        return nucpathsuffstat;
    }

    //! \brief return log prob of the current substitution mapping, as a function
    //! of the current codon substitution process
    //!
    //! Calculated using pathsuffstat (which summarizes all information about the
    //! substitution mapping) and the codonmatrix. Both pathsuffstat and
    //! codonmatrix are assumed to be updated.
    double PathSuffStatLogProb() const  {
        return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
    }

    //! \brief return log prob of current branch lengths, as a function of branch
    //! lengths hyperparameter lambda
    double LambdaHyperSuffStatLogProb() const   {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    //! \brief return log prob of current substitution mapping, as a function of
    //! nucleotide parameters (nucrelrate and nucstat)
    //!
    //! Calculated using nucpathsuffstat
    //! (which summarizes all information about how the probability of the
    //! substitution mapping depends on nucleotide mutation rates) and the
    //! nucmatrix. Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesSuffStatLogProb() const  {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    //! \brief return log prob of current substitution mapping, as a function of
    //! omega mixture configuration
    //!
    //! Calculated using siteomegapathsuffstatarray
    double OmegaPathSuffStatLogProb() const {
        return componentomegaarray->GetPostProbArray(*siteomegapathsuffstatarray, sitepostprobarray);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    //! \brief log prob factor to be recomputed when moving branch lengths
    //! hyperparameters (here, lambda)
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    //! \brief log prob factor to be recomputed when moving parameters of omega
    //! mixture
    double OmegaLogProb() const { return OmegaLogPrior() + OmegaPathSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() override  {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac)   {
        // UpdateMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep)   {
        for (int rep = 0; rep < nrep; rep++) {
            if (!FixedBranchLengths()) {
                MoveBranchLengths();
            }

            CollectPathSuffStat();

            MoveOmega();

            if (!FixedNucRates()) {
                MoveNucRates();
            }
        }
    }

    //
    // Branch Lengths and hyperparam lambda
    //

    //! overall schedule branch length updatdes
    void MoveBranchLengths()    {
        ResampleBranchLengths();
        if (blmode == 0) {
            MoveLambda();
        }
    }

    //! Gibbs resample branch lengths (based on sufficient statistics and current
    //! value of lambda)
    void ResampleBranchLengths()    {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! collect sufficient statistics for moving branch lengths (directly from the
    //! substitution mappings)
    void CollectLengthSuffStat()    {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda()   {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &CodonM2aModel::LambdaHyperLogProb, &CodonM2aModel::NoUpdate,
                    this);
        ScalingMove(lambda, 0.3, 10, &CodonM2aModel::LambdaHyperLogProb, &CodonM2aModel::NoUpdate,
                    this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat()  {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    //
    // Omega mixture
    //

    //! collect sufficient statistics per component of omega mixtures
    void CollectComponentPathSuffStat() {
        componentpathsuffstatarray->Clear();
        componentpathsuffstatarray->Add(*sitepathsuffstatarray, *sitealloc);
    }

    //! complete move schedule for omega mixture parameters
    void MoveOmega()    {
        CollectOmegaPathSuffStat();

        if (puromhyperinvconc)  {
            SlidingMove(purom, 0.1, 10, 0, 1, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture, this);
        }
        if (purwhyperinvconc)   {
            SlidingMove(purw, 1.0, 10, 0, 1, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture, this);
        }
        if (pi != 0) {
            if (dposomhyperinvshape)    {
                ScalingMove(dposom, 1.0, 10, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture, this);
            }
            if (poswhyperinvconc)   {
                SlidingMove(posw, 1.0, 10, 0, 1, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture,
                        this);
            }
        }
        if ((pi != 0) && (pi != 1)) {
            SwitchPosWeight(10);
        }
        ResampleAlloc();
        UpdateCodonMatrices();
    }

    //! collect sufficient statistics as a function of omega (per site)
    void CollectOmegaPathSuffStat() {
        siteomegapathsuffstatarray->Clear();
        siteomegapathsuffstatarray->AddSuffStat(*sitecodonmatrixarray, *sitepathsuffstatarray);
    }

    //! resample site allocations of omega mixture
    void ResampleAlloc()    {
        OmegaPathSuffStatLogProb();
        sitealloc->GibbsResample(sitepostprobarray);
    }

    //! reversible jump move on posw
    double SwitchPosWeight(int nrep)    {
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
            UpdateOmegaMixture();
            deltalogprob += PosSwitchLogPrior() + OmegaPathSuffStatLogProb();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                posw = bkposw;
                UpdateOmegaMixture();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    //
    // nucleotide parameters
    //

    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat()   {
        UpdateMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray, *componentpathsuffstatarray);
    }

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectComponentPathSuffStat();
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &CodonM2aModel::NucRatesLogProb,
                    &CodonM2aModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &CodonM2aModel::NucRatesLogProb,
                    &CodonM2aModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &CodonM2aModel::NucRatesLogProb,
                    &CodonM2aModel::UpdateNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &CodonM2aModel::NucRatesLogProb,
                    &CodonM2aModel::UpdateNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &CodonM2aModel::NucRatesLogProb,
                    &CodonM2aModel::UpdateNucMatrix, this);

        UpdateMatrices();
    }

    void SetLengthsFromTree() {
        cerr << "before set lengths: " << branchlength->GetTotalLength() << '\n';
        RecursiveSetLengthsFromTree(tree->GetRoot());
        cerr << "after set lengths: " << branchlength->GetTotalLength() << '\n';
    }

    void RecursiveSetLengthsFromTree(const Link *from) {
        if (!from->isRoot()) {
            double tmp = atof(from->GetBranch()->GetName().c_str());
            if (tmp <= 0) {
                cerr << "error: branch length is not positive: " << tmp << '\n';
                exit(1);
            }
            (*branchlength)[from->GetBranch()->GetIndex()] = tmp;
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveSetLengthsFromTree(link->Out());
        }
    }

    void FromStreamCodeML(istream &is)  {
        is >> purom >> dposom >> purw >> posw;
        cerr << "purw: " << purw << '\n';
        cerr << "posw: " << posw << '\n';
        cerr << "w0  : " << purw * (1 - posw) << '\n';
        double kappa;
        is >> kappa;
        double tot = 4 + 2 * kappa;
        nucrelrate[0] = 1.0 / tot;
        nucrelrate[1] = kappa / tot;
        nucrelrate[2] = 1.0 / tot;
        nucrelrate[3] = 1.0 / tot;
        nucrelrate[4] = kappa / tot;
        nucrelrate[5] = 1.0 / tot;
        cerr << "in FromStreamCodeML\n";
        exit(1);
        // nucstat = data->GetEmpiricalFreq();
    }

    //-------------------
    // Data structures
    // ------------------

  private:
    const Tree *tree;
    SequenceAlignment *data;
    const TaxonSet *taxonset;
    const CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    double lambda;
    BranchIIDGamma *blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise *branchlength;

    //
    // parameters of the distribution of omega across sites
    //

    // omega0 < 1: weight purw * (1 - posw)
    // omega1 = 1: weight (1-purw) * (1-posw)
    // omega2 > 1: weight posw

    // omega0 = purom
    double purom;

    // omega2 = 1 + dposom
    double dposom;

    double posw;
    double purw;

    M2aMix *componentomegaarray;
    MultinomialAllocationVector *sitealloc;
    mutable vector<vector<double>> sitepostprobarray;

    //
    // hyperparameters of the priors over the mixture parameters
    //

    // prior probability for the gene to be under positive selection (i.e. prior
    // prob that posw > 0)
    double pi;

    // Beta prior for purom (with hypermean and hyper inverse concentration)
    double puromhypermean;
    double puromhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape
    // parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    // Beta prior for purw
    double purwhypermean;
    double purwhyperinvconc;

    // Beta prior for posw (assuming posw>0)
    double poswhypermean;
    double poswhyperinvconc;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    // nucleotide rate parameters
    vector<double> nucrelrate;
    vector<double> nucstat;
    GTRSubMatrix *nucmatrix;

    // array of matrices across components of the mixture
    MGOmegaCodonSubMatrixArray *componentcodonmatrixarray;

    // arrays of matrices across sites (such as determined by the site allocations
    // to the mixture components) two versions differing only by their exact type

    // used for collecting omega suffstats: need to have access to the *codon*
    // matrix for each site
    MixtureSelector<MGOmegaCodonSubMatrix> *sitecodonmatrixarray;

    // used by PhyloProcess: has to be a Selector<SubMatrix>
    MixtureSelector<SubMatrix> *sitesubmatrixarray;

    PhyloProcess *phyloprocess;

    // suffstats

    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStat hyperlengthsuffstat;
    OmegaPathSuffStatArray *siteomegapathsuffstatarray;
    PathSuffStatArray *sitepathsuffstatarray;
    PathSuffStatArray *componentpathsuffstatarray;

    NucPathSuffStat nucpathsuffstat;

    int blmode;
    int nucmode;
};
