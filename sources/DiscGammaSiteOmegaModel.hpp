
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "DiscGamma.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "MultinomialAllocationVector.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

class DiscGammaSiteOmegaModel : public ProbModel {
    // tree and data
    Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;
    int Ncat;

    // branch lengths iid expo (gamma of shape 1 and scale lambda)
    // where lambda is a hyperparameter
    double lambda;
    BranchIIDGamma *branchlength;

    // nucleotide exchange rates and equilibrium frequencies (stationary
    // probabilities)
    std::vector<double> nucrelrate;
    std::vector<double> nucstat;

    // a nucleotide matrix (parameterized by nucrelrate and nucstat)
    GTRSubMatrix *nucmatrix;

    // omega across sites: discretized Gamma distribution
    // of mean omegamean and inverse shape parameter omegainvshape
    double omegamean;
    double omegainvshape;
    DiscGamma *componentomegaarray;

    // multinomial allocation of sites to components of omega distribution
    MultinomialAllocationVector *sitealloc;
    mutable vector<vector<double>> sitepostprobarray;

    // an array of codon matrices (one for each discrete value of omega)
    MGOmegaCodonSubMatrixArray *componentcodonmatrixarray;

    // 2 arrays of matrix pointers across sites
    // obtained from component matrix array and site allocations
    // the 2 arrays are identical in value, they just differ by their type

    // this one is used by PhyloProcess: has to be a Selector<SubMatrix>
    MixtureSelector<SubMatrix> *sitesubmatrixarray;

    // this one is used for collecting omega suffstats: need to have access to the
    // *codon* matrix for each site
    MixtureSelector<MGOmegaCodonSubMatrix> *sitecodonmatrixarray;

    PhyloProcess *phyloprocess;

    // suffstats

    // generic suff stats for substitution paths
    // per site
    PathSuffStatArray *sitepathsuffstatarray;
    // per component of the mixture
    PathSuffStatArray *componentpathsuffstatarray;

    // which can be collected across all sites and branches
    // and summarized in terms of 4x4 suff stats, as a function of nucleotide
    // rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function
    // of omega per site
    OmegaPathSuffStatArray *siteomegapathsuffstatarray;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    // suff stats for branch lengths, as a function of their hyper parameter
    // lambda (bl are iid gamma, of scale parameter lambda)
    GammaSuffStat hyperlengthsuffstat;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    DiscGammaSiteOmegaModel(string datafile, string treefile, int inNcat) {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);
        Ncat = inNcat;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        std::cerr << "-- Number of sites: " << Nsite << std::endl;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();
    }

    void Allocate() {
        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);

        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));

        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        omegamean = 1.0;
        omegainvshape = 1.0;
        componentomegaarray = new DiscGamma(Ncat, omegamean, omegainvshape);
        sitealloc = new MultinomialAllocationVector(Nsite, componentomegaarray->GetWeights());
        sitepostprobarray.assign(Nsite, vector<double>(Ncat, 0));

        componentcodonmatrixarray = new MGOmegaCodonSubMatrixArray(
            (CodonStateSpace *)codondata->GetStateSpace(), nucmatrix, componentomegaarray);

        sitesubmatrixarray = new MixtureSelector<SubMatrix>(componentcodonmatrixarray, sitealloc);
        sitecodonmatrixarray =
            new MixtureSelector<MGOmegaCodonSubMatrix>(componentcodonmatrixarray, sitealloc);

        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, sitesubmatrixarray);

        // suff stats
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        sitepathsuffstatarray = new PathSuffStatArray(Nsite);
        componentpathsuffstatarray = new PathSuffStatArray(Ncat);
        siteomegapathsuffstatarray = new OmegaPathSuffStatArray(Nsite);
    }

    void Unfold() {
        cerr << "-- unfold\n";
        phyloprocess->Unfold();
        cerr << phyloprocess->GetLogLikelihood() << '\n';
        std::cerr << "-- mapping substitutions\n";
        phyloprocess->ResampleSub();
    }

    void Update() {
        componentomegaarray->SetParameters(omegamean, omegainvshape);
        UpdateMatrices();
        ResampleAlloc();
    }

    //-------------------
    // Accessors
    // ------------------

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //-------------------
    // Setting and updating
    // ------------------

    void SetOmega(double inomegamean, double inomegainvshape) {
        omegamean = inomegamean;
        omegainvshape = inomegainvshape;
        componentomegaarray->SetParameters(omegamean, omegainvshape);
        UpdateMatrices();
    }

    void SetBranchLengths(const BranchSelector<double> &inbranchlength) {
        branchlength->Copy(inbranchlength);
    }

    void SetNucRates(const std::vector<double> &innucrelrate,
                     const std::vector<double> &innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void UpdateCodonMatrices() { componentcodonmatrixarray->UpdateCodonMatrices(); }

    void UpdateMatrices() {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    void NoUpdate() {}

    //-------------------
    // Priors and likelihood
    //-------------------

    double GetLogPrior() const {
        double total = 0;
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    // conditional on site allocations
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    // integrated over site allocations
    double GetIntegratedLogLikelihood() const {
        double total = 0;
        double logp[Ncat];
        const vector<double> &w = componentomegaarray->GetWeights();
        double max = 0;
        for (int i = 0; i < Nsite; i++) {
            int bkalloc = sitealloc->GetVal(i);

            for (int k = 0; k < Ncat; k++) {
                (*sitealloc)[i] = k;
                logp[k] = phyloprocess->SiteLogLikelihood(i);
                if ((!k) || (max < logp[k])) {
                    max = logp[k];
                }
            }

            double p = 0;
            for (int k = 0; k < Ncat; k++) {
                p += w[k] * exp(logp[k] - max);
            }
            double logl = log(p) + max;
            total += logl;

            (*sitealloc)[i] = bkalloc;
        }
        return total;
    }

    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    double BranchLengthsHyperLogPrior() const {
        // exponential of mean 10
        return -lambda / 10;
    }

    double BranchLengthsLogPrior() const { return branchlength->GetLogProb(); }

    // uniform prior
    double NucRatesLogPrior() const { return 0; }

    double OmegaLogPrior() const {
        double total = 0;
        total -= omegamean;
        total -= omegainvshape;
        return total;
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    const NucPathSuffStat &GetNucPathSuffStat() const { return nucpathsuffstat; }

    double PathSuffStatLogProb() const {
        return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
    }

    double BranchLengthsHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    double OmegaPathSuffStatLogProb() const {
        componentomegaarray->SetParameters(omegamean, omegainvshape);
        return componentomegaarray->GetPostProbArray(*siteomegapathsuffstatarray,
                                                     sitepostprobarray);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // for moving branch lengths hyperparameter lambda
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // for moving nuc rates
    double NucRatesLogProb() const { return NucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    // for moving omegamean and invshape
    double OmegaLogProb() const { return OmegaLogPrior() + OmegaPathSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    double Move() {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    void ResampleSub(double frac) {
        UpdateMatrices();
        phyloprocess->Move(frac);
    }

    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();

            CollectPathSuffStat();

            MoveOmega();
            MoveNucRates();
        }
    }

    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    void MoveBranchLengthsHyperParameter() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &DiscGammaSiteOmegaModel::BranchLengthsHyperLogProb,
                    &DiscGammaSiteOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &DiscGammaSiteOmegaModel::BranchLengthsHyperLogProb,
                    &DiscGammaSiteOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    // per site
    void CollectPathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    // per component of the mixture
    void CollectComponentPathSuffStat() {
        componentpathsuffstatarray->Clear();
        componentpathsuffstatarray->Add(*sitepathsuffstatarray, *sitealloc);
    }

    void CollectOmegaPathSuffStat() {
        siteomegapathsuffstatarray->Clear();
        siteomegapathsuffstatarray->AddSuffStat(*sitecodonmatrixarray, *sitepathsuffstatarray);
    }

    void ResampleAlloc() {
        OmegaPathSuffStatLogProb();
        sitealloc->GibbsResample(sitepostprobarray);
    }

    void MoveOmega() {
        CollectOmegaPathSuffStat();

        ScalingMove(omegamean, 1.0, 10, &DiscGammaSiteOmegaModel::OmegaLogProb,
                    &DiscGammaSiteOmegaModel::NoUpdate, this);
        ScalingMove(omegamean, 0.3, 10, &DiscGammaSiteOmegaModel::OmegaLogProb,
                    &DiscGammaSiteOmegaModel::NoUpdate, this);
        ScalingMove(omegainvshape, 1.0, 10, &DiscGammaSiteOmegaModel::OmegaLogProb,
                    &DiscGammaSiteOmegaModel::NoUpdate, this);
        ScalingMove(omegainvshape, 0.3, 10, &DiscGammaSiteOmegaModel::OmegaLogProb,
                    &DiscGammaSiteOmegaModel::NoUpdate, this);

        ResampleAlloc();
    }

    void CollectNucPathSuffStat() {
        UpdateMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray, *componentpathsuffstatarray);
    }

    void MoveNucRates() {
        CollectComponentPathSuffStat();
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &DiscGammaSiteOmegaModel::NucRatesLogProb,
                    &DiscGammaSiteOmegaModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &DiscGammaSiteOmegaModel::NucRatesLogProb,
                    &DiscGammaSiteOmegaModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &DiscGammaSiteOmegaModel::NucRatesLogProb,
                    &DiscGammaSiteOmegaModel::UpdateNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &DiscGammaSiteOmegaModel::NucRatesLogProb,
                    &DiscGammaSiteOmegaModel::UpdateNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &DiscGammaSiteOmegaModel::NucRatesLogProb,
                    &DiscGammaSiteOmegaModel::UpdateNucMatrix, this);

        UpdateMatrices();
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(std::ostream &os) const {
        os << "#logprior\tlnL\tlength\tlambda\t";
        os << "omegamean\tinvshape\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << lambda << '\t';
        os << omegamean << '\t';
        os << omegainvshape << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
        cerr << *componentomegaarray << '\n';
    }

    void Monitor(ostream &os) const {}

    void ToStream(ostream &os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << omegamean << '\t' << omegainvshape << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

    void FromStream(istream &is) {
        is >> lambda;
        is >> *branchlength;
        is >> omegamean >> omegainvshape;
        is >> nucrelrate;
        is >> nucstat;
    }
};
