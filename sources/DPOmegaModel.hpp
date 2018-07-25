
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "MultinomialAllocationVector.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "StickBreakingProcess.hpp"
#include "Tree.hpp"

class DPOmegaModel : public ProbModel {
    // tree and data
    Tree *tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;
    int Ncat;  // number of categories of the stick breaking

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

    // omega across sites: Dirichlet process
    // made in three steps

    // (1) component weights (truncated stick breaking)
    double kappa;
    StickBreakingProcess *weight;
    OccupancySuffStat *occupancy;

    // (2) component values: independent omega values from Gamma distribution
    // of mean omegamean and inverse shape parameter omegainvshape
    double omegamean;
    double omegainvshape;
    IIDGamma *componentomegaarray;

    // (3) site allocations: multinomial given the weights
    // multinomial allocation of sites to components of omega distribution
    MultinomialAllocationVector *sitealloc;

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
    // per component
    OmegaPathSuffStatArray *componentomegapathsuffstatarray;

    // for moving omegamean and invshape
    GammaSuffStat omegahypersuffstat;

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

    DPOmegaModel(string datafile, string treefile, int inNcat) {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);
        Ncat = inNcat;

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        if (Ncat == -1) {
            Ncat = Nsite;
            if (Ncat > 100) {
                Ncat = 100;
            }
        }

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

        kappa = 1.0;
        weight = new StickBreakingProcess(Ncat, kappa);
        occupancy = new OccupancySuffStat(Ncat);

        sitealloc = new MultinomialAllocationVector(Nsite, weight->GetArray());

        omegamean = 1.0;
        omegainvshape = 1.0;
        componentomegaarray = new IIDGamma(Ncat, omegamean, omegainvshape);

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
        componentomegapathsuffstatarray = new OmegaPathSuffStatArray(Ncat);
    }

    void Unfold() {
        cerr << "-- unfold\n";
        phyloprocess->Unfold();
        cerr << phyloprocess->GetLogLikelihood() << '\n';
        std::cerr << "-- mapping substitutions\n";
        phyloprocess->ResampleSub();
    }

    void Update() {
        UpdateOccupancies();
        UpdateMatrices();
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

    void UpdateCodonMatrices() { componentcodonmatrixarray->UpdateCodonMatrices(*occupancy); }

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
        total += StickBreakingHyperLogPrior();
        total += StickBreakingLogPrior();
        total += OmegaHyperLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    // conditional on site allocations
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    double BranchLengthsHyperLogPrior() const {
        // exponential of mean 10
        return -lambda / 10;
    }

    double BranchLengthsLogPrior() const { return branchlength->GetLogProb(); }

    // uniform prior
    double NucRatesLogPrior() const { return 0; }

    double StickBreakingHyperLogPrior() const { return -kappa / 10; }

    double StickBreakingLogPrior() const { return weight->GetLogProb(kappa); }

    double OmegaHyperLogPrior() const {
        double total = 0;
        total -= omegamean;
        total -= omegainvshape;
        return total;
    }

    double OmegaLogPrior() const { return componentomegaarray->GetLogProb(); }

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

    double OmegaHyperSuffStatLogProb() const {
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        return omegahypersuffstat.GetLogProb(alpha, beta);
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
    double OmegaHyperLogProb() const { return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb(); }

    double StickBreakingHyperLogProb() const {
        return StickBreakingHyperLogPrior() + StickBreakingLogPrior();
    }

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

            CollectSitePathSuffStat();

            MoveOmega();
            MoveOmegaHyperParameters();
            MoveKappa();
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
        ScalingMove(lambda, 1.0, 10, &DPOmegaModel::BranchLengthsHyperLogProb,
                    &DPOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &DPOmegaModel::BranchLengthsHyperLogProb,
                    &DPOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    // per site
    void CollectSitePathSuffStat() {
        sitepathsuffstatarray->Clear();
        sitepathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    // per component of the mixture
    void CollectComponentPathSuffStat() {
        componentpathsuffstatarray->Clear();
        componentpathsuffstatarray->Add(*sitepathsuffstatarray, *sitealloc);
    }

    void MoveOmega() {
        CollectSiteOmegaPathSuffStat();
        for (int rep = 0; rep < 10; rep++) {
            ResampleOmega();
            ResampleAlloc();
            LabelSwitchingMove();
            ResampleWeights();
        }
    }

    void UpdateOccupancies() {
        occupancy->Clear();
        occupancy->AddSuffStat(*sitealloc);
    }

    void CollectSiteOmegaPathSuffStat() {
        siteomegapathsuffstatarray->Clear();
        siteomegapathsuffstatarray->AddSuffStat(*sitecodonmatrixarray, *sitepathsuffstatarray);
    }

    void CollectComponentOmegaPathSuffStat() {
        componentomegapathsuffstatarray->Clear();
        componentomegapathsuffstatarray->Add(*siteomegapathsuffstatarray, *sitealloc);
    }

    void ResampleOmega() {
        CollectComponentOmegaPathSuffStat();
        componentomegaarray->GibbsResample(*componentomegapathsuffstatarray);
    }

    void ResampleAlloc() {
        vector<double> postprob(Ncat, 0);
        for (int i = 0; i < Nsite; i++) {
            componentomegaarray->GetAllocPostProb(siteomegapathsuffstatarray->GetVal(i),
                                                  weight->GetArray(), postprob);
            sitealloc->GibbsResample(i, postprob);
        }
        UpdateOccupancies();
    }

    void LabelSwitchingMove() {
        MoveOccupiedCompAlloc(5);
        MoveAdjacentCompAlloc(5);
    }

    double MoveOccupiedCompAlloc(int k0) {
        const vector<double> &w = weight->GetArray();

        int nrep = (int)(k0 * kappa);
        ResampleWeights();
        double total = 0.0;
        int Nocc = GetNcluster();
        if (Nocc != 1) {
            for (int i = 0; i < nrep; i++) {
                int occupiedComponentIndices[Nocc];
                int j = 0;
                for (int k = 0; k < Ncat; k++) {
                    if ((*occupancy)[k] != 0) {
                        occupiedComponentIndices[j] = k;
                        j++;
                    }
                }
                if (j != Nocc) {
                    cerr << "error in MoveOccupiedCompAlloc.\n";
                    exit(1);
                }
                int indices[2];
                Random::DrawFromUrn(indices, 2, Nocc);
                int cat1 = occupiedComponentIndices[indices[0]];
                int cat2 = occupiedComponentIndices[indices[1]];
                double logMetropolis =
                    ((*occupancy)[cat2] - (*occupancy)[cat1]) * log(w[cat1] / w[cat2]);
                int accepted = (log(Random::Uniform()) < logMetropolis);
                if (accepted) {
                    total += 1.0;
                    componentomegaarray->Swap(cat1, cat2);
                    sitealloc->SwapComponents(cat1, cat2);
                    occupancy->Swap(cat1, cat2);
                }
            }
            return total /= nrep;
        }
        return 0;
    }

    double MoveAdjacentCompAlloc(int k0) {
        ResampleWeights();
        int nrep = (int)(k0 * kappa);

        double total = 0;

        const vector<double> &V = weight->GetBetaVariates();

        for (int i = 0; i < nrep; i++) {
            int cat1 = (int)(Random::Uniform() * (Ncat - 2));
            int cat2 = cat1 + 1;
            double logMetropolis =
                ((*occupancy)[cat1] * log(1 - V[cat2])) - ((*occupancy)[cat2] * log(1 - V[cat1]));
            int accepted = (log(Random::Uniform()) < logMetropolis);
            if (accepted) {
                total += 1.0;
                componentomegaarray->Swap(cat1, cat2);
                sitealloc->SwapComponents(cat1, cat2);
                weight->SwapComponents(cat1, cat2);
                occupancy->Swap(cat1, cat2);
            }
        }

        return total /= nrep;
    }

    void ResampleWeights() { weight->GibbsResample(*occupancy); }

    void MoveOmegaHyperParameters() {
        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*componentomegaarray, *occupancy);
        ScalingMove(omegamean, 1.0, 10, &DPOmegaModel::OmegaHyperLogProb, &DPOmegaModel::NoUpdate,
                    this);
        ScalingMove(omegamean, 0.3, 10, &DPOmegaModel::OmegaHyperLogProb, &DPOmegaModel::NoUpdate,
                    this);
        ScalingMove(omegainvshape, 1.0, 10, &DPOmegaModel::OmegaHyperLogProb,
                    &DPOmegaModel::NoUpdate, this);
        ScalingMove(omegainvshape, 0.3, 10, &DPOmegaModel::OmegaHyperLogProb,
                    &DPOmegaModel::NoUpdate, this);
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        componentomegaarray->SetShape(alpha);
        componentomegaarray->SetScale(beta);
        componentomegaarray->PriorResample(*occupancy);
    }

    void MoveKappa() {
        ScalingMove(kappa, 1.0, 10, &DPOmegaModel::StickBreakingHyperLogProb,
                    &DPOmegaModel::NoUpdate, this);
        ScalingMove(kappa, 0.3, 10, &DPOmegaModel::StickBreakingHyperLogProb,
                    &DPOmegaModel::NoUpdate, this);
    }

    void CollectNucPathSuffStat() {
        UpdateMatrices();
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray, *componentpathsuffstatarray);
    }

    void MoveNucRates() {
        CollectComponentPathSuffStat();
        CollectNucPathSuffStat();

        ProfileMove(nucrelrate, 0.1, 1, 3, &DPOmegaModel::NucRatesLogProb,
                    &DPOmegaModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &DPOmegaModel::NucRatesLogProb,
                    &DPOmegaModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &DPOmegaModel::NucRatesLogProb,
                    &DPOmegaModel::UpdateNucMatrix, this);

        ProfileMove(nucstat, 0.1, 1, 3, &DPOmegaModel::NucRatesLogProb,
                    &DPOmegaModel::UpdateNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 3, &DPOmegaModel::NucRatesLogProb,
                    &DPOmegaModel::UpdateNucMatrix, this);

        UpdateMatrices();
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    int GetNcluster() const {
        int n = 0;
        for (int i = 0; i < Ncat; i++) {
            if (occupancy->GetVal(i)) {
                n++;
            }
        }
        return n;
    }

    double GetMeanOmega() const {
        double mean = 0;
        for (int i = 0; i < Ncat; i++) {
            mean += occupancy->GetVal(i) * componentomegaarray->GetVal(i);
        }
        mean /= Nsite;
        return mean;
    }

    double GetOmegaRelVar() const {
        double mean = 0;
        double var = 0;
        for (int i = 0; i < Ncat; i++) {
            mean += occupancy->GetVal(i) * componentomegaarray->GetVal(i);
            var += occupancy->GetVal(i) * componentomegaarray->GetVal(i) *
                   componentomegaarray->GetVal(i);
        }
        mean /= Nsite;
        var /= Nsite;
        var -= mean * mean;
        var /= mean * mean;
        return var;
    }

    double GetEmpiricalPosFrac() const {
        double tot = 0;
        for (int i = 0; i < Nsite; i++) {
            if ((*componentomegaarray)[(*sitealloc)[i]] > 1) {
                tot++;
            }
        }
        return tot / Nsite;
    }

    void TraceHeader(std::ostream &os) const {
        os << "#logprior\tlnL\tlength\t";
        os << "meanomega\t";
        os << "relvar\t";
        os << "posfrac\t";
        os << "ncluster\t";
        os << "omegamean\tinvshape\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream &os) const {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << GetMeanOmega() << '\t';
        os << GetOmegaRelVar() << '\t';
        os << GetEmpiricalPosFrac() << '\t';
        os << GetNcluster() << '\t';
        os << omegamean << '\t';
        os << omegainvshape << '\t';
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void TraceSiteOmega(ostream &os) const {
        for (int i = 0; i < Nsite; i++) {
            os << (*componentomegaarray)[(*sitealloc)[i]] << '\t';
        }
        os << '\n';
        os.flush();
    }

    void Monitor(ostream &os) const {}

    void ToStream(ostream &os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << kappa << '\n';
        os << omegamean << '\t' << omegainvshape << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

    void FromStream(istream &is) {
        is >> lambda;
        is >> *branchlength;
        is >> kappa;
        is >> omegamean >> omegainvshape;
        is >> nucrelrate;
        is >> nucstat;
    }
};
