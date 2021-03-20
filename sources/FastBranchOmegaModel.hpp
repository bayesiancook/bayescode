#include "WhiteNoise.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "ProbModel.hpp"
#include "Product.hpp"
#include "SimpleSubMatrixSelector.hpp"
#include "IIDGamma.hpp"
#include "Tree.hpp"
#include "CodonSuffStat.hpp"
#include "dSOmegaPathSuffStat.hpp"

class FastBranchOmegaModel : public ProbModel {

    // tree and data
    const Tree *tree;
    const TaxonSet *taxonset;
    string dsomsuffstatfile;

    int Ntaxa;
    int Nbranch;

    // Branch lengths

    double lambda;
    BranchIIDGamma *blhypermean;
    double blhyperinvshape;
    GammaWhiteNoise *branchlength;

    double omegamean;
    double omegainvshape;
    BranchIIDGamma *omegabrancharray;

    dSOmegaPathSuffStatBranchArray* dsompathsuffstatarray;

    PoissonSuffStatBranchArray *lengthpathsuffstatbrancharray;
    GammaSuffStat hyperlengthsuffstat;

    OmegaPathSuffStatBranchArray *omegapathsuffstatbrancharray;
    GammaSuffStat omegahypersuffstat;

  public:
    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    FastBranchOmegaModel(string taxonfile, string treefile, string indsomsuffstatfile)  {

        taxonset = new TaxonSet(taxonfile);
        Ntaxa = taxonset->GetNtaxa();

        // get tree from file (newick format)
        Tree* tmptree = new Tree(treefile);
        // check whether tree and data fits together
        tmptree->RegisterWith(taxonset);
        tmptree->SetIndices();
        tree = tmptree;

        Nbranch = tree->GetNbranch();

        dsomsuffstatfile = indsomsuffstatfile;
    }


    //! model allocation
    void Allocate() {

        // Branch lengths

        lambda = 10.0;
        blhypermean = new BranchIIDGamma(*tree, 1.0, lambda);
        blhypermean->SetAllBranches(1.0 / lambda);
        blhyperinvshape = 1.0;
        branchlength = new GammaWhiteNoise(*tree, *blhypermean, 1.0 / blhyperinvshape);

        // omega across branches

        omegamean = 0.1;
        omegainvshape = 0.3;
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegabrancharray = new BranchIIDGamma(*tree, alpha, beta);

        dsompathsuffstatarray = new dSOmegaPathSuffStatBranchArray(*tree);
        ifstream is(dsomsuffstatfile.c_str());
        is >> *dsompathsuffstatarray;

        lengthpathsuffstatbrancharray = new PoissonSuffStatBranchArray(*tree);
        omegapathsuffstatbrancharray = new OmegaPathSuffStatBranchArray(*tree);
    }

    //-------------------
    // Accessors
    // ------------------

    const Tree& GetTree() const {return *tree;}

    //-------------------
    // Setting and updating
    // ------------------

    void Update() override {
        blhypermean->SetAllBranches(1.0 / lambda);
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegabrancharray->SetShape(alpha);
        omegabrancharray->SetScale(beta);
    }

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(string name) override {
        Update();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in ProbModel), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        total += BranchLengthsLogPrior();
        total += OmegaLogPrior();
        total += OmegaHyperLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { 
        return dSOmPathSuffStatLogProb();
    }

    double dSOmPathSuffStatLogProb() const {
        double total = 0;
        for (int index=0; index<tree->GetNbranch(); index++)    {
            total += dsompathsuffstatarray->GetVal(index).GetLogProb(branchlength->GetVal(index), omegabrancharray->GetVal(index));
        }
        return total;
    }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override { return GetLogPrior() + GetLogLikelihood(); }

    //! log prior over branch lengths (iid exponential of rate lambda)
    double BranchLengthsLogPrior() const {
        double total = 0;
        total += LambdaHyperLogPrior();
        total += branchlength->GetLogProb();
        return total;
    }

    //! \brief log prior over hyperparameter of prior over branch lengths (here,
    //! lambda ~ exponential of rate 10)
    double LambdaHyperLogPrior() const { return -lambda / 10; }

    //! log prior over omega
    double OmegaLogPrior() const { return omegabrancharray->GetLogProb(); }

    double OmegaHyperLogPrior() const   {
        double total = 0;
        total -= omegamean;
        total -= omegainvshape;
        return total;
    }

    void UpdateOmega()  {
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegabrancharray->SetShape(alpha);
        omegabrancharray->SetScale(beta);
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

    double OmegaHyperSuffStatLogProb() const {
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        return omegahypersuffstat.GetLogProb(alpha, beta);
    }

    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    double OmegaHyperLogProb() const    {
        return OmegaHyperLogPrior() + OmegaLogPrior();
    }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() override {
        MoveParameters(30);
        return 1.0;
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveBranchLengths();
            MoveOmega();
        }
    }

    //! overall schedule branch length updatdes
    void MoveBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatbrancharray);
        MoveLambda();
    }

    void CollectLengthSuffStat() {
        lengthpathsuffstatbrancharray->Clear();
        dsompathsuffstatarray->TodSSuffStat(*lengthpathsuffstatbrancharray, *omegabrancharray);
    }

    //! MH move on branch lengths hyperparameters (here, scaling move on lambda,
    //! based on suffstats for branch lengths)
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &FastBranchOmegaModel::LambdaHyperLogProb, &FastBranchOmegaModel::NoUpdate,
                    this);
        ScalingMove(lambda, 0.3, 10, &FastBranchOmegaModel::LambdaHyperLogProb, &FastBranchOmegaModel::NoUpdate,
                    this);
        blhypermean->SetAllBranches(1.0 / lambda);
    }

    void MoveOmega()    {
        CollectOmegaSuffStat();
        omegabrancharray->GibbsResample(*omegapathsuffstatbrancharray);
        MoveOmegaHyperParameters();
    }

    void CollectOmegaSuffStat() {
        omegapathsuffstatbrancharray->Clear();
        dsompathsuffstatarray->ToOmSuffStat(*omegapathsuffstatbrancharray, *branchlength);
    }

    void MoveOmegaHyperParameters() {
        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegabrancharray);
        ScalingMove(omegamean, 1.0, 10, &FastBranchOmegaModel::OmegaHyperLogProb,
                    &FastBranchOmegaModel::NoUpdate, this);
        ScalingMove(omegamean, 0.3, 10, &FastBranchOmegaModel::OmegaHyperLogProb,
                    &FastBranchOmegaModel::NoUpdate, this);
        ScalingMove(omegainvshape, 1.0, 10, &FastBranchOmegaModel::OmegaHyperLogProb,
                    &FastBranchOmegaModel::NoUpdate, this);
        ScalingMove(omegainvshape, 0.3, 10, &FastBranchOmegaModel::OmegaHyperLogProb,
                    &FastBranchOmegaModel::NoUpdate, this);
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / omegamean;
        omegabrancharray->SetShape(alpha);
        omegabrancharray->SetScale(beta);
    }

    //-------------------
    // Traces and Monitors
    // ------------------

    double GetMeanOmega() const {
        return omegabrancharray->GetMean();
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "meanom\tvarom\n";
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << omegabrancharray->GetMean() << '\t' << omegabrancharray->GetVar() << '\n';
    }

    void TraceBranchOmega(ostream& os) const {
        for (int i=0; i<Nbranch; i++) {
            os << omegabrancharray->GetVal(i) << '\t';
        }
        os << '\n';
        os.flush();
    }

    const BranchSelector<double>& GetOmegaTree() const  {
        return *omegabrancharray;
    }

    void Monitor(ostream &os) const override {}

    void ToStream(ostream &os) const override {
        os << lambda << '\t';
        os << *branchlength << '\t';
        os << omegamean << '\t' << omegainvshape << '\t';
        os << *omegabrancharray << '\n';
    }

    void FromStream(istream &is) override {
        is >> lambda;
        is >> *branchlength;
        is >> omegamean >> omegainvshape;
        is >> *omegabrancharray;
    }
};
