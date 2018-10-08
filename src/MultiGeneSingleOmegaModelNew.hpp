// this is a multigene version of singleomegamodel
//
// - branch lengths are shared across genes, and are iid Exponential of rate
// lambda
// - nucleotide relative exchangeabilities and stationaries are also shared
// across genes (uniform Dirichlet)
// - the array of gene-specific omega's are iid gamma with hyperparameters
// omegahypermean and omegahyperinvshape
//
// the sequence of MCMC moves is as follows:
// - genes resample substitution histories, gather path suff stats and move
// their omega's
// - master receives the array of omega's across genes, moves their
// hyperparameters and then broadcast the new value of these hyperparams
// - master collects branch length suff stats across genes, moves branch lengths
// and broadcasts their new value
// - master collects nuc path suffstats across genes, moves nuc rates and
// broadcasts their new value

#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "SingleOmegaModel.hpp"
#include "MultiGeneMPIModule.hpp"
#include "components/ChainComponent.hpp"

class MultiGeneSingleOmegaModelShared {
public:
    MultiGeneSingleOmegaModelShared(string datafile, string intreefile)
        : treefile(intreefile), nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {
        blmode = 1;
        nucmode = 1;
        omegamode = 1;
    }

protected:
    std::unique_ptr<const Tree> tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;

    string treefile;

    int Ntaxa;
    int Nbranch;

    int blmode;
    int nucmode;
    int omegamode;

    // Branch lengths

    double lambda;
    BranchIIDGamma *branchlength;
    GammaSuffStat hyperlengthsuffstat;

    double blhyperinvshape;
    GammaWhiteNoiseArray *branchlengtharray;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStatBranchArray *lengthhypersuffstatarray;

    // Nucleotide rates

    // shared nuc rates
    GTRSubMatrix *nucmatrix;
    NucPathSuffStat nucpathsuffstat;

    // gene-specific nuc rates
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    IIDDirichlet *nucrelratearray;
    DirichletSuffStat nucrelratesuffstat;

    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    IIDDirichlet *nucstatarray;
    DirichletSuffStat nucstatsuffstat;

    // Omega

    // each gene has its own omega
    // omegaarray[gene], for gene=1..Ngene
    // iid gamma, with hyperparameters omegahypermean and hyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
    IIDGamma *omegaarray;

    // suffstat for gene-specific omega's
    // as a function of omegahypermean and omegahyperinvshape
    GammaSuffStat omegahypersuffstat;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;
};

class MultiGeneSingleOmegaModelMaster : public MultiGeneSingleOmegaModelShared,
                                        public ProbModel,
                                        public ChainComponent {
    MultiGeneMPIModule mpi;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneSingleOmegaModelMaster(string datafile, string intreefile, int inmyid, int innprocs)
        : MultiGeneSingleOmegaModelShared(datafile, intreefile),
          mpi(inmyid, innprocs) {


        mpi.AllocateAlignments(datafile);

        refcodondata = new CodonSequenceAlignment(mpi.refdata, true);
        taxonset = mpi.refdata->GetTaxonSet();
        Ntaxa = mpi.refdata->GetNtaxa();

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        // // check whether tree and data fits together
        // tree->RegisterWith(taxonset);

        // tree->SetIndices();
        // Nbranch = tree->GetNbranch();

        cerr << "number of taxa : " << Ntaxa << '\n';
        cerr << "number of branches : " << Nbranch << '\n';
        cerr << "tree and data fit together\n";
    }

    void start() override {}
    void move(int) override {
        SendRunningStatus(1);
        Move();
    }
    void savepoint(int) override {}
    void end() override {}

    void SendRunningStatus(int status) {
        MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    void Allocate() {
        // Branch lengths

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
        blhyperinvshape = 0.1;
        if (blmode == 2) {
            lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
            lengthhypersuffstatarray = 0;
        } else {
            branchlength->SetAllBranches(1.0 / lambda);
            branchlengtharray = new GammaWhiteNoiseArray(
                mpi.GetLocalNgene(), *tree, *branchlength, 1.0 / blhyperinvshape);
            lengthpathsuffstatarray = 0;
            lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
        }

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;

        if (nucmode == 2) {
            nucrelratearray =
                new IIDDirichlet(1, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(1, nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = new GTRSubMatrix(Nnuc, (*nucrelratearray)[0], (*nucstatarray)[0], true);
        } else {
            nucrelratearray = new IIDDirichlet(
                mpi.GetLocalNgene(), nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray =
                new IIDDirichlet(mpi.GetLocalNgene(), nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = 0;
        }

        // Omega

        /*
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        */
        omegaarray = new IIDGamma(mpi.GetLocalNgene(), omegahypermean, omegahyperinvshape);

        // Gene processes

        lnL = 0;
        GeneLogPrior = 0;
    }

    // called upon constructing the model
    // mode == 2: global
    // mode == 1: gene specific, with hyperparameters estimated across genes
    // mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode, int inomegamode) {
        blmode = inblmode;
        nucmode = innucmode;
        omegamode = inomegamode;
    }

    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape) {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    void FastUpdate() {
        branchlength->SetScale(lambda);
        if (blmode == 1) { branchlengtharray->SetShape(1.0 / blhyperinvshape); }
        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    void Update() {
        FastUpdate();

        if (mpi.GetNprocs() > 1) {
            SendBranchLengthsHyperParameters();
            SendNucRatesHyperParameters();

            if (blmode == 2) {
                SendGlobalBranchLengths();
            } else {
                SendGeneBranchLengths();
            }

            if (nucmode == 2) {
                SendGlobalNucRates();
            } else {
                SendGeneNucRates();
            }

            SendOmegaHyperParameters();
            SendOmega();
            ReceiveLogProbs();
        }
    }

    void PostPred(string name) {
        FastUpdate();
        if (mpi.GetNprocs() > 1) {
            SendBranchLengthsHyperParameters();
            SendNucRatesHyperParameters();

            if (blmode == 2) {
                SendGlobalBranchLengths();
            } else {
                SendGeneBranchLengths();
            }

            if (nucmode == 2) {
                SendGlobalNucRates();
            } else {
                SendGeneNucRates();
            }

            SendOmegaHyperParameters();
            SendOmega();
        }
    }

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    const vector<double> &GetOmegaArray() const { return omegaarray->GetArray(); }

    //-------------------
    // Traces and Monitors
    //-------------------

    void TraceHeader(ostream &os) const {
        os << "#logprior\tlnL";
        if (blmode == 2) {
            os << "\tlength";
        } else {
            os << "\tmeanlength\tstdev";
        }
        os << "\tmeanomega";
        os << "\tvaromega";
        os << "\tomegahypermean\tinvshape";
        os << "\tstatent";
        os << "\trrent";
        if (nucmode != 2) {
            os << "\tstdevrr\tcenter\thyperinvconc";
            os << "\tstdevstat\tcenter\thyperinvconc";
        }
        os << '\n';
    }

    void Trace(ostream &os) const {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood();

        if (blmode == 2) {
            os << '\t' << GetMeanTotalLength();
        } else {
            os << '\t' << GetMeanLength();
            os << '\t' << sqrt(GetVarLength());
        }
        os << '\t' << omegaarray->GetMean();
        os << '\t' << omegaarray->GetVar();
        os << '\t' << omegahypermean << '\t' << omegahyperinvshape;

        os << '\t' << nucstatarray->GetMeanEntropy();
        os << '\t' << nucrelratearray->GetMeanEntropy();
        if (nucmode != 2) {
            os << '\t' << sqrt(GetVarNucRelRate()) << '\t'
               << Random::GetEntropy(nucrelratehypercenter) << '\t' << nucrelratehyperinvconc;
            os << '\t' << sqrt(GetVarNucStat()) << '\t' << Random::GetEntropy(nucstathypercenter)
               << '\t' << nucstathyperinvconc;
        }
        os << '\n';
        os.flush();
    }

    // Branch lengths

    double GetMeanTotalLength() const {
        double tot = 0;
        for (int j = 0; j < Nbranch; j++) { tot += branchlength->GetVal(j); }
        return tot;
    }

    double GetMeanLength() const {
        if (blmode == 2) {
            cerr << "error: in getvarlength\n";
            exit(1);
        }

        return branchlengtharray->GetMeanLength();
    }

    double GetVarLength() const {
        if (blmode == 2) {
            cerr << "error: in getvarlength\n";
            exit(1);
        }

        return branchlengtharray->GetVarLength();
    }

    // Nucleotide rates

    double GetVarNucRelRate() const {
        if (nucmode == 2) {
            cerr << "error in getvarnucrelrate\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nrr; j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < mpi.GetNgene(); g++) {
                double tmp = (*nucrelratearray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= mpi.GetNgene();
            var /= mpi.GetNgene();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nrr;
        return tot;
    }

    double GetVarNucStat() const {
        if (nucmode == 2) {
            cerr << "error in getvarnucstat\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nnuc; j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < mpi.GetNgene(); g++) {
                double tmp = (*nucstatarray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= mpi.GetNgene();
            var /= mpi.GetNgene();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nnuc;
        return tot;
    }

    void Monitor(ostream &os) const {}

    void FromStream(istream &is) {
        if (blmode == 2) {
            is >> lambda;
            is >> *branchlength;
        } else {
            is >> lambda;
            is >> *branchlength;
            is >> blhyperinvshape;
            is >> *branchlengtharray;
        }

        is >> nucrelratehypercenter;
        is >> nucrelratehyperinvconc;
        is >> nucstathypercenter;
        is >> nucstathyperinvconc;
        is >> *nucrelratearray;
        is >> *nucstatarray;

        is >> omegahypermean;
        is >> omegahyperinvshape;
        is >> *omegaarray;
    }

    void ToStream(ostream &os) const {
        if (blmode == 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        } else {
            os << lambda << '\t';
            os << *branchlength << '\t';
            os << blhyperinvshape << '\t';
            os << *branchlengtharray << '\t';
        }

        os << nucrelratehypercenter << '\t';
        os << nucrelratehyperinvconc << '\t';
        os << nucstathypercenter << '\t';
        os << nucstathyperinvconc << '\t';
        os << *nucrelratearray << '\t';
        os << *nucstatarray << '\t';

        os << omegahypermean << '\t';
        os << omegahyperinvshape << '\t';
        os << *omegaarray << '\n';
    }

    void TraceOmega(ostream &os) const {
        for (int gene = 0; gene < mpi.GetNgene(); gene++) { os << omegaarray->GetVal(gene) << '\t'; }
        os << '\n';
        os.flush();
    }

    //-------------------
    // Updates
    //-------------------

    void UpdateNucMatrix() {
        nucmatrix->CopyStationary((*nucstatarray)[0]);
        nucmatrix->CorruptMatrix();
    }

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
        // gene contributions
        double total = GeneLogPrior;

        // branch lengths
        if (blmode == 2) {
            total += GlobalBranchLengthsLogPrior();
        } else if (blmode == 1) {
            total += GeneBranchLengthsHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        // nuc rates
        if (nucmode == 2) {
            total += GlobalNucRatesLogPrior();
        } else if (nucmode == 1) {
            total += GeneNucRatesHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        if (omegamode == 1) { total += OmegaHyperLogPrior(); }
        // already accounted for in GeneLogPrior
        // total += OmegaLogPrior();
        return total;
    }

    // Branch lengths

    double LambdaHyperLogPrior() const { return -lambda / 10; }

    double GlobalBranchLengthsLogPrior() const {
        return LambdaHyperLogPrior() + branchlength->GetLogProb();
    }

    // exponential of mean 1 for blhyperinvshape
    double BranchLengthsHyperInvShapeLogPrior() const { return -blhyperinvshape; }

    double GeneBranchLengthsHyperLogPrior() const {
        return BranchLengthsHyperInvShapeLogPrior() + branchlength->GetLogProb();
    }

    // Nucleotide rates

    double GlobalNucRatesLogPrior() const {
        return nucrelratearray->GetLogProb() + nucstatarray->GetLogProb();
    }

    // exponential of mean 1 for nucrelrate and nucstat hyper inverse
    // concentration
    double GeneNucRatesHyperLogPrior() const {
        double total = 0;
        if (nucmode == 1) {
            total -= nucrelratehyperinvconc;
            total -= nucstathyperinvconc;
        }
        return total;
    }

    // Omega

    double OmegaHyperLogPrior() const {
        double total = 0;
        total -= omegahypermean;
        total -= omegahyperinvshape;
        return total;
    }

    double OmegaLogPrior() const { return omegaarray->GetLogProb(); }

    double GetLogLikelihood() const { return lnL; }

    //-------------------
    // Suff Stat Log Probs
    //-------------------

    // Branch lengths

    // suff stat for global branch lengths, as a function of lambda
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    // suff stat for gene-specific branch lengths, as a function of bl
    // hyperparameters
    double BranchLengthsHyperSuffStatLogProb() const {
        return lengthhypersuffstatarray->GetLogProb(*branchlength, blhyperinvshape);
    }

    // Nucleotide rates

    // suff stat for global nuc rates, as a function of nucleotide matrix
    // (which itself depends on nucstat and nucrelrate)
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    // suff stat for gene-specific nuc rates, as a function of nucrate
    // hyperparameters
    double NucRatesHyperSuffStatLogProb() const {
        double total = 0;
        total += nucrelratesuffstat.GetLogProb(nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total += nucstatsuffstat.GetLogProb(nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    // Omega

    // suff stats for moving omega hyper parameters
    double OmegaHyperSuffStatLogProb() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return omegahypersuffstat.GetLogProb(alpha, beta);
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // Branch lengths

    // logprob for moving lambda
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    // logprob for moving hyperparameters of gene-specific branchlengths
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperInvShapeLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // Nucleotide rates

    // log prob for moving nuc rates hyper params
    double NucRatesHyperLogProb() const {
        return GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates
    double NucRatesLogProb() const { return GlobalNucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    // Omega

    // log prob for moving omega hyperparameters
    double OmegaHyperLogProb() const { return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb(); }

    //-------------------
    // Moves
    //-------------------

    // all methods starting with Master are called only by master
    // for each such method, there is a corresponding method called by slave, and
    // starting with Slave
    //
    // all methods starting with Gene are called only be slaves, and do some work
    // across all genes allocated to that slave

    double Move() override {
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            if (omegamode == 1) {
                ReceiveOmega();
                MoveOmegaHyperParameters();
                SendOmegaHyperParameters();
            }

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                ReceiveBranchLengthsSuffStat();
                ResampleBranchLengths();
                MoveLambda();
                SendGlobalBranchLengths();
            } else if (blmode == 1) {
                ReceiveBranchLengthsHyperSuffStat();
                MoveBranchLengthsHyperParameters();
                SendBranchLengthsHyperParameters();
            }

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                ReceiveNucPathSuffStat();
                MoveNucRates();
                SendGlobalNucRates();
            } else if (nucmode == 1) {
                ReceiveNucRatesHyperSuffStat();
                MoveNucRatesHyperParameters();
                SendNucRatesHyperParameters();
            }
        }

        // collect current state
        if (blmode != 2) { ReceiveGeneBranchLengths(); }
        if (nucmode != 2) { ReceiveGeneNucRates(); }
        ReceiveOmega();
        ReceiveLogProbs();
        return 1;
    }

    // Branch lengths

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneSingleOmegaModelMaster::LambdaHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneSingleOmegaModelMaster::LambdaHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);

        ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneSingleOmegaModelMaster::BranchLengthsHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneSingleOmegaModelMaster::BranchLengthsHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);

        branchlengtharray->SetShape(1.0 / blhyperinvshape);
        MoveLambda();
    }

    double BranchLengthsHyperScalingMove(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int j = 0; j < Nbranch; j++) {
                double deltalogprob =
                    -branchlength->GetLogProb(j) -
                    lengthhypersuffstatarray->GetVal(j).GetLogProb(
                        1.0 / blhyperinvshape, 1.0 / blhyperinvshape / branchlength->GetVal(j));
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*branchlength)[j] *= e;
                deltalogprob +=
                    branchlength->GetLogProb(j) +
                    lengthhypersuffstatarray->GetVal(j).GetLogProb(
                        1.0 / blhyperinvshape, 1.0 / blhyperinvshape / branchlength->GetVal(j));
                deltalogprob += m;
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    (*branchlength)[j] /= e;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    // Nucleotide rates

    void MoveNucRatesHyperParameters() {
        ProfileMove(nucrelratehypercenter, 1.0, 1, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
            &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelMaster::NoUpdate,
            this);
        ScalingMove(nucstathyperinvconc, 1.0, 10, &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10, &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10, &MultiGeneSingleOmegaModelMaster::NucRatesHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveNucRates() {
        vector<double> &nucrelrate = (*nucrelratearray)[0];
        ProfileMove(nucrelrate, 0.1, 1, 10, &MultiGeneSingleOmegaModelMaster::NucRatesLogProb,
            &MultiGeneSingleOmegaModelMaster::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &MultiGeneSingleOmegaModelMaster::NucRatesLogProb,
            &MultiGeneSingleOmegaModelMaster::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &MultiGeneSingleOmegaModelMaster::NucRatesLogProb,
            &MultiGeneSingleOmegaModelMaster::UpdateNucMatrix, this);

        vector<double> &nucstat = (*nucstatarray)[0];
        ProfileMove(nucstat, 0.1, 1, 10, &MultiGeneSingleOmegaModelMaster::NucRatesLogProb,
            &MultiGeneSingleOmegaModelMaster::UpdateNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 10, &MultiGeneSingleOmegaModelMaster::NucRatesLogProb,
            &MultiGeneSingleOmegaModelMaster::UpdateNucMatrix, this);
    }

    // Omega

    void MoveOmegaHyperParameters() {
        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegaarray);

        ScalingMove(omegahypermean, 1.0, 10, &MultiGeneSingleOmegaModelMaster::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        ScalingMove(omegahypermean, 0.3, 10, &MultiGeneSingleOmegaModelMaster::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        ScalingMove(omegahyperinvshape, 1.0, 10, &MultiGeneSingleOmegaModelMaster::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);
        ScalingMove(omegahyperinvshape, 0.3, 10, &MultiGeneSingleOmegaModelMaster::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelMaster::NoUpdate, this);

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    //-------------------
    // MPI send / receive
    //-------------------

    // Branch lengths

    void SendGlobalBranchLengths() { mpi.MasterSendGlobal(*branchlength); }

    void SendBranchLengthsHyperParameters() {
        mpi.MasterSendGlobal(*branchlength, blhyperinvshape);
    }

    void SendGeneBranchLengths() { mpi.MasterSendGeneArray(*branchlengtharray); }

    void ReceiveGeneBranchLengths() { mpi.MasterReceiveGeneArray(*branchlengtharray); }

    void ReceiveBranchLengthsSuffStat() {
        lengthpathsuffstatarray->Clear();
        mpi.MasterReceiveAdditive(*lengthpathsuffstatarray);
    }

    void ReceiveBranchLengthsHyperSuffStat() {
        lengthhypersuffstatarray->Clear();
        mpi.MasterReceiveAdditive(*lengthhypersuffstatarray);
    }

    // Nucleotide Rates

    void SendGlobalNucRates() {
        mpi.MasterSendGlobal(nucrelratearray->GetVal(0), nucstatarray->GetVal(0));
    }

    void SendGeneNucRates() { mpi.MasterSendGeneArray(*nucrelratearray, *nucstatarray); }

    void ReceiveGeneNucRates() { mpi.MasterReceiveGeneArray(*nucrelratearray, *nucstatarray); }

    void SendNucRatesHyperParameters() {
        mpi.MasterSendGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        mpi.MasterSendGlobal(nucstathypercenter, nucstathyperinvconc);
    }

    void ReceiveNucRatesHyperSuffStat() {
        nucrelratesuffstat.Clear();
        mpi.MasterReceiveAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        mpi.MasterReceiveAdditive(nucstatsuffstat);
    }

    void ReceiveNucPathSuffStat() {
        nucpathsuffstat.Clear();
        mpi.MasterReceiveAdditive(nucpathsuffstat);
    }

    // omega (and hyperparameters)

    void ReceiveOmega() { mpi.MasterReceiveGeneArray(*omegaarray); }

    void SendOmega() { mpi.MasterSendGeneArray(*omegaarray); }

    // omega hyperparameters

    void SendOmegaHyperParameters() { mpi.MasterSendGlobal(omegahypermean, omegahyperinvshape); }

    // log probs

    void ReceiveLogProbs() {
        GeneLogPrior = 0;
        mpi.MasterReceiveAdditive(GeneLogPrior);
        lnL = 0;
        mpi.MasterReceiveAdditive(lnL);
    }
};


class MultiGeneSingleOmegaModelSlave : public ChainComponent,
                                       public ProbModel,
                                       public MultiGeneSingleOmegaModelShared {
  private:
    MultiGeneMPIModule mpi;
    // each gene defines its own SingleOmegaModel
    std::vector<SingleOmegaModel *> geneprocess;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneSingleOmegaModelSlave(string datafile, string intreefile, int inmyid, int innprocs)
        : MultiGeneSingleOmegaModelShared(datafile, intreefile),
          mpi(inmyid, innprocs) {
        mpi.AllocateAlignments(datafile);

        refcodondata = new CodonSequenceAlignment(mpi.refdata, true);
        taxonset = mpi.refdata->GetTaxonSet();
        Ntaxa = mpi.refdata->GetNtaxa();

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        // // check whether tree and data fits together
        // tree->RegisterWith(taxonset);

        // tree->SetIndices();
        // Nbranch = tree->GetNbranch();
    }

    void start() override {}
    void move(int) override {}
    void savepoint(int) override {}
    void end() override {}

    void Allocate() {
        // Branch lengths

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
        blhyperinvshape = 0.1;
        if (blmode == 2) {
            lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
            lengthhypersuffstatarray = 0;
        } else {
            branchlength->SetAllBranches(1.0 / lambda);
            branchlengtharray = new GammaWhiteNoiseArray(
                mpi.GetLocalNgene(), *tree, *branchlength, 1.0 / blhyperinvshape);
            lengthpathsuffstatarray = 0;
            lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
        }

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;

        if (nucmode == 2) {
            nucrelratearray =
                new IIDDirichlet(1, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(1, nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = new GTRSubMatrix(Nnuc, (*nucrelratearray)[0], (*nucstatarray)[0], true);
        } else {
            nucrelratearray = new IIDDirichlet(
                mpi.GetLocalNgene(), nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray =
                new IIDDirichlet(mpi.GetLocalNgene(), nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = 0;
        }

        // Omega

        /*
        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
        */
        omegaarray = new IIDGamma(mpi.GetLocalNgene(), omegahypermean, omegahyperinvshape);

        // Gene processes

        lnL = 0;
        GeneLogPrior = 0;

        geneprocess.assign(mpi.GetLocalNgene(), (SingleOmegaModel *)0);
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene] = new SingleOmegaModel(mpi.GetLocalGeneName(gene), treefile);
            geneprocess[gene]->SetAcrossGenesModes(blmode, nucmode);
            geneprocess[gene]->Allocate();
        }
    }

    // called upon constructing the model
    // mode == 2: global
    // mode == 1: gene specific, with hyperparameters estimated across genes
    // mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode, int inomegamode) {
        blmode = inblmode;
        nucmode = innucmode;
        omegamode = inomegamode;
    }

    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape) {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    void FastUpdate() {
        branchlength->SetScale(lambda);
        if (blmode == 1) { branchlengtharray->SetShape(1.0 / blhyperinvshape); }
        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    void Update() {
        ReceiveBranchLengthsHyperParameters();
        ReceiveNucRatesHyperParameters();

        if (blmode == 2) {
            ReceiveGlobalBranchLengths();
        } else {
            ReceiveGeneBranchLengths();
        }
        if (nucmode == 2) {
            ReceiveGlobalNucRates();
        } else {
            ReceiveGeneNucRates();
        }

        ReceiveOmegaHyperParameters();
        ReceiveOmega();
        GeneUpdate();
        SendLogProbs();
    }

    void GeneUpdate() {
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) { geneprocess[gene]->Update(); }
    }

    void PostPred(string name) {
        ReceiveBranchLengthsHyperParameters();
        ReceiveNucRatesHyperParameters();

        if (blmode == 2) {
            ReceiveGlobalBranchLengths();
        } else {
            ReceiveGeneBranchLengths();
        }
        if (nucmode == 2) {
            ReceiveGlobalNucRates();
        } else {
            ReceiveGeneNucRates();
        }

        ReceiveOmegaHyperParameters();
        ReceiveOmega();
        GenePostPred(name);
    }

    void GenePostPred(string name) {
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->PostPred(name + mpi.GetLocalGeneName(gene));
        }
    }

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    const vector<double> &GetOmegaArray() const { return omegaarray->GetArray(); }

    //-------------------
    // Traces and Monitors
    //-------------------

    void TraceHeader(ostream &os) const {
        os << "#logprior\tlnL";
        if (blmode == 2) {
            os << "\tlength";
        } else {
            os << "\tmeanlength\tstdev";
        }
        os << "\tmeanomega";
        os << "\tvaromega";
        os << "\tomegahypermean\tinvshape";
        os << "\tstatent";
        os << "\trrent";
        if (nucmode != 2) {
            os << "\tstdevrr\tcenter\thyperinvconc";
            os << "\tstdevstat\tcenter\thyperinvconc";
        }
        os << '\n';
    }

    void Trace(ostream &os) const {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood();

        if (blmode == 2) {
            os << '\t' << GetMeanTotalLength();
        } else {
            os << '\t' << GetMeanLength();
            os << '\t' << sqrt(GetVarLength());
        }
        os << '\t' << omegaarray->GetMean();
        os << '\t' << omegaarray->GetVar();
        os << '\t' << omegahypermean << '\t' << omegahyperinvshape;

        os << '\t' << nucstatarray->GetMeanEntropy();
        os << '\t' << nucrelratearray->GetMeanEntropy();
        if (nucmode != 2) {
            os << '\t' << sqrt(GetVarNucRelRate()) << '\t'
               << Random::GetEntropy(nucrelratehypercenter) << '\t' << nucrelratehyperinvconc;
            os << '\t' << sqrt(GetVarNucStat()) << '\t' << Random::GetEntropy(nucstathypercenter)
               << '\t' << nucstathyperinvconc;
        }
        os << '\n';
        os.flush();
    }

    // Branch lengths

    double GetMeanTotalLength() const {
        double tot = 0;
        for (int j = 0; j < Nbranch; j++) { tot += branchlength->GetVal(j); }
        return tot;
    }

    double GetMeanLength() const {
        if (blmode == 2) {
            cerr << "error: in getvarlength\n";
            exit(1);
        }

        return branchlengtharray->GetMeanLength();
    }

    double GetVarLength() const {
        if (blmode == 2) {
            cerr << "error: in getvarlength\n";
            exit(1);
        }

        return branchlengtharray->GetVarLength();
    }

    // Nucleotide rates

    double GetVarNucRelRate() const {
        if (nucmode == 2) {
            cerr << "error in getvarnucrelrate\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nrr; j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < mpi.GetNgene(); g++) {
                double tmp = (*nucrelratearray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= mpi.GetNgene();
            var /= mpi.GetNgene();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nrr;
        return tot;
    }

    double GetVarNucStat() const {
        if (nucmode == 2) {
            cerr << "error in getvarnucstat\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nnuc; j++) {
            double mean = 0;
            double var = 0;
            for (int g = 0; g < mpi.GetNgene(); g++) {
                double tmp = (*nucstatarray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= mpi.GetNgene();
            var /= mpi.GetNgene();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nnuc;
        return tot;
    }

    void Monitor(ostream &os) const {}

    void TraceOmega(ostream &os) const {
        for (int gene = 0; gene < mpi.GetNgene(); gene++) { os << omegaarray->GetVal(gene) << '\t'; }
        os << '\n';
        os.flush();
    }

    //-------------------
    // Updates
    //-------------------

    void UpdateNucMatrix() {
        nucmatrix->CopyStationary((*nucstatarray)[0]);
        nucmatrix->CorruptMatrix();
    }

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
        // gene contributions
        double total = GeneLogPrior;

        // branch lengths
        if (blmode == 2) {
            total += GlobalBranchLengthsLogPrior();
        } else if (blmode == 1) {
            total += GeneBranchLengthsHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        // nuc rates
        if (nucmode == 2) {
            total += GlobalNucRatesLogPrior();
        } else if (nucmode == 1) {
            total += GeneNucRatesHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        if (omegamode == 1) { total += OmegaHyperLogPrior(); }
        // already accounted for in GeneLogPrior
        // total += OmegaLogPrior();
        return total;
    }

    // Branch lengths

    double LambdaHyperLogPrior() const { return -lambda / 10; }

    double GlobalBranchLengthsLogPrior() const {
        return LambdaHyperLogPrior() + branchlength->GetLogProb();
    }

    // exponential of mean 1 for blhyperinvshape
    double BranchLengthsHyperInvShapeLogPrior() const { return -blhyperinvshape; }

    double GeneBranchLengthsHyperLogPrior() const {
        return BranchLengthsHyperInvShapeLogPrior() + branchlength->GetLogProb();
    }

    // Nucleotide rates

    double GlobalNucRatesLogPrior() const {
        return nucrelratearray->GetLogProb() + nucstatarray->GetLogProb();
    }

    // exponential of mean 1 for nucrelrate and nucstat hyper inverse
    // concentration
    double GeneNucRatesHyperLogPrior() const {
        double total = 0;
        if (nucmode == 1) {
            total -= nucrelratehyperinvconc;
            total -= nucstathyperinvconc;
        }
        return total;
    }

    // Omega

    double OmegaHyperLogPrior() const {
        double total = 0;
        total -= omegahypermean;
        total -= omegahyperinvshape;
        return total;
    }

    double OmegaLogPrior() const { return omegaarray->GetLogProb(); }

    double GetLogLikelihood() const { return lnL; }

    //-------------------
    // Suff Stat Log Probs
    //-------------------

    // Branch lengths

    // suff stat for global branch lengths, as a function of lambda
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    // suff stat for gene-specific branch lengths, as a function of bl
    // hyperparameters
    double BranchLengthsHyperSuffStatLogProb() const {
        return lengthhypersuffstatarray->GetLogProb(*branchlength, blhyperinvshape);
    }

    // Nucleotide rates

    // suff stat for global nuc rates, as a function of nucleotide matrix
    // (which itself depends on nucstat and nucrelrate)
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    // suff stat for gene-specific nuc rates, as a function of nucrate
    // hyperparameters
    double NucRatesHyperSuffStatLogProb() const {
        double total = 0;
        total += nucrelratesuffstat.GetLogProb(nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total += nucstatsuffstat.GetLogProb(nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    // Omega

    // suff stats for moving omega hyper parameters
    double OmegaHyperSuffStatLogProb() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        return omegahypersuffstat.GetLogProb(alpha, beta);
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // Branch lengths

    // logprob for moving lambda
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    // logprob for moving hyperparameters of gene-specific branchlengths
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperInvShapeLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // Nucleotide rates

    // log prob for moving nuc rates hyper params
    double NucRatesHyperLogProb() const {
        return GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates
    double NucRatesLogProb() const { return GlobalNucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    // Omega

    // log prob for moving omega hyperparameters
    double OmegaHyperLogProb() const { return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb(); }

    //-------------------
    // Moves
    //-------------------


    // slave move
    double Move() override {
        GeneResampleSub(1.0);

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            MoveGeneParameters(1.0);

            if (omegamode == 1) {
                SendOmega();
                ReceiveOmegaHyperParameters();
            }

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                SendBranchLengthsSuffStat();
                ReceiveGlobalBranchLengths();
            } else if (blmode == 1) {
                SendBranchLengthsHyperSuffStat();
                ReceiveBranchLengthsHyperParameters();
            }

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                SendNucPathSuffStat();
                ReceiveGlobalNucRates();
            } else if (nucmode == 1) {
                SendNucRatesHyperSuffStat();
                ReceiveNucRatesHyperParameters();
            }
        }

        // collect current state
        if (blmode != 2) { SendGeneBranchLengths(); }
        if (nucmode != 2) { SendGeneNucRates(); }
        SendOmega();
        SendLogProbs();
        return 1;
    }

    void GeneResampleSub(double frac) {
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) { geneprocess[gene]->ResampleSub(frac); }
    }

    void MoveGeneParameters(int nrep) {
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveParameters(nrep);

            (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
            if (blmode != 2) { geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]); }
            if (nucmode != 2) {
                geneprocess[gene]->GetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
            }
        }
    }

    // Branch lengths

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void ResampleGeneBranchLengths() {
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleBranchLengths();
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
    }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneSingleOmegaModelSlave::LambdaHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneSingleOmegaModelSlave::LambdaHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);

        ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneSingleOmegaModelSlave::BranchLengthsHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneSingleOmegaModelSlave::BranchLengthsHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);

        branchlengtharray->SetShape(1.0 / blhyperinvshape);
        MoveLambda();
    }

    double BranchLengthsHyperScalingMove(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int j = 0; j < Nbranch; j++) {
                double deltalogprob =
                    -branchlength->GetLogProb(j) -
                    lengthhypersuffstatarray->GetVal(j).GetLogProb(
                        1.0 / blhyperinvshape, 1.0 / blhyperinvshape / branchlength->GetVal(j));
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*branchlength)[j] *= e;
                deltalogprob +=
                    branchlength->GetLogProb(j) +
                    lengthhypersuffstatarray->GetVal(j).GetLogProb(
                        1.0 / blhyperinvshape, 1.0 / blhyperinvshape / branchlength->GetVal(j));
                deltalogprob += m;
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    (*branchlength)[j] /= e;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    // Nucleotide rates

    void MoveNucRatesHyperParameters() {
        ProfileMove(nucrelratehypercenter, 1.0, 1, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
            &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb, &MultiGeneSingleOmegaModelSlave::NoUpdate,
            this);
        ScalingMove(nucstathyperinvconc, 1.0, 10, &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10, &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10, &MultiGeneSingleOmegaModelSlave::NucRatesHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveNucRates() {
        vector<double> &nucrelrate = (*nucrelratearray)[0];
        ProfileMove(nucrelrate, 0.1, 1, 10, &MultiGeneSingleOmegaModelSlave::NucRatesLogProb,
            &MultiGeneSingleOmegaModelSlave::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &MultiGeneSingleOmegaModelSlave::NucRatesLogProb,
            &MultiGeneSingleOmegaModelSlave::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &MultiGeneSingleOmegaModelSlave::NucRatesLogProb,
            &MultiGeneSingleOmegaModelSlave::UpdateNucMatrix, this);

        vector<double> &nucstat = (*nucstatarray)[0];
        ProfileMove(nucstat, 0.1, 1, 10, &MultiGeneSingleOmegaModelSlave::NucRatesLogProb,
            &MultiGeneSingleOmegaModelSlave::UpdateNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 10, &MultiGeneSingleOmegaModelSlave::NucRatesLogProb,
            &MultiGeneSingleOmegaModelSlave::UpdateNucMatrix, this);
    }

    // Omega

    void MoveOmegaHyperParameters() {
        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegaarray);

        ScalingMove(omegahypermean, 1.0, 10, &MultiGeneSingleOmegaModelSlave::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        ScalingMove(omegahypermean, 0.3, 10, &MultiGeneSingleOmegaModelSlave::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        ScalingMove(omegahyperinvshape, 1.0, 10, &MultiGeneSingleOmegaModelSlave::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);
        ScalingMove(omegahyperinvshape, 0.3, 10, &MultiGeneSingleOmegaModelSlave::OmegaHyperLogProb,
            &MultiGeneSingleOmegaModelSlave::NoUpdate, this);

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    //-------------------
    // MPI send / receive
    //-------------------

    // Branch lengths

    void ReceiveGlobalBranchLengths() {
        mpi.SlaveReceiveGlobal(*branchlength);
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengths(*branchlength);
        }
    }

    void ReceiveBranchLengthsHyperParameters() {
        mpi.SlaveReceiveGlobal(*branchlength, blhyperinvshape);
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
        }
    }

    void ReceiveGeneBranchLengths() {
        mpi.SlaveReceiveGeneArray(*branchlengtharray);
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengths(branchlengtharray->GetVal(gene));
        }
    }

    void SendGeneBranchLengths() {
        // in principle, redundant..
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
        mpi.SlaveSendGeneArray(*branchlengtharray);
    }

    void SendBranchLengthsSuffStat() {
        lengthpathsuffstatarray->Clear();
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectLengthSuffStat();
            lengthpathsuffstatarray->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
        }
        mpi.SlaveSendAdditive(*lengthpathsuffstatarray);
    }

    void SendBranchLengthsHyperSuffStat() {
        lengthhypersuffstatarray->Clear();
        lengthhypersuffstatarray->AddSuffStat(*branchlengtharray);
        mpi.SlaveSendAdditive(*lengthhypersuffstatarray);
    }

    // Nucleotide Rates

    void ReceiveGlobalNucRates() {
        mpi.SlaveReceiveGlobal((*nucrelratearray)[0], (*nucstatarray)[0]);

        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[0], (*nucstatarray)[0]);
        }
    }

    void ReceiveGeneNucRates() {
        mpi.SlaveReceiveGeneArray(*nucrelratearray, *nucstatarray);

        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
        }
   }

    void SendGeneNucRates() { mpi.SlaveSendGeneArray(*nucrelratearray, *nucstatarray); }

    void ReceiveNucRatesHyperParameters() {
        mpi.SlaveReceiveGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        mpi.SlaveReceiveGlobal(nucstathypercenter, nucstathyperinvconc);

        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter,
                nucrelratehyperinvconc, nucstathypercenter, nucstathyperinvconc);
        }
    }

    void SendNucRatesHyperSuffStat() {
        nucrelratesuffstat.Clear();
        nucrelratearray->AddSuffStat(nucrelratesuffstat);
        mpi.SlaveSendAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        nucstatarray->AddSuffStat(nucstatsuffstat);
        mpi.SlaveSendAdditive(nucstatsuffstat);
    }

    void SendNucPathSuffStat() {
        nucpathsuffstat.Clear();
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
        }

        mpi.SlaveSendAdditive(nucpathsuffstat);
    }

    // omega (and hyperparameters)

    void SendOmega() { mpi.SlaveSendGeneArray(*omegaarray); }

    void ReceiveOmega() {
        mpi.SlaveReceiveGeneArray(*omegaarray);
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmega((*omegaarray)[gene]);
        }
    }

    // omega hyperparameters

    void ReceiveOmegaHyperParameters() {
        mpi.SlaveReceiveGlobal(omegahypermean, omegahyperinvshape);
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape);
        }
    }

    // log probs

    void SendLogProbs() {
        GeneLogPrior = 0;
        lnL = 0;
        for (int gene = 0; gene < mpi.GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
        }
        mpi.SlaveSendAdditive(GeneLogPrior);
        mpi.SlaveSendAdditive(lnL);
    }
};
