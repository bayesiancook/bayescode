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
#include "datafile_parsers.hpp"
#include "mpi_components/Process.hpp"
#include "mpi_components/broadcast.hpp"
#include "mpi_components/gather.hpp"
#include "mpi_components/partition.hpp"
#include "mpi_components/reduce.hpp"
#include "mpi_components/scatter.hpp"
#include "operations/proxies.hpp"

struct omega_param_t {
    bool variable{false};
    double hypermean{1.0}, hyperinvshape{1.0};
};

std::istream &operator>>(std::istream &is, omega_param_t &t) {
    is >> t.variable >> t.hypermean >> t.hyperinvshape;
    return is;
}

std::ostream &operator<<(std::ostream &os, omega_param_t &t) {
    os << t.variable << '\t' << t.hypermean << '\t' << t.hyperinvshape;
    return os;
}

class MultiGeneSingleOmegaModelShared : public ChainComponent {
  public:
    MultiGeneSingleOmegaModelShared(std::string datafile, std::string treefile, param_mode_t blmode,
        param_mode_t nucmode, omega_param_t omega_param)
        : datafile(datafile),
          treefile(treefile),
          codonstatespace(Universal),
          gene_set(parse_datafile(datafile)),
          partition(gene_set, MPI::p->size - 1, 1),
          mpi(MPI::p->rank, MPI::p->size),
          blmode(blmode),
          nucmode(nucmode),
          omega_param(omega_param),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {
        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        mpi.AllocateAlignments(datafile);  // marked for deletion

        // GeneLengths gene_lengths = parse_geneset_alignments(gene_set);

        // Branch lengths

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
        blhyperinvshape = 0.1;
        if (blmode == shared) {
            lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
            lengthhypersuffstatarray = 0;
        } else {
            branchlength->SetAllBranches(1.0 / lambda);
            branchlengtharray = new GammaWhiteNoiseArray(
                partition.my_allocation_size(), *tree, *branchlength, 1.0 / blhyperinvshape);
            lengthpathsuffstatarray = 0;
            lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
        }

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;

        if (nucmode == shared) {
            nucrelratearray =
                new IIDDirichlet(1, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(1, nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = new GTRSubMatrix(Nnuc, (*nucrelratearray)[0], (*nucstatarray)[0], true);
        } else {
            nucrelratearray = new IIDDirichlet(partition.my_allocation_size(),
                nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(
                partition.my_allocation_size(), nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = 0;
        }

        // Omega
        omegaarray = new IIDGamma(
            partition.my_allocation_size(), omega_param.hypermean, omega_param.hyperinvshape);

        if (MPI::p->rank > 0) {
            geneprocess.assign(partition.my_partition_size(), (SingleOmegaModel *)0);
            for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                geneprocess[gene] =
                    new SingleOmegaModel(mpi.GetLocalGeneName(gene), treefile, blmode, nucmode);
                geneprocess[gene]->Update();
            }
        }

        declare_groups();
    }

    // MPI communication groups
    void declare_groups() {

        // branch lengths
        if (blmode == shared) {
            // clang-format off
            mpibranchlengths = make_group(
                broadcast(*branchlength),

                slave_acquire([this]() {
                    for (auto gene : geneprocess) {
                        gene->SetBranchLengths(*branchlength);
                    }
                }),

                slave_release([this]() {
                    lengthpathsuffstatarray->Clear();
                    for (auto gene : geneprocess) {
                        gene->CollectLengthSuffStat();
                        lengthpathsuffstatarray->Add(*gene->GetLengthPathSuffStatArray());
                    }
                }),

                reduce(*lengthpathsuffstatarray)
            );
            //clang-format on
        } else {
            // clang-format off
            mpibranchlengths = make_group(
                broadcast(*branchlength, blhyperinvshape),

                slave_acquire([this]() {
                    for (auto gene : geneprocess) {
                        gene->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
                    }
                }),

                slave_release([this]() {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
                    }
                    lengthhypersuffstatarray->Clear();
                    lengthhypersuffstatarray->AddSuffStat(*branchlengtharray);
                }),

                reduce(*lengthhypersuffstatarray)
            );
            //clang-format on
        }

        // nucleotide rates
        if (nucmode == shared)  {
            // clang-format off
            mpinucrates = make_group(
                broadcast(*nucrelratearray,*nucstatarray),

                slave_acquire([this]()  {
                    for (auto gene : geneprocess)   {
                        gene->SetNucRates((*nucrelratearray)[0], (*nucstatarray)[0]);
                    }
                }),

                slave_release([this]()  {
                    nucpathsuffstat.Clear();
                    for (auto gene : geneprocess)   {
                        gene->CollectNucPathSuffStat();
                        nucpathsuffstat += gene->GetNucPathSuffStat();
                    }
                }),

                reduce(nucpathsuffstat)
            );
            //clang-format on
        }
        else    {
            // clang-format off
            mpinucrates = make_group(
                broadcast(nucrelratehypercenter, nucrelratehyperinvconc,
                          nucstathypercenter, nucstathyperinvconc),

                slave_acquire([this]()  {
                    for (auto gene : geneprocess)   {
                        gene->SetNucRatesHyperParameters(nucrelratehypercenter,
                            nucrelratehyperinvconc, nucstathypercenter, nucstathyperinvconc);
                    }
                }),

                slave_release([this]()  {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->GetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
                    }
                    nucrelratesuffstat.Clear();
                    nucrelratearray->AddSuffStat(nucrelratesuffstat);
                    nucstatsuffstat.Clear();
                    nucstatarray->AddSuffStat(nucstatsuffstat);
                }),

                reduce(nucrelratesuffstat,nucstatsuffstat)
            );
            //clang-format on
        }

        // omega
        // clang-format off
        mpiomega = make_group(
                broadcast(omega_param.hypermean, omega_param.hyperinvshape),

                slave_acquire([this]()  {
                    for (auto gene : geneprocess)   {
                        gene->SetOmegaHyperParameters(
                            omega_param.hypermean, omega_param.hyperinvshape);
                    }
                }),

                slave_release([this]()  {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
                    }
                    omegahypersuffstat.Clear();
                    omegahypersuffstat.AddSuffStat(*omegaarray);
                }),

                reduce(omegahypersuffstat)
            );
            //clang-format on

        // all gene-specific variables; used when first setting up the model before starting the moves
        // clang-format off
        mpigenevariables = make_group(
                scatter(partition, *omegaarray),

                (blmode != shared) ?
                    scatter(partition, *branchlengtharray) : 
                    nullptr,


                (nucmode != shared) ? 
                    scatter(partition, *nucrelratearray, *nucstatarray) : 
                    nullptr

            );
            //clang-format on
        
        syncgenevariables = make_group(
                slave_acquire([this]()  {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->SetOmega((*omegaarray)[gene]);
                    }
                }),

                (blmode != shared) ?
                    slave_acquire([this]()  {
                        for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                            geneprocess[gene]->SetBranchLengths((*branchlengtharray)[gene]);
                        }
                    }) :
                    nullptr,

                (nucmode != shared) ? 
                    slave_acquire([this]()  {
                        for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                            geneprocess[gene]->SetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
                        }
                    }) :
                    nullptr
            );
            //clang-format on

        // clang-format off
        mpitrace = make_group(
            gather(partition, *omegaarray),

            (blmode != shared) ?
                gather(partition, *branchlengtharray) :
                nullptr,

            (nucmode != shared) ?
                gather(partition, *nucrelratearray, *nucstatarray) :
                nullptr,

            slave_release([this]()  {
                GeneLogPrior = 0;
                GeneLogLikelihood = 0;
                for (auto gene : geneprocess)   {
                    GeneLogPrior += gene->GetLogPrior();
                    GeneLogLikelihood += gene->GetLogLikelihood();
                }
            }),

            reduce(GeneLogPrior, GeneLogLikelihood)
        );
        //clang-format on
    }

    int GetNbranch() const { return tree->nb_nodes() - 1; }

    virtual ~MultiGeneSingleOmegaModelShared() = default;

    void FastUpdate() {
        branchlength->SetScale(lambda);
        if (blmode == shrunken) { branchlengtharray->SetShape(1.0 / blhyperinvshape); }
        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);

        double alpha = 1.0 / omega_param.hyperinvshape;
        double beta = alpha / omega_param.hypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
        // gene contributions
        double total = GeneLogPrior;

        // branch lengths
        if (blmode == shared) {
            total += GlobalBranchLengthsLogPrior();
        } else if (blmode == shrunken) {
            total += GeneBranchLengthsHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        // nuc rates
        if (nucmode == shared) {
            total += GlobalNucRatesLogPrior();
        } else if (nucmode == shrunken) {
            total += GeneNucRatesHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        if (omega_param.variable) { total += OmegaHyperLogPrior(); }
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
        total -= omega_param.hypermean;
        total -= omega_param.hyperinvshape;
        return total;
    }

    double OmegaLogPrior() const { return omegaarray->GetLogProb(); }

    double GetLogLikelihood() const { return GeneLogLikelihood; }

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
        return nucpathsuffstat.GetLogProb(*nucmatrix, codonstatespace);
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
        double alpha = 1.0 / omega_param.hyperinvshape;
        double beta = alpha / omega_param.hypermean;
        return omegahypersuffstat.GetLogProb(alpha, beta);
    }

    // Branch lengths

    double GetMeanTotalLength() const {
        double tot = 0;
        for (int j = 0; j < GetNbranch(); j++) { tot += branchlength->GetVal(j); }
        return tot;
    }

    double GetMeanLength() const {
        if (blmode == shared) {
            std::cerr << "error: in getvarlength\n";
            exit(1);
        }

        return branchlengtharray->GetMeanLength();
    }

    double GetVarLength() const {
        if (blmode == shared) {
            std::cerr << "error: in getvarlength\n";
            exit(1);
        }

        return branchlengtharray->GetVarLength();
    }

    // Nucleotide rates

    double GetVarNucRelRate() const {
        if (nucmode == shared) {
            std::cerr << "error in getvarnucrelrate\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nrr; j++) {
            double mean = 0;
            double var = 0;
            for (size_t g = 0; g < gene_set.size(); g++) {
                double tmp = (*nucrelratearray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= gene_set.size();
            var /= gene_set.size();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nrr;
        return tot;
    }

    double GetVarNucStat() const {
        if (nucmode == shared) {
            std::cerr << "error in getvarnucstat\n";
            exit(1);
        }

        double tot = 0;
        for (int j = 0; j < Nnuc; j++) {
            double mean = 0;
            double var = 0;
            for (size_t g = 0; g < gene_set.size(); g++) {
                double tmp = (*nucstatarray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= gene_set.size();
            var /= gene_set.size();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nnuc;
        return tot;
    }

    template <class C>
    void declare_model(C &t) {
        if (blmode == shared) {
            t.add("lambda", lambda);
            t.add("branchlength", *branchlength);
        } else {
            t.add("lambda", lambda);
            t.add("branchlength", *branchlength);
            t.add("blhyperinvshape", blhyperinvshape);
            t.add("branchlengtharray", *branchlengtharray);
        }

        t.add("nucrelratehypercenter", nucrelratehypercenter);
        t.add("nucstathyperinvconc", nucrelratehyperinvconc);
        t.add("nucstathypercenter", nucstathypercenter);
        t.add("nucstathyperinvconc", nucstathyperinvconc);
        t.add("nucrelratearray", *nucrelratearray);
        t.add("nucstatarray", *nucstatarray);

        t.add("omegahypermean", omega_param.hypermean);
        t.add("omegahyperinvshape", omega_param.hyperinvshape);
        t.add("omegaarray", *omegaarray, partition);

        t.add("genelogprior", GeneLogPrior);
        t.add("geneloglikelihood", GeneLogLikelihood);
    }

    template <class C>
    void declare_stats(C &t) {
        t.add("logprior", this, &MultiGeneSingleOmegaModelShared::GetLogPrior);
        t.add("GeneLogLikelihood", this, &MultiGeneSingleOmegaModelShared::GetLogLikelihood);

        if (blmode == shared) {
            t.add("length", this, &MultiGeneSingleOmegaModelShared::GetMeanTotalLength);
        } else {
            t.add("mean_length", this, &MultiGeneSingleOmegaModelShared::GetMeanLength);
            t.add("sd_length", [this]() { return sqrt(GetVarLength()); });
        }
        t.add("omegaarray_mean", [this]() { return omegaarray->GetMean(); });
        t.add("omegaarray_var", [this]() { return omegaarray->GetVar(); });
        t.add("omegahypermean", omega_param.hypermean);
        t.add("omegahyperinvshape", omega_param.hyperinvshape);

        t.add("nucstatarray_mean_entropy", [this]() { return nucstatarray->GetMeanEntropy(); });
        t.add(
            "nucrelratearray_mean_entropy", [this]() { return nucrelratearray->GetMeanEntropy(); });
        if (nucmode != shared) {
            t.add("sqrt_varnucrelrate", [this]() { return sqrt(GetVarNucRelRate()); });
            t.add("nucrelratehypercenter_entropy",
                [this]() { return Random::GetEntropy(nucrelratehypercenter); });
            t.add("nucrelratehyperinvconc", nucrelratehyperinvconc);
            t.add("sd_nucstat", [this]() { return sqrt(GetVarNucStat()); });
            t.add("nucstathypercenter_entropy",
                [this]() { return Random::GetEntropy(nucstathypercenter); });
            t.add("nucstathyperinvconc", nucstathyperinvconc);
        }
    }

    // ============================================================================================
    // ============================================================================================
    //   MASTER THINGS
    // ============================================================================================
    // ============================================================================================

    using M = MultiGeneSingleOmegaModelShared;

    friend std::ostream &operator<<(std::ostream &os, MultiGeneSingleOmegaModelShared &m);

    //-------------------
    // Construction and allocation
    //-------------------

    void start() override {}

    void move(int) override {
        MPI::p->message("Moving!");
        if (!MPI::p->rank){
            SendRunningStatus(1);
            MoveMaster();
        } else {
            MoveSlave();
        }
    }

    void savepoint(int) override {}

    void end() override {
        if (!MPI::p->rank) {
            SendRunningStatus(0);
        }
    }

    void SendRunningStatus(int status) { MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD); }

    void Update() { // TEMPORARY
        if (!MPI::p->rank) {
            UpdateMaster();
        } else {
            UpdateSlave();
        }
    }

    void UpdateSlave() {
        mpibranchlengths->acquire();
        mpinucrates->acquire();
        mpiomega->acquire();
        mpigenevariables->acquire();
        syncgenevariables->acquire();
        for (auto gene : geneprocess)   {
            gene->Update();
        }
        mpitrace->release();
    }

    void UpdateMaster() {
        FastUpdate();
        mpibranchlengths->release();
        mpinucrates->release();
        mpiomega->release();
        mpigenevariables->release();
        // syncgenevariables->release();
        mpitrace->acquire();
    }

    void PostPredSlave(std::string name) {
        mpibranchlengths->acquire();
        mpinucrates->acquire();
        mpiomega->acquire();
        mpigenevariables->acquire();
        syncgenevariables->acquire();
        for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
            geneprocess[gene]->PostPred(name + mpi.GetLocalGeneName(gene));
        }
    }

    void PostPredMaster(std::string name) {
        FastUpdate();
        mpibranchlengths->release();
        mpinucrates->release();
        mpiomega->release();
        mpigenevariables->release();
        // syncgenevariables->release();
    }

    //-------------------
    // Moves
    //-------------------

    // slave move
    double MoveSlave() {
        for (auto gene : geneprocess)   {
            gene->ResampleSub(1.0);
        }

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            for (auto gene : geneprocess)   {
                gene->MoveParameters(1);
            }

            if (omega_param.variable) {
                mpiomega->release();
                mpiomega->acquire();
            }

            mpibranchlengths->release();
            mpibranchlengths->acquire();

            mpinucrates->release();
            mpinucrates->acquire();
        }

        mpitrace->release();
        return 1;
    }

    double MoveMaster() {
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            if (omega_param.variable) {
                mpiomega->acquire();
                MoveOmegaHyperParameters();
                mpiomega->release();
            }

            mpibranchlengths->acquire();
            if (blmode == shared) {
                ResampleBranchLengths();
                MoveLambda();
            } else if (blmode == shrunken) {
                MoveBranchLengthsHyperParameters();
            }
            mpibranchlengths->release();

            // global nucrates, or gene nucrates hyperparameters
            mpinucrates->acquire();
            if (nucmode == shared) {
                MoveNucRates();
            } else if (nucmode == shrunken) {
                MoveNucRatesHyperParameters();
            }
            mpinucrates->release();
        }
        mpitrace->acquire();
        return 1;
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

    // Branch lengths
    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        Move::Scaling(lambda, 1.0, 10, &M::LambdaHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(lambda, 0.3, 10, &M::LambdaHyperLogProb, &M::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);

        Move::Scaling(blhyperinvshape, 1.0, 10, &M::BranchLengthsHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(blhyperinvshape, 0.3, 10, &M::BranchLengthsHyperLogProb, &M::NoUpdate, this);

        branchlengtharray->SetShape(1.0 / blhyperinvshape);
        MoveLambda();
    }

    double BranchLengthsHyperScalingMove(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int j = 0; j < GetNbranch(); j++) {
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
        Move::Profile(
            nucrelratehypercenter, 1.0, 1, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(
            nucrelratehypercenter, 0.3, 1, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(
            nucrelratehypercenter, 0.1, 3, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucrelratehyperinvconc, 1.0, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucrelratehyperinvconc, 0.3, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucrelratehyperinvconc, 0.03, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);

        Move::Profile(nucstathypercenter, 1.0, 1, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(nucstathypercenter, 0.3, 1, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(nucstathypercenter, 0.1, 2, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(nucstathyperinvconc, 1.0, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(nucstathyperinvconc, 0.3, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(nucstathyperinvconc, 0.03, 10, &M::NucRatesHyperLogProb, &M::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveNucRates() {
        std::vector<double> &nucrelrate = (*nucrelratearray)[0];
        Move::Profile(nucrelrate, 0.1, 1, 10, &M::NucRatesLogProb, &M::UpdateNucMatrix, this);
        Move::Profile(nucrelrate, 0.03, 3, 10, &M::NucRatesLogProb, &M::UpdateNucMatrix, this);
        Move::Profile(nucrelrate, 0.01, 3, 10, &M::NucRatesLogProb, &M::UpdateNucMatrix, this);

        std::vector<double> &nucstat = (*nucstatarray)[0];
        Move::Profile(nucstat, 0.1, 1, 10, &M::NucRatesLogProb, &M::UpdateNucMatrix, this);
        Move::Profile(nucstat, 0.01, 1, 10, &M::NucRatesLogProb, &M::UpdateNucMatrix, this);
    }

    // Omega

    void MoveOmegaHyperParameters() {
        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegaarray);

        Move::Scaling(omega_param.hypermean, 1.0, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(omega_param.hypermean, 0.3, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            omega_param.hyperinvshape, 1.0, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            omega_param.hyperinvshape, 0.3, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);

        double alpha = 1.0 / omega_param.hyperinvshape;
        double beta = alpha / omega_param.hypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void ToStream(std::ostream &os) { os << *this; }

    // ============================================================================================
    // ============================================================================================
    //   END OF THINGS
    // ============================================================================================
    // ============================================================================================

  protected:
    // each gene defines its own SingleOmegaModel
    std::vector<SingleOmegaModel *> geneprocess;

    std::string datafile, treefile;
    CodonStateSpace codonstatespace;
    GeneSet gene_set;
    Partition partition;
    MultiGeneMPIModule mpi;
    std::unique_ptr<const Tree> tree;

    param_mode_t blmode;
    param_mode_t nucmode;
    omega_param_t omega_param;

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
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    IIDDirichlet *nucrelratearray;
    DirichletSuffStat nucrelratesuffstat;

    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    IIDDirichlet *nucstatarray;
    DirichletSuffStat nucstatsuffstat;

    // Omega

    // each gene has its own omega
    // omegaarray[gene], for gene=1..Ngene
    // iid gamma, with hyperparameters omegahypermean and hyperinvshape
    IIDGamma *omegaarray;

    // suffstat for gene-specific omega's
    // as a function of omegahypermean and omegahyperinvshape
    GammaSuffStat omegahypersuffstat;

    // total log likelihood (summed across all genes)
    double GeneLogLikelihood;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

    std::unique_ptr<Proxy> mpiomega, mpibranchlengths, mpinucrates, mpigenevariables, syncgenevariables, mpitrace;
};

template <class M>
std::istream &operator>>(std::istream &is, std::unique_ptr<M> &m) {
    std::string model_name, datafile, treefile;
    int blmode, nucmode;
    omega_param_t omega_param;

    is >> model_name;
    if (model_name != "MultiGeneSingleOmega") {
        std::cerr << "Expected MultiGeneSingleOmega for model name, got " << model_name << "\n";
        exit(1);
    }
    is >> datafile;
    is >> treefile;
    is >> blmode >> nucmode;
    is >> omega_param;
    m.reset(new M(datafile, treefile, param_mode_t(blmode), param_mode_t(nucmode), omega_param));
    Tracer tracer{*m, &M::declare_model};
    tracer.read_line(is);
    return is;
}

std::ostream &operator<<(std::ostream &os, MultiGeneSingleOmegaModelShared &m) {
    Tracer tracer{m, &MultiGeneSingleOmegaModelShared::declare_model};
    os << "MultiGeneSingleOmega"
       << "\t";
    os << m.datafile << '\t' << m.treefile << '\t';
    os << m.blmode << '\t' << m.nucmode << '\t';
    os << m.omega_param << '\t';
    tracer.write_line(os);
    return os;
}
