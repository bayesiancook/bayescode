// this is a multigene version of singleomegamodel
//
// - the array of gene-specific omega's are iid gamma with hyperparameters
// omegahypermean and omegahyperinvshape
//

#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "SingleOmegaModel.hpp"
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

class MultiGeneSingleOmegaModel : public ChainComponent {
  public:
    MultiGeneSingleOmegaModel(std::string datafile, std::string treefile, param_mode_t blmode,
        param_mode_t nucmode, omega_param_t omega_param)
        : datafile(datafile),
          treefile(treefile),
          codonstatespace(Universal),
          gene_vector(parse_datafile(datafile)),
          partition(gene_vector, MPI::p->size - 1, 1),
          blmode(blmode),
          nucmode(nucmode),
          omega_param(omega_param),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {
        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        // GeneLengths gene_lengths = parse_geneset_alignments(gene_vector);

        // Branch lengths

        meanoverbranches = 0.1;
        one = 1.0;

        // if shared, then these are the global branch lengths
        // otherwise, they are the bector of branch lengths hypermeans across branches
        branchlength = new BranchIIDGamma(*tree, meanoverbranches, one);

        // if branch lengths are not shared, this is the inv shape for the variation of branch
        // lengths across GENES
        blhyperinvshape = 1.0;

        if (blmode == shared) {
            branchlengtharray = nullptr;
            // an array across branches of suff stats (for substitution mappings across genes)
            lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
            lengthhypersuffstatarray = nullptr;
        } else {
            // just for ensuring reasonable initialization -- could possibly be removed
            branchlength->SetAllBranches(meanoverbranches);
            branchlengtharray = new GammaWhiteNoiseArray(
                partition.my_allocation_size(), *tree, *branchlength, blhyperinvshape);
            lengthpathsuffstatarray = nullptr;
            // an array across branches of GammaSuffStats (for gene-specific branch lengths)
            lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
        }

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;

        if (nucmode == shared) {
            nucrelratearray = new IIDDirichlet(1, nucrelratehypercenter, nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(1, nucstathypercenter, nucstathyperinvconc);
            nucmatrix = new GTRSubMatrix(Nnuc, (*nucrelratearray)[0], (*nucstatarray)[0], true);
        } else {
            nucrelratearray = new IIDDirichlet(
                partition.my_allocation_size(), nucrelratehypercenter, nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(
                partition.my_allocation_size(), nucstathypercenter, nucstathyperinvconc);
            nucmatrix = 0;
        }

        // Omega
        omegaarray = new IIDGamma(
            partition.my_allocation_size(), omega_param.hypermean, omega_param.hyperinvshape);

        if (MPI::p->rank > 0) {
            size_t nb_genes = partition.my_partition_size();
            geneprocess.reserve(nb_genes);
            for (size_t gene_i = 0; gene_i < nb_genes; gene_i++) {
                geneprocess.emplace_back(
                    new SingleOmegaModel(gene_vector.at(gene_i), treefile, blmode, nucmode));
                geneprocess.back()->Update();
            }
        }

        declare_groups();
    }

    // MPI communication groups
    void declare_groups() {
        // clang-format off

        // mpi comm, gene syncing, suffstat collection and reduction over genes when moving branch
        // lengths used in update and move
        if (blmode == shared) {
            globalbl = make_group(

                // forward mode (from master to slaves)

                // broadcast global branch lengths
                broadcast(*branchlength),

                // sync gene processes
                slave_acquire([this]() {
                    for (auto &gene : geneprocess) { gene->SetBranchLengths(*branchlength); }
                }),

                // backward mode (from slaves to master)

                // collect bl suff stats across genes
                slave_release([this]() {
                    lengthpathsuffstatarray->Clear();
                    for (auto &gene : geneprocess) {
                        gene->CollectLengthSuffStat();
                        lengthpathsuffstatarray->Add(*gene->GetLengthPathSuffStatArray());
                    }
                }),

                // reduce suffstats on master
                reduce(*lengthpathsuffstatarray));
        } else {
            globalbl = make_group(

                // forward mode (from master to slaves)

                // broadcast gene branch lengths hyperparam
                broadcast(*branchlength, blhyperinvshape),

                // sync gene processes
                slave_acquire([this]() {
                    for (auto &gene : geneprocess) {
                        gene->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
                    }
                }),

                // backward mode (from slaves to master)

                // update arrays of gene branch lengths for each slave and collect hyper suffstats
                slave_release([this]() {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
                    }
                    lengthhypersuffstatarray->Clear();
                    lengthhypersuffstatarray->AddSuffStat(*branchlengtharray);
                }),

                // reduce suffstats on master
                reduce(*lengthhypersuffstatarray));
        }

        // mpi comm, gene syncing, suffstat collection and reduction over genes when moving nuc
        // rates used in update and move
        if (nucmode == shared) {
            globalnuc = make_group(

                // forward mode (from master to slaves)

                // broadcast global nuc rates
                broadcast(*nucrelratearray, *nucstatarray),

                // sync gene processes
                slave_acquire([this]() {
                    for (auto &gene : geneprocess) {
                        gene->SetNucRates((*nucrelratearray)[0], (*nucstatarray)[0]);
                    }
                }),

                // backward mode (from slaves to master)

                // collect nuc path suffstat across genes, for each slave
                slave_release([this]() {
                    nucpathsuffstat.Clear();
                    for (auto &gene : geneprocess) {
                        gene->CollectNucPathSuffStat();
                        nucpathsuffstat += gene->GetNucPathSuffStat();
                    }
                }),

                // reduce suff stats on master
                reduce(nucpathsuffstat));
        } else {
            globalnuc = make_group(

                // forward mode (from master to slaves)

                // broadcast nucrates hyperparams
                broadcast(nucrelratehypercenter, nucrelratehyperinvconc, nucstathypercenter,
                    nucstathyperinvconc),

                // sync gene processes
                slave_acquire([this]() {
                    for (auto &gene : geneprocess) {
                        gene->SetNucRatesHyperParameters(nucrelratehypercenter,
                            nucrelratehyperinvconc, nucstathypercenter, nucstathyperinvconc);
                    }
                }),

                // backward mode (from slaves to master)

                // collect gene-specific nuc rates and hyper suff stats
                slave_release([this]() {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->GetNucRates(
                            (*nucrelratearray)[gene], (*nucstatarray)[gene]);
                    }
                    nucrelratesuffstat.Clear();
                    nucrelratearray->AddSuffStat(nucrelratesuffstat);
                    nucstatsuffstat.Clear();
                    nucstatarray->AddSuffStat(nucstatsuffstat);
                }),

                // reduce suffstats on master
                reduce(nucrelratesuffstat, nucstatsuffstat));
        }

        // mpi comm, gene syncing, suffstat collection and reduction over genes when moving omega
        // hyper params used in update and move
        globalomega = make_group(

            // forward mode (from master to slaves)

            // broadcase omega hyperparams
            broadcast(omega_param.hypermean, omega_param.hyperinvshape),

            // sync gene processes
            slave_acquire([this]() {
                for (auto &gene : geneprocess) {
                    gene->SetOmegaHyperParameters(omega_param.hypermean, omega_param.hyperinvshape);
                }
            }),

            // backward mode (from slaves to master)

            // collect gene-specific omegas and suffstats
            slave_release([this]() {
                for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                    (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
                }
                omegahypersuffstat.Clear();
                omegahypersuffstat.AddSuffStat(*omegaarray);
            }),

            // reduce suffstats on master
            reduce(omegahypersuffstat));

        // mpi comm and gene syncing for gene-specific omega array
        // used in update and trace
        geneomega = make_group(

            // forward mode (from master to slaves)

            // scatter gene omegas across genes
            scatter(partition, *omegaarray),

            // sync gene processes
            slave_acquire([this]() {
                for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                    geneprocess[gene]->SetOmega((*omegaarray)[gene]);
                }
            }),

            // backward mode (from slaves to master)

            // collect gene omegas
            slave_release([this]() {
                for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                    (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
                }
            }),

            // gather them onto master
            gather(partition, *omegaarray));

        // mpi comm and gene syncing for gene-specific branch lengths
        // used in update and trace
        if (blmode != shared) {
            genebl = make_group(

                // forward mode (from master to slaves)

                scatter(partition, *branchlengtharray),

                slave_acquire([this]() {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->SetBranchLengths((*branchlengtharray)[gene]);
                    }
                }),

                // backward mode (from slaves to master)

                slave_release([this]() {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
                    }
                }),

                gather(partition, *branchlengtharray));
        }

        // mpi comm and gene syncing for gene-specific nuc rates
        // used in update and trace
        if (nucmode != shared) {
            genenuc = make_group(

                // forward mode (from master to slaves)

                scatter(partition, *nucrelratearray, *nucstatarray),

                slave_acquire([this]() {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->SetNucRates(
                            (*nucrelratearray)[gene], (*nucstatarray)[gene]);
                    }
                }),

                // backward mode (from slaves to master)

                slave_release([this]() {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->GetNucRates(
                            (*nucrelratearray)[gene], (*nucstatarray)[gene]);
                    }
                }),

                gather(partition, *nucrelratearray, *nucstatarray));
        }

        // reduction of summary statistics across genes
        // used in trace
        genestats = make_group(

            slave_release([this]() {
                GeneLogPrior = 0;
                GeneLogLikelihood = 0;
                for (auto &gene : geneprocess) {
                    GeneLogPrior += gene->GetLogPrior();
                    GeneLogLikelihood += gene->GetLogLikelihood();
                }
            }),

            reduce(GeneLogPrior, GeneLogLikelihood));

        // clang-format on
    }

    int GetNbranch() const { return tree->nb_nodes() - 1; }

    virtual ~MultiGeneSingleOmegaModel() = default;

    // void FastUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
        // gene contributions
        double total = GeneLogPrior;

        // branch lengths
        if (blmode == shared) {
            total += BranchLengthsMeanLogPrior();
            total += BranchLengthsLogPrior();
        } else if (blmode == shrunken) {
            total += BranchLengthsMeanLogPrior();
            total += BranchLengthsLogPrior();
            total += GeneBranchLengthsHyperInvShapeLogPrior();
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

    double BranchLengthsMeanLogPrior() const { return -meanoverbranches / 0.1; }

    double GeneBranchLengthsHyperInvShapeLogPrior() const { return -blhyperinvshape; }

    double BranchLengthsLogPrior() const { return branchlength->GetLogProb(); }

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

    // Branch lengths

    double GetMeanTotalLength() const {
        double tot = 0;
        for (int j = 0; j < GetNbranch(); j++) { tot += branchlength->GetVal(j); }
        return tot;
    }

    double GetMeanLength() const {
        assert(blmode != shared);
        return branchlengtharray->GetMeanLength();
    }

    double GetVarLength() const {
        assert(blmode != shared);
        return branchlengtharray->GetVarLength();
    }

    // Nucleotide rates

    double GetVarNucRelRate() const {
        assert(nucmode != shared);
        double tot = 0;
        for (int j = 0; j < Nrr; j++) {
            double mean = 0;
            double var = 0;
            for (size_t g = 0; g < gene_vector.size(); g++) {
                double tmp = (*nucrelratearray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= gene_vector.size();
            var /= gene_vector.size();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nrr;
        return tot;
    }

    double GetVarNucStat() const {
        assert(nucmode != shared);
        double tot = 0;
        for (int j = 0; j < Nnuc; j++) {
            double mean = 0;
            double var = 0;
            for (size_t g = 0; g < gene_vector.size(); g++) {
                double tmp = (*nucstatarray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= gene_vector.size();
            var /= gene_vector.size();
            var -= mean * mean;
            tot += var;
        }
        tot /= Nnuc;
        return tot;
    }

    template <class C>
    void declare_model(C &t) {
        if (blmode == shared) {
            t.add("meanoverbranches", meanoverbranches);
            t.add("branchlength", *branchlength);
        } else {
            t.add("meanoverbranches", meanoverbranches);
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
        t.add("logprior", this, &MultiGeneSingleOmegaModel::GetLogPrior);
        t.add("GeneLogLikelihood", this, &MultiGeneSingleOmegaModel::GetLogLikelihood);

        if (blmode == shared) {
            t.add("length", this, &MultiGeneSingleOmegaModel::GetMeanTotalLength);
        } else {
            t.add("mean_length", this, &MultiGeneSingleOmegaModel::GetMeanLength);
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

    using M = MultiGeneSingleOmegaModel;

    friend std::ostream &operator<<(std::ostream &os, MultiGeneSingleOmegaModel &m);

    //-------------------
    // Construction and allocation
    //-------------------

    void start() override {}

    void move(int) override {
        MPI::p->message("Moving!");
        if (!MPI::p->rank) {
            SendRunningStatus(1);
            MoveMaster();
        } else {
            MoveSlave();
        }
    }

    void savepoint(int) override {}

    void end() override {
        if (!MPI::p->rank) { SendRunningStatus(0); }
    }

    void SendRunningStatus(int status) { MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD); }

    void Update() {  // TEMPORARY
        MPI::p->message("Updating");
        if (!MPI::p->rank) {
            UpdateMaster();
        } else {
            UpdateSlave();
        }
    }

    void UpdateSlave() {
        globalbl->acquire();
        globalnuc->acquire();
        globalomega->acquire();
        if (blmode != shared) { genebl->acquire(); }
        if (nucmode != shared) { genenuc->acquire(); }
        geneomega->acquire();

        for (auto &gene : geneprocess) { gene->Update(); }

        genestats->release();
    }

    void UpdateMaster() {
        // FastUpdate();
        globalbl->release();
        globalnuc->release();
        globalomega->release();
        if (blmode != shared) { genebl->release(); }
        if (nucmode != shared) { genenuc->release(); }
        geneomega->release();

        genestats->acquire();
    }

    void PostPredSlave(std::string name) {
        globalbl->acquire();
        globalnuc->acquire();
        globalomega->acquire();
        if (blmode != shared) { genebl->acquire(); }
        if (nucmode != shared) { genenuc->acquire(); }
        geneomega->acquire();
        for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
            geneprocess[gene]->PostPred(name + gene_vector.at(gene));
        }
    }

    void PostPredMaster(std::string name) {
        // FastUpdate();
        globalbl->release();
        globalnuc->release();
        globalomega->release();
        if (blmode != shared) { genebl->release(); }
        if (nucmode != shared) { genenuc->release(); }
        geneomega->release();
    }

    //-------------------
    // Moves
    //-------------------

    // slave move
    double MoveSlave() {
        for (auto &gene : geneprocess) {
            // resample substitution history
            gene->ResampleSub(1.0);
        }

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            for (auto &gene : geneprocess) {
                // move whatever is gene-specific
                gene->MoveParameters(1);
            }

            if (omega_param.variable) {
                // release suffstats
                globalomega->release();
                // here, master moves omega hyperparams
                // get new hyperparams (and sync genes)
                globalomega->acquire();
            }

            if (blmode != independent) {
                // release suffstats
                globalbl->release();
                // here, master moves global bl or genebl-hyperparams
                // get new global bl or gene bl-hyperparams (and sync genes)
                globalbl->acquire();
            }

            if (nucmode != independent) {
                // release suffstats
                globalnuc->release();
                // here, master moves global nucrates or genenuc-hyperparams
                // get new global nucrates or genenuc-hyperparams (and sync genes)
                globalnuc->acquire();
            }
        }

        // gather gene-specific variables on master
        geneomega->release();
        if (blmode != shared) { genebl->release(); }
        if (blmode != shared) { genenuc->release(); }

        // reduce summary statistics across genes on master
        genestats->release();
        return 1;
    }

    double MoveMaster() {
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            if (omega_param.variable) {
                // get suffstats from slaves
                globalomega->acquire();
                // move omega hyperparams
                MoveOmegaHyperParameters();
                // broadcast new hyperparams to slaves
                globalomega->release();
            }

            if (blmode != independent) {
                // get suffstats for global bl or gene-bl hyperparams
                globalbl->acquire();
                // move global variables
                if (blmode == shared) {
                    ResampleGlobalBranchLengths();
                    MoveGlobalBranchLengthsHyperMean();
                } else if (blmode == shrunken) {
                    MoveGeneBranchLengthsHyperParameters();
                }
                // broadcast new global bl or gene-bl hyperparams to slaves
                globalbl->release();
            }

            if (nucmode != independent) {
                // get suffstats for global nucrates, or gene nucrates hyperparameters
                globalnuc->acquire();
                // move global variables
                if (nucmode == shared) {
                    MoveGlobalNucRates();
                } else if (nucmode == shrunken) {
                    MoveGeneNucRatesHyperParameters();
                }
                globalnuc->release();
            }
        }

        // gather gene-specific variables on master
        geneomega->acquire();
        if (blmode != shared) { genebl->acquire(); }
        if (nucmode != shared) { genenuc->acquire(); }

        // reduce summary statistics across genes on master
        genestats->acquire();
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
    void VectorNoUpdate(int i) {}

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // Branch lengths

    // logprob for moving meanoverbranches
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsMeanLogPrior() + hyperlengthsuffstat.GetLogProb(meanoverbranches, 1.0);
    }

    // logprob for moving blhyperinvshape (hyper invshape of gene-specific branchlengths)
    double GeneBranchLengthsHyperInvShapeLogProb() const {
        return GeneBranchLengthsHyperInvShapeLogPrior() +
               lengthhypersuffstatarray->GetLogProb(*branchlength, blhyperinvshape);
    }

    // logprob for moving hyperparameters of gene-specific branchlengths (branchlength array, in bl
    // shrunken mode)
    double GeneBranchLengthsHyperMeanLogProb(int j) const {
        return branchlength->GetLogProb(j) + lengthhypersuffstatarray->GetVal(j).GetLogProb(
                                                 branchlength->GetVal(j), blhyperinvshape);
    }

    // Nucleotide rates

    // log prob for moving nuc rates
    double GlobalNucRatesLogProb() const {
        return GlobalNucRatesLogPrior() + nucpathsuffstat.GetLogProb(*nucmatrix, codonstatespace);
    }

    // log prob for moving nuc rates hyper params
    double GeneNucRatesHyperLogProb() const {
        double total = GeneNucRatesHyperLogPrior();
        total += nucrelratesuffstat.GetLogProb(nucrelratehypercenter, nucrelratehyperinvconc);
        total += nucstatsuffstat.GetLogProb(nucstathypercenter, nucstathyperinvconc);
        return total;
    }

    // Omega

    // log prob for moving omega hyperparameters
    double OmegaHyperLogProb() const {
        return OmegaHyperLogPrior() +
               omegahypersuffstat.GetLogProb(omega_param.hypermean, omega_param.hyperinvshape);
    }

    //-------------------
    // MCMC moves
    //-------------------

    // Branch lengths

    void ResampleGlobalBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void MoveGlobalBranchLengthsHyperMean() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        Move::Scaling(meanoverbranches, 1.0, 10, &M::BranchLengthsHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(meanoverbranches, 0.3, 10, &M::BranchLengthsHyperLogProb, &M::NoUpdate, this);
    }

    void MoveGeneBranchLengthsHyperParameters() {
        Move::VectorScaling(branchlength->GetArray(), 1.0, 10,
            &M::GeneBranchLengthsHyperMeanLogProb, &M::VectorNoUpdate, this);
        Move::VectorScaling(branchlength->GetArray(), 0.3, 10,
            &M::GeneBranchLengthsHyperMeanLogProb, &M::VectorNoUpdate, this);
        Move::Scaling(blhyperinvshape, 1.0, 10, &M::GeneBranchLengthsHyperInvShapeLogProb,
            &M::NoUpdate, this);
        Move::Scaling(blhyperinvshape, 0.3, 10, &M::GeneBranchLengthsHyperInvShapeLogProb,
            &M::NoUpdate, this);

        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        Move::Scaling(meanoverbranches, 1.0, 10, &M::BranchLengthsHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(meanoverbranches, 0.3, 10, &M::BranchLengthsHyperLogProb, &M::NoUpdate, this);
    }

    // Nucleotide rates

    void MoveGeneNucRatesHyperParameters() {
        Move::Profile(
            nucrelratehypercenter, 1.0, 1, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(
            nucrelratehypercenter, 0.3, 1, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(
            nucrelratehypercenter, 0.1, 3, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucrelratehyperinvconc, 1.0, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucrelratehyperinvconc, 0.3, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucrelratehyperinvconc, 0.03, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);

        Move::Profile(
            nucstathypercenter, 1.0, 1, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(
            nucstathypercenter, 0.3, 1, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Profile(
            nucstathypercenter, 0.1, 2, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucstathyperinvconc, 1.0, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucstathyperinvconc, 0.3, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            nucstathyperinvconc, 0.03, 10, &M::GeneNucRatesHyperLogProb, &M::NoUpdate, this);
    }

    void MoveGlobalNucRates() {
        std::vector<double> &nucrelrate = (*nucrelratearray)[0];
        Move::Profile(nucrelrate, 0.1, 1, 10, &M::GlobalNucRatesLogProb, &M::UpdateNucMatrix, this);
        Move::Profile(
            nucrelrate, 0.03, 3, 10, &M::GlobalNucRatesLogProb, &M::UpdateNucMatrix, this);
        Move::Profile(
            nucrelrate, 0.01, 3, 10, &M::GlobalNucRatesLogProb, &M::UpdateNucMatrix, this);

        std::vector<double> &nucstat = (*nucstatarray)[0];
        Move::Profile(nucstat, 0.1, 1, 10, &M::GlobalNucRatesLogProb, &M::UpdateNucMatrix, this);
        Move::Profile(nucstat, 0.01, 1, 10, &M::GlobalNucRatesLogProb, &M::UpdateNucMatrix, this);
    }

    // Omega

    void MoveOmegaHyperParameters() {
        Move::Scaling(omega_param.hypermean, 1.0, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(omega_param.hypermean, 0.3, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            omega_param.hyperinvshape, 1.0, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);
        Move::Scaling(
            omega_param.hyperinvshape, 0.3, 10, &M::OmegaHyperLogProb, &M::NoUpdate, this);
    }

    void ToStream(std::ostream &os) { os << *this; }

    // ============================================================================================
    // ============================================================================================
    //   END OF THINGS
    // ============================================================================================
    // ============================================================================================

  protected:
    // each gene defines its own SingleOmegaModel
    std::vector<std::unique_ptr<SingleOmegaModel>> geneprocess;

    std::string datafile, treefile;
    CodonStateSpace codonstatespace;
    GeneSet gene_vector;
    Partition partition;
    std::unique_ptr<const Tree> tree;

    param_mode_t blmode;
    param_mode_t nucmode;
    omega_param_t omega_param;

    // Branch lengths

    double meanoverbranches;
    double one;
    BranchIIDGamma *branchlength;
    double blhyperinvshape;
    GammaWhiteNoiseArray *branchlengtharray;

    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStatBranchArray *lengthhypersuffstatarray;
    GammaSuffStat hyperlengthsuffstat;

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

    std::unique_ptr<Proxy> globalomega, globalbl, globalnuc, geneomega, genebl, genenuc, genestats;
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
    Tracer tracer{*m};
    tracer.read_line(is);
    return is;
}

std::ostream &operator<<(std::ostream &os, MultiGeneSingleOmegaModel &m) {
    Tracer tracer{m};
    os << "MultiGeneSingleOmega"
       << "\t";
    os << m.datafile << '\t' << m.treefile << '\t';
    os << m.blmode << '\t' << m.nucmode << '\t';
    os << m.omega_param << '\t';
    tracer.write_line(os);
    return os;
}
