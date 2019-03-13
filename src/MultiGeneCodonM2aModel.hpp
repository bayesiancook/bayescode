
/**
 * \brief A multi-gene version of CodonM2aModel
 *
 * Branch lengths and nucrates can be either:
 * - (2): shared across all genes
 * - (1): gene specific, with hyperparameters estimated across genes (shrinkage)
 * - (0): gene-specific, with fixed hyperparameters (no shrinkage); in that
 * case, the hyperparameters are set up so as to implement vague priors
 *
 * These three alternative modes can be tuned separately for branch lengths and
 * nuc rates by setting the variables blmode and nucmode to 0, 1 or 2. By
 * default, bl and nucrates are shared across genes.
 *
 * For each gene, the 3-component mixture of omega's across sites has four
 * parameters (see CodonM2aModel.hpp):
 * - 0 < purom < 1
 * - 0 < dposom < +infty
 * - 0 < purw < 1
 * - 0 <= posw < 1
 *
 * These 4 parameters are always gene-specific (they cannot be shared by all
 * genes); the priors over these parameters are as follows:
 * - purom ~ Beta(puromhypermean,puromhyperinvconc)
 * - dposom ~ Gamma(dposomhypermean,dposomhyperinvshape)
 * - purw ~ Beta purwhypermean,purwhyperinvconc)
 * - posw ~  mixture of point mass at 0 with prob 1-pi, and
 * Beta(poswhypermean,poswhyperinvconc) with prob pi
 *
 * In turn, the 9 hyperparameters of these priors can be either:
 * - (1) : estimated across genes
 * - (0): fixed (such as given by the user)
 *
 * again, separately for each set of hyperparameters, by setting the following
 * variables to 0 or 1: purommode, dposommode, purwmode, poswmode. By default,
 * the 9 mixture hyperparameters are estimated (shrunken) across genes (mode 1).
 * The hyperpriors over these parameters are then as follows:
 * - Uniform(0,1) for puromhypermean, purwhypermean and poswhypermean
 * - Exponential(1) for puromhyperinvconc, dposomhypermean, dposomhyperinvshape,
 * purwhyperinvconc, poswhyperinvconc
 * - Beta(pihypermean=0.1,pihyperinvconc=0.2) for pi
 *
 * Note that, when poswmode == 0, only poswhypermean and poswhyperinvconc are
 * fixed, on the other hand, pi is free. pi can be fixed by setting
 * pihyperinvconc to 0 (in which case pi is set equal to pihypermean).
 */

#include "Chrono.hpp"
#include "CodonM2aModel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "IIDDirichlet.hpp"

#include "components/ChainComponent.hpp"
#include "datafile_parsers.hpp"
#include "mpi_components/Process.hpp"
#include "mpi_components/broadcast.hpp"
#include "mpi_components/gather.hpp"
#include "mpi_components/partition.hpp"
#include "mpi_components/reduce.hpp"
#include "mpi_components/scatter.hpp"
#include "operations/proxies.hpp"

class MultiGeneCodonM2aModel : public ChainComponent {
  public:
    //-------------------
    // Constructors
    // ------------------

    MultiGeneCodonM2aModel(
        std::string datafile, std::string treefile, param_mode_t blmode, param_mode_t nucmode)
        : datafile(datafile),
          treefile(treefile),
          codonstatespace(Universal),
          gene_vector(parse_datafile(datafile)),
          partition(gene_vector, MPI::p->size - 1, 1),
          blmode(blmode),
          nucmode(nucmode),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {
        modalprior = 1;
        pihypermean = 0.1;
        pihyperinvconc = 0.2;

        puromhypermean = 0.5;
        puromhyperinvconc = 0.5;
        purommode = 1;

        purwhypermean = 0.5;
        purwhyperinvconc = 0.5;
        purwmode = 1;

        poswhypermean = 0.5;
        poswhyperinvconc = 0.1;
        poswmode = 1;

        dposomhypermean = 1.0;
        dposomhyperinvshape = 0.5;
        dposommode = 1;

        writegenedata = 1;

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);

        Allocate();
    }

    void Allocate() {
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

        puromarray = new IIDBeta(partition.my_allocation_size(), puromhypermean, puromhyperinvconc);
        dposomarray =
            new IIDGamma(partition.my_allocation_size(), dposomhypermean, dposomhyperinvshape);
        purwarray = new IIDBeta(partition.my_allocation_size(), purwhypermean, purwhyperinvconc);
        poswarray = new IIDBernoulliBeta(
            partition.my_allocation_size(), pi, poswhypermean, poswhyperinvconc);

        if (MPI::p->rank > 0) {
            size_t nb_genes = partition.my_partition_size();
            geneprocess.reserve(nb_genes);
            for (size_t gene_i = 0; gene_i < nb_genes; gene_i++) {
                geneprocess.emplace_back(new CodonM2aModel(gene_vector.at(gene_i), treefile, 0.1));
                // new SingleOmegaModel(gene_vector.at(gene_i), treefile, blmode, nucmode));
                geneprocess.back()->Update();
            }
        }
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

        globalomega = make_group();
        /*
        // mpi comm, gene syncing, suffstat collection and reduction over genes when moving omega
        hyper params
        // used in update and move
        globalomega = make_group(

                // forward mode (from master to slaves)

                // broadcase omega hyperparams
                broadcast(omega_param.hypermean, omega_param.hyperinvshape),

                // sync gene processes
                slave_acquire([this]()  {
                    for (auto& gene : geneprocess)   {
                        gene->SetOmegaHyperParameters(
                            omega_param.hypermean, omega_param.hyperinvshape);
                    }
                }),

                // backward mode (from slaves to master)

                // collect gene-specific omegas and suffstats
                slave_release([this]()  {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
                    }
                    omegahypersuffstat.Clear();
                    omegahypersuffstat.AddSuffStat(*omegaarray);
                }),

                // reduce suffstats on master
                reduce(omegahypersuffstat)
            );
        */

        geneomega = make_group();
        /*
        // mpi comm and gene syncing for gene-specific omega array
        // used in update and trace
        geneomega = make_group(

                // forward mode (from master to slaves)

                // scatter gene omegas across genes
                scatter(partition, *omegaarray),

                // sync gene processes
                slave_acquire([this]()  {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        geneprocess[gene]->SetOmega((*omegaarray)[gene]);
                    }
                }),

                // backward mode (from slaves to master)

                // collect gene omegas
                slave_release([this]()  {
                    for (size_t gene = 0; gene < partition.my_partition_size(); gene++) {
                        (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
                    }
                }),

                // gather them onto master
                gather(partition, *omegaarray)
            );
        */

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

    virtual ~MultiGeneCodonM2aModel() = default;

    int GetNbranch() const { return tree->nb_nodes() - 1; }

    void FastUpdate() {}

    //-------------------
    // Traces and Monitors
    // ------------------

    // summary statistics for tracing MCMC
    int GetNpos() const { return poswarray->GetSize() - poswarray->GetNullSet(); }

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

    /*
    void TracePosWeight(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << poswarray->GetVal(gene) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TracePosOm(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << 1 + dposomarray->GetVal(gene) << '\t';
        }
        os << '\n';
        os.flush();
    }
    */

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

        // add mixture

        t.add("genelogprior", GeneLogPrior);
        t.add("geneloglikelihood", GeneLogLikelihood);
    }

    template <class C>
    void declare_stats(C &t) {
        t.add("logprior", this, &MultiGeneCodonM2aModel::GetLogPrior);
        t.add("GeneLogLikelihood", this, &MultiGeneCodonM2aModel::GetLogLikelihood);

        if (blmode == shared) {
            t.add("length", this, &MultiGeneCodonM2aModel::GetMeanTotalLength);
        } else {
            t.add("mean_length", this, &MultiGeneCodonM2aModel::GetMeanLength);
            t.add("sd_length", [this]() { return sqrt(GetVarLength()); });
        }

        // add summary stats for omega mixture

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

    using M = MultiGeneCodonM2aModel;

    friend std::ostream &operator<<(std::ostream &os, MultiGeneCodonM2aModel &m);

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

            // if (omega_param.variable) {
            if (true) {
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
            if (true) {
                // if (omega_param.variable) {
                // get suffstats from slaves
                globalomega->acquire();
                // move omega hyperparams
                MoveOmegaMixture();
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
    // Log Priors
    // ------------------

    // total log prior
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

        total += MixtureHyperLogPrior();

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

    // mixture
    double MixtureHyperLogPrior() const {
        double total = 0;
        if (pi) {
            // beta distribution for pi, if not 0
            double pialpha = pihypermean / pihyperinvconc;
            double pibeta = (1 - pihypermean) / pihyperinvconc;
            total += (pialpha - 1) * log(1.0 - pi) + (pibeta - 1) * log(pi);
        }
        // exponential of mean 1 for purom and purw hyperinvconc
        total -= puromhyperinvconc;
        total -= purwhyperinvconc;

        // exponential of mean 0.1 for poswhyperinvconc
        total -= 10 * poswhyperinvconc;
        // exponential of mean 1 for dposomhypermean
        total -= dposomhypermean;
        // exponential of mean 0.1 for dposomhyperinvshape
        total -= 10 * dposomhyperinvshape;

        // dposom:
        // distribution across genes should be modal
        if (modalprior && (dposomhyperinvshape > 1.0)) {
            total += log(0);
            // total += Random::INFPROB;
        }
        // distribution mean should not be too close to 0 (hypermean>0.5)
        /*
        if (dposomhypermean < 0.5)  {
            total += log(0);
            // total += Random::INFPROB;
        }
        */
        // posw:
        // distribution across genes should be modal
        double alpha = poswhypermean / poswhyperinvconc;
        double beta = (1 - poswhypermean) / poswhyperinvconc;
        if (modalprior && ((alpha < 1) || (beta < 1))) {
            total += log(0);
            // total += Random::INFPROB;
        }
        // distribution mean should not be too close to 0 (hypermean>0.1)
        /*
        if (poswhypermean < 0.1)    {
            total += log(0);
            // total += Random::INFPROB;
        }
        */
        return total;
    }

    //-------------------
    // Log Likelihood
    // ------------------

    double GetLogLikelihood() const { return GeneLogLikelihood; }

    // suff stat for gene-specific mixture parameters, as a function of mixture
    // hyperparameters
    double MixtureHyperSuffStatLogProb() const {
        double total = 0;
        total += puromsuffstat.GetLogProb(puromhypermean, puromhyperinvconc);
        total += dposomsuffstat.GetLogProb(dposomhypermean, dposomhyperinvshape);
        total += purwsuffstat.GetLogProb(purwhypermean, purwhyperinvconc);
        total += poswsuffstat.GetLogProb(pi, poswhypermean, poswhyperinvconc);
        return total;
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // Update functions for MH moves
    void NoUpdate() {}
    void VectorNoUpdate(int i) {}
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary((*nucstatarray)[0]);
        nucmatrix->CorruptMatrix();
    }

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

    // log prob for moving mixture hyper params
    double MixtureHyperLogProb() const {
        return MixtureHyperLogPrior() + MixtureHyperSuffStatLogProb();
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

    void MoveOmegaMixture() {}

    /*
    // moving mixture hyper params
    void MoveMixtureHyperParameters();
    void MovePoswHyper();
    int PoswCompMove(double tuning);
    // special function for moving pi
    void ResamplePi();
    */

    void ToStream(std::ostream &os) { os << *this; }

  protected:
    int modalprior;

    int purommode;
    int dposommode;
    int purwmode;
    int poswmode;

    int writegenedata;

    std::vector<std::unique_ptr<CodonM2aModel>> geneprocess;

    std::string datafile, treefile;
    CodonStateSpace codonstatespace;
    GeneSet gene_vector;
    Partition partition;
    std::unique_ptr<const Tree> tree;

    param_mode_t blmode;
    param_mode_t nucmode;

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
    double puromhypermean;
    double puromhyperinvconc;
    IIDBeta *puromarray;
    BetaSuffStat puromsuffstat;

    double dposomhypermean;
    double dposomhyperinvshape;
    IIDGamma *dposomarray;
    GammaSuffStat dposomsuffstat;

    double purwhypermean;
    double purwhyperinvconc;
    IIDBeta *purwarray;
    BetaSuffStat purwsuffstat;

    double poswhypermean;
    double poswhyperinvconc;
    IIDBernoulliBeta *poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double pi;

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

    is >> model_name;
    if (model_name != "MultiGeneCodonM2a") {
        std::cerr << "Expected MultiGeneCodonM2a for model name, got " << model_name << "\n";
        exit(1);
    }
    is >> datafile;
    is >> treefile;
    is >> blmode >> nucmode;
    m.reset(new M(datafile, treefile, param_mode_t(blmode), param_mode_t(nucmode)));
    Tracer tracer{*m};
    tracer.read_line(is);
    return is;
}

std::ostream &operator<<(std::ostream &os, MultiGeneCodonM2aModel &m) {
    Tracer tracer{m};
    os << "MultiGeneCodonM2a"
       << "\t";
    os << m.datafile << '\t' << m.treefile << '\t';
    os << m.blmode << '\t' << m.nucmode << '\t';
    tracer.write_line(os);
    return os;
}
