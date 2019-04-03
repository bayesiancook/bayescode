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

#pragma once

#include "AADiffSelCodonMatrixBidimArray.hpp"
#include "CodonSequenceAlignment.hpp"
#include "DiffSelSparseFitnessArray.hpp"
#include "DiffselDoubleSparseConfig.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "IIDMultiBernoulli.hpp"
#include "IIDMultiGamma.hpp"
#include "IIDProfileMask.hpp"
#include "Move.hpp"
#include "MultiGammaSuffStat.hpp"
#include "PathSuffStat.hpp"
#include "PhyloProcess.hpp"
#include "SubMatrixSelector.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"
#include "components/mh_utils.hpp"
#include "components/probnode_utils.hpp"
#include "tree/implem.hpp"

/**
 * \brief A doubly-sparse version of the differential selection model (see also
 * DiffSelModel and DiffSelSparseModel)
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
 * conditions. The model is based on the following system of random variables:
 * - BidimIIDMultiGamma *fitness (G_kia): an array of site- and condition
 * specific pre-fitness parameters, for condition k=0..K, site i and amino-acid
 * a.
 * - IIDProfileMask *sitemaskarray (m_ia): an array of masks (for amino-acid a
 * and site i)
 * - BidimIIDMultiBernoulli *toggle (d_kia): an array of site and
 * condition-specific toggles (only for alternative conditions, k=1..K).
 *
 * As in AAMutSelSparseModel, the site-specific fitness masks specify which
 * amino-acids have a high (1) or a low (0) fitness. High-fitness amino-acids
 * take their fitness value from the G_kia vector; for low-fitness amino-acids,
 * the fitness is equal to some background level (maskepsilon) -- and this,
 * across all conditions. On the top of this masking system (which infuences all
 * conditions uniformly), the toggles determine whether the fitness for each
 * amino-acid and at each site should change, upon going from condition 0 to
 * condition k. Quantitatively, the fitness vector at site i under condition k,
 * (DiffSelDoublySparseFitnessArray *fitnessprofile in the code, mathematically
 * denoted F_kia in the following) is defined as follows:
 * - F_0ia = G_0ia * m_ia + epsilon * (1-m_ia)
 * - F_kia = F_0ia^(1-d_kia) * G_kia^(d_kia), for k=1..K-1
 *
 * Statistical support for a differential effect between conditions, for a given
 * site i and a given amino acid a, is quantified by the posterior probability
 * that the corresponding toggle is equal to 1.
 */

class DiffSelDoublySparseModel : public ChainComponent {
    // shift/mask counts for an array of sites
    class MaskCounts {
        struct mask_counts_t {
            int nshift{0};
            int nmask{0};
        };

        vector<mask_counts_t> counts;
        mask_counts_t totals;

        void add_shift(int site) {
            assert(nb_active(site) > 1);
            totals.nshift++;
            counts.at(site).nshift++;
        }
        void remove_shift(int site) {
            assert(nb_active(site) > 1);
            totals.nshift--;
            counts.at(site).nshift--;
        }

      public:
        MaskCounts(DiffSelDoublySparseModel &model, int condition) {
            counts.reserve(model.Nsite);  // pre-allocating for model.Nsite sites

            for (int site = 0; site < model.Nsite; site++) {  // for all sites
                counts.emplace_back();
                auto &current_count = counts.back();

                for (int aa = 0; aa < Naa; aa++) {  // for all aas
                    int active = model.sitemaskarray->GetVal(site).at(aa);
                    int convergent = active * model.get_toggle(condition, site, aa);
                    current_count.nmask += active;
                    current_count.nshift += convergent;
                }

                if (current_count.nmask > 1) {
                    totals.nmask += current_count.nmask;
                    totals.nshift += current_count.nshift;
                }
            }
        }

        int nshift() const { return totals.nshift; }
        int nmask() const { return totals.nmask; }
        int nb_active(int site) { return counts.at(site).nmask; }
        void update_toggle(int site, int new_value) {
            if (new_value == 1) {
                add_shift(site);
            } else {
                remove_shift(site);
            }
        }

        bool check(DiffSelDoublySparseModel &model, int condition) const {
            MaskCounts other(model, condition);
            if (other.counts.size() != counts.size()) { return false; }
            if (totals.nshift != other.totals.nshift) { return false; }
            if (totals.nmask != other.totals.nmask) { return false; }
            for (size_t i = 0; i < counts.size(); i++) {
                if (other.counts.at(i).nshift != counts.at(i).nshift) { return false; }
                if (other.counts.at(i).nmask != counts.at(i).nmask) { return false; }
            }
            return true;
        }
    };

    // -----
    // model selectors
    // -----

    const DiffselDoubleSparseConfig config;

    // -----
    // external parameters
    // -----

    unique_ptr<const Tree> tree;
    unique_ptr<FileSequenceAlignment> data;
    unique_ptr<CodonSequenceAlignment> codondata;

    // number of sites
    int Nsite;
    int Ntaxa;
    int Nbranch;

    int Ncond, Nlevel;  // copies of config, to simplify things

    // which branch is under which condition
    unique_ptr<SimpleBranchArray<int>> branchalloc;

    // -----
    //  model structure
    // -----

    per_cond<MaskCounts> mask_counts;

    // Site-wise nodes and toggles
    bool site_wise{false};
    // double sw_toggle_hypermean{0.0};
    // double sw_toggle_hyperinvshape{0.0};
    per_cond<double> sw_toggle_prob;
    per_cond<per_site<indicator_t>> sw_toggles;

    // Branch lengths
    double blhypermean;
    double blhyperinvshape;
    unique_ptr<SimpleBranchArray<double>> blhypermeanarray;
    unique_ptr<GammaWhiteNoise> branchlength;
    unique_ptr<PoissonSuffStatBranchArray> lengthpathsuffstatarray;

    // nucleotide exchange rates and equilibrium frequencies (stationary
    // probabilities) hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    // parameters
    vector<double> nucrelrate;
    vector<double> nucstat;
    unique_ptr<GTRSubMatrix> nucmatrix;

    double fitnessshape;
    vector<double> fitnesscenter;
    unique_ptr<BidimIIDMultiGamma> fitness;

    int maskepsilonmode;
    mask_mode_t mask_mode;
    double maskprob;
    double maskepsilon;
    unique_ptr<IIDProfileMask> sitemaskarray;

    // shiftprob (across conditions):
    // either Beta(shiftprobhypermean,shiftprobhyperinvconc), estimated across
    // genes or mixture 1-pi * 0 + pi * Beta: and this, for each condition
    // separately
    double pihypermean;
    double shiftprobmean;
    double shiftprobinvconc;
    vector<double> pi;
    vector<double> shiftprobhypermean;
    vector<double> shiftprobhyperinvconc;
    vector<double> shiftprob;

    unique_ptr<BidimIIDMultiBernoulli> toggle;

    // fitness profiles (combinations of baseline and delta)
    // across conditions and across sites
    unique_ptr<DiffSelDoublySparseFitnessArray> fitnessprofile;

    // codon substitution matrices
    // across conditions and sites
    unique_ptr<AADiffSelCodonMatrixBidimArray> condsubmatrixarray;

    // branch- and site-substitution matrices (for phyloprocess)
    unique_ptr<SubMatrixSelector> submatrixarray;
    // and for root (condition 0)
    unique_ptr<RootSubMatrixSelector> rootsubmatrixarray;

    // phyloprocess
    unique_ptr<PhyloProcess> phyloprocess;

    // suff stats

    // path suff stats across conditions and sites
    unique_ptr<PathSuffStatBidimArray> suffstatarray;

    MultiGammaSuffStat hyperfitnesssuffstat;

    int gammanullcount;

  public:
    friend ostream &operator<<(ostream &os, DiffSelDoublySparseModel &m);

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
    //! - inepsilon: background fitness for low-fitness amino-acids: if
    //! 0<inepsilon<1, then epsilon is fixed, if epsilon == 1, then this is the
    //! model without masks, if epsilon == -1, then epsilon is estimated from the
    //! data
    //! - inshape: shape parameter of the Gamma distribution of pre-fitness
    //! parameters. If inshape>0, shape parameter is fixed, if inshape == -1,
    //! shape parameter is estimated
    //! - withtoggle: false toggles all fixed to 0, true : random toggles
    DiffSelDoublySparseModel(DiffselDoubleSparseConfig config)
        : config(config), Ncond(config.Ncond), Nlevel(config.Nlevel), hyperfitnesssuffstat(Naa) {
        /* -- */

        if (config.fitness_shape.mode() == param_mode_t::fixed) {
            fitnessshape = config.fitness_shape.value();
        } else {
            fitnessshape = 20.0;
        }

        if (config.epsilon >= 0) {
            maskepsilon = config.epsilon;
            maskepsilonmode = 3;
            mask_mode = gene_spec_mask_fixed_hyper;
        } else {
            maskepsilonmode = 0;
            mask_mode = gene_spec_mask_fixed_hyper;
            maskepsilon = 0.01;
        }

        if (inepsilon == 1) {
            maskepsilon = 1;
            mask_mode = no_mask;
            maskepsilonmode = 3;
        } else {
            maskepsilonmode = 0;
            mask_mode = gene_spec_mask_fixed_hyper;
            maskepsilon = 0.01;
        }

        assert(Nlevel == 1 or Nlevel == 2);
        if (Ncond <= 2 and Nlevel == 2) {
            WARNING("Nlevel set to 2 although there are no more than two conditions");
            Nlevel = 1;
        }

        ReadFiles(config.datafile, config.treefile);
        Allocate();
    }

    DiffSelDoublySparseModel(const DiffSelDoublySparseModel &) = delete;

    ~DiffSelDoublySparseModel() = default;

    //! read files (and read out the distribution of conditions across branches,
    //! based on the tree read from treefile)
    void ReadFiles(string datafile, string treefile) {
        INFO("Parsing nucleotide sequence alignment...");
        data = std::make_unique<FileSequenceAlignment>(datafile);

        INFO("Translating to codons...");
        codondata = std::make_unique<CodonSequenceAlignment>(data.get(), true);
        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();
        INFO("Alignment has {} sites and {} taxa", Nsite, Ntaxa);

        INFO("Parsing tree...");
        std::ifstream file(treefile);
        NHXParser parser{file};
        tree = make_from_parser(parser);
        Nbranch = tree->nb_nodes() - 1;

        INFO("Building branch alloc...");
        per_cond<int> leaf_counts(Ncond, 0);
        auto condition_array = branch_container_from_parser<int>(
            parser, [this, &leaf_counts](int i, const AnnotatedTree &t) -> int {
                auto condition = atoi(t.tag(i, "Condition").c_str());
                if (t.children(i).size() == 0) { leaf_counts.at(condition)++; }
                if (condition >= Ncond) {
                    WARNING(
                        "Found condition {} in tree although model has Ncond={}. Setting condition "
                        "to {}",
                        condition, Ncond, Ncond - 1);
                    return Ncond - 1;
                } else {
                    return condition;
                }
            });
        if (Ncond == 2) {
            INFO("Ncond=2, found {} taxa in condition 0 and {} in condition 1", leaf_counts.at(0),
                leaf_counts.at(1));
        }
        // vector<int> iv(v.size(), 0);
        // for (size_t i = 0; i < v.size(); i++) {
        //     iv[i] = atoi(v[i].c_str());
        //     if (iv[i] >= Ncond) { iv[i] = Ncond - 1; }
        // }
        branchalloc = std::make_unique<SimpleBranchArray<int>>(*tree, condition_array);
    }

    //! allocate the model (data structures)
    void Allocate() {
        // ----------
        // construction of the model
        // ----------
        // allocating data structures and sampling initial configuration

        // Branch lengths
        blhypermean = 0.1;
        blhyperinvshape = 1.0;
        blhypermeanarray = std::make_unique<SimpleBranchArray<double>>(*tree, blhypermean);
        branchlength = std::make_unique<GammaWhiteNoise>(*tree, *blhypermeanarray, blhyperinvshape);
        lengthpathsuffstatarray = std::make_unique<PoissonSuffStatBranchArray>(*tree);

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;

        // nucleotide mutation matrix
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));
        nucmatrix = std::make_unique<GTRSubMatrix>(Nnuc, nucrelrate, nucstat, true);

        // fitness parameters: IID Gamma, across all conditions, sites, and
        // amino-acids those are not the final fitness values (depends on the system
        // of masks and toggles, specified below)
        fitnesscenter.assign(Naa, 1.0 / Naa);
        fitness =
            std::make_unique<BidimIIDMultiGamma>(Ncond, Nsite, Naa, fitnessshape, fitnesscenter);

        // profiles across sites are masked:
        // each site has a 20-dim mask, iid bernoulli of prob maskprob, conditional
        // on at least one entry being 1 if mask[a] == 0 for amino-acid a, then its
        // fitness is equal to maskepsilon, across all conditions, for that site
        maskprob = 0.1;
        sitemaskarray = std::make_unique<IIDProfileMask>(Nsite, Naa, maskprob);

        if (site_wise) {
            set_all_to(sw_toggle_prob, Ncond, 0.2);
            draw_bernoulli_iid(sw_toggles, Ncond, Nsite, sw_toggle_prob);
        } else {
            // hyperparameters for the system of toggles
            // shiftprob is a vector of Ncond-1 probabilities;
            // for each non-baseline condition, k=1..Ncond:
            //     - with probability 1-pi[k-1], shiftprob[k-1] = 0
            //     - with probability pi[k-1]  , shiftprob[k-1] ~
            //     Beta(shiftprobhypermean, shiftprobhyperinvconc)
            pi.assign(Ncond - 1, pihypermean);
            shiftprobhypermean.assign(Ncond - 1, shiftprobmean);
            shiftprobhyperinvconc.assign(Ncond - 1, shiftprobinvconc);
            shiftprob.assign(Ncond - 1, shiftprobmean);

            // toggles specifying the sites and amino-acids displaying fitness
            // modulations across conditions for each k=1..Ncond; all toggles across all
            // sites and amino-acids are iid Bernoulli of parameter shiftprob[k-1]
            toggle = std::make_unique<BidimIIDMultiBernoulli>(Ncond - 1, Nsite, Naa, shiftprob);
            if (not config.withtoggle) { toggle->Reset(); }
        }

        // final amino-acid fitness profiles across sites
        // (*fitnessprofile)(k,i)[a]: fitness of amino-acid a for site i under
        // condition k deterministic functions of the gamma fitness parameters, the
        // system of masks across sites and toggles across sites and conditions for
        // site i and amino-acid a: if sitemaskarray[i][a] == 0 :
        // fitnessprofile(k,i)[a] = maskepsilon across all conditions if
        // sitemaskarray[i][a] == 1 : fitnessprofile(k,i)[a] determined by the
        // system of fitness(0:Ncond,i)[a] and toggle(1:Ncond,i)[a] as in the simple
        // DiffSelSparseModel
        fitnessprofile = std::make_unique<DiffSelDoublySparseFitnessArray>(*fitness, *sitemaskarray,
            [this](int k, int i, int aa) { return get_toggle(k, i, aa); }, Nlevel, maskepsilon);

        // codon matrices
        // per condition and per site
        condsubmatrixarray = std::make_unique<AADiffSelCodonMatrixBidimArray>(
            *fitnessprofile, *GetCodonStateSpace(), *nucmatrix);

        // sub matrices per branch and per site
        submatrixarray = std::make_unique<SubMatrixSelector>(*condsubmatrixarray, *branchalloc);
        // sub matrices for root, across sites
        rootsubmatrixarray = std::make_unique<RootSubMatrixSelector>(*condsubmatrixarray);

        // create phyloprocess
        // TODO: FIX this mysterious nullptr?
        phyloprocess = std::make_unique<PhyloProcess>(tree.get(), codondata.get(),
            branchlength.get(), nullptr, submatrixarray.get(), rootsubmatrixarray.get());
        phyloprocess->Unfold();

        // create suffstat arrays
        suffstatarray = std::make_unique<PathSuffStatBidimArray>(Ncond, Nsite);

        // gathering mask counts
        for (int k = 1; k < Ncond; k++) {
            assert(mask_counts.size() == 0);
            mask_counts.emplace_back(*this, k);
        }
    }

    // --------------------------------------------------------
    // toggle accessors (should be removed ideally)
    // --------------------------------------------------------
    indicator_t &get_toggle(int condition, int site, int aa) {
        assert(condition >= 0 and condition < Ncond);
        assert_warn(condition != 0, "Querying toggle for condition 0");
        return site_wise ? sw_toggles[condition][site] : (*toggle)(condition - 1, site)[aa];
    }

    const indicator_t &get_toggle(int condition, int site, int aa) const {
        assert(condition >= 0 and condition < Ncond);
        assert_warn(condition != 0, "Querying toggle for condition 0");
        return site_wise ? sw_toggles[condition][site] : (*toggle)(condition - 1, site)[aa];
    }

    indicator_t &get_toggle(int condition, int site) {
        assert(condition >= 0 and condition < Ncond);
        assert(site_wise);
        assert_warn(condition != 0, "Querying toggle for condition 0");
        return sw_toggles[condition][site];
    }

    const indicator_t &get_toggle(int condition, int site) const {
        assert(condition >= 0 and condition < Ncond);
        assert(site_wise);
        assert_warn(condition != 0, "Querying toggle for condition 0");
        return sw_toggles[condition][site];
    }

    // --------------------------------------------------------
    // mask count functions
    // --------------------------------------------------------
    MaskCounts &get_mask_counts(int condition) { return mask_counts.at(condition - 1); }

    void update_mask_counts(int condition) {
        get_mask_counts(condition) = MaskCounts(*this, condition);
    }

    bool check_mask_counts(int condition) {
        return get_mask_counts(condition).check(*this, condition);
    }

    // --------------------------------------------------------
    // resample fitness
    // --------------------------------------------------------
    void resample_fitness(int cond, int site, int aa) {
        auto &fitness_ref = (*fitness)(cond, site)[aa];
        fitness_ref = Random::sGamma(fitnessshape * fitnesscenter[aa]);
        if (fitness_ref == 0) {
            gammanullcount++;
            fitness_ref = 1e-8;
        }
    }

    void resample_fitness(int cond, int site) {
        for (int aa = 0; aa < Naa; aa++) {
            if (sitemaskarray->GetVal(site).at(aa)) { resample_fitness(cond, site, aa); }
        }
    }

    // ------------------
    // Update system
    // ------------------
    void Update() {
        UpdateMask();
        fitness->SetShape(fitnessshape);
        UpdateAll();
        ResampleSub(1.0);
    }

    void PostPred(string name) {
        UpdateMask();
        fitness->SetShape(fitnessshape);
        UpdateAll();
        phyloprocess->PostPredSample(name);
    }

    //! update mask array
    void UpdateMask() { sitemaskarray->SetPi(maskprob); }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of Move::Scaling, Move::Sliding and Move::Profile
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

    //! update fitness profiles and matrices across all sites and conditions
    void UpdateAll() {
        fitnessprofile->Update();
        CorruptMatrices();
    }

    //! update fitness profiles and matrices across all conditions for site i
    void UpdateSite(int i) {
        fitnessprofile->UpdateColumn(i);
        condsubmatrixarray->CorruptColumn(i);
    }

    // ---------------
    // log priors
    // ---------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = 0;
        if (config.branch_lengths.mode() != shared) { total += BranchLengthsLogPrior(); }
        if (config.nucmode != shared) { total += NucRatesLogPrior(); }
        if (config.fitness_shape.resampled() or config.fitness_center.resampled()) {
            total += FitnessHyperLogPrior();
        }
        // not updated at all times
        // total += FitnessLogPrior();
        if (gene_specific_mask_mode(config.mask_prob.mode)) {
            total += MaskHyperLogPrior();
            total += MaskLogPrior();
        }
        total += ToggleHyperLogPrior();
        // not updated at all times
        // total += ToggleLogPrior();
        return total;
    }

    //! log prior over branch lengths
    double BranchLengthsLogPrior() const {
        double ret = branchlength->GetLogProb();
        return ret;
    }

    //! log prior over nuc rates rho and pi (uniform)
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(
            nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //! log prior over fitness hyperparameters
    double FitnessHyperLogPrior() const {
        double ret = 0;
        if (config.fitness_shape.resampled()) { ret += FitnessShapeLogPrior(); }
        if (config.fitness_center.resampled()) { ret += FitnessCenterLogPrior(); }
        return ret;
    }

    //! log prior over fitness center hyperparameter
    double FitnessCenterLogPrior() const { return 0; }

    //! log prior over fitness shape hyperparameter
    double FitnessShapeLogPrior() const {
        // exponential on shape
        return -fitnessshape;
    }

    //! log prior over input fitness parameters
    double FitnessLogPrior() const { return fitness->GetLogProb(); }

    //! log prior over mask array hyperparameter (maskprob: uniform between 0 and
    //! 1 -- could be hyperparameterized)
    double MaskHyperLogPrior() const {
        if (maskepsilon > 1) { return Random::INFPROB; }
        return 0;
    }

    //! log prior over mask arrays across sites (iid Bernoulli conditional on at
    //! least one entry being 1)
    double MaskLogPrior() const { return sitemaskarray->GetLogProb(); }

    //! log prior over mask array for site i
    double MaskLogPrior(int i) const { return sitemaskarray->GetLogProb(i); }

    //! log prior over toggle array hyperparameters (shiftprob vector)
    double ToggleHyperLogPrior() const {
        if (not site_wise) {
            double total = 0;
            for (int k = 1; k < Ncond; k++) {
                if (shiftprobhyperinvconc[k - 1]) {
                    double alpha = shiftprobhypermean[k - 1] / shiftprobhyperinvconc[k - 1];
                    double beta = (1 - shiftprobhypermean[k - 1]) / shiftprobhyperinvconc[k - 1];
                    if (shiftprob[k - 1] != 0) {
                        total +=
                            log(pi[k - 1]) + Random::logBetaDensity(shiftprob[k - 1], alpha, beta);
                    } else {
                        assert(pi[k - 1] != 1);
                        total += log(1 - pi[k - 1]);
                    }
                }
            }
            return total;
        } else {
            return 0;
        }
    }

    //! log prior over toggle array (IID bernoulli)
    double ToggleLogPrior() const {
        assert(not site_wise);
        return toggle->GetLogProb();
    }

    //! return log likelihood
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    // ---------------
    // collecting suff stats
    // ---------------

    // //! \brief const access to array of length-pathsuffstats across branches
    // const PoissonSuffStatBranchArray *GetLengthPathSuffStatArray() const {
    //     return lengthpathsuffstatarray;
    // }

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

    //! return log prob of current fitness parameters, conditional on their
    //! hyperparameters
    double FitnessHyperSuffStatLogProb() const {
        return hyperfitnesssuffstat.GetLogProb(fitnessshape, fitnesscenter);
    }

    // ---------------
    // log probs for MH moves
    // ---------------

    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    double NucRatesLogProb() const { return NucRatesLogPrior() + SuffStatLogProb(); }

    //! \brief log prob factor to be recomputed when moving fitness
    //! hyperparameters
    double FitnessHyperLogProb() const {
        return FitnessHyperLogPrior() + FitnessHyperSuffStatLogProb();
    }

    //! \brief log prob factor to be recomputed when moving mask hyperparameter pi
    double MaskLogProb() const { return MaskHyperLogPrior() + MaskLogPrior(); }

    //! \brief log prob factor to be recomputed when moving maskepsilon
    double MaskEpsilonLogProb() const { return MaskHyperLogPrior() + SuffStatLogProb(); }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
    void move(int) override {
        gammanullcount = 0;
        ResampleSub(1.0);
        MoveParameters(3, 20);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep0, int nrep) {
        for (int rep0 = 0; rep0 < nrep0; rep0++) {
            if (config.branch_lengths.mode() != shared) { MoveBranchLengths(); }

            CollectPathSuffStat();
            UpdateAll();

            int weight = 10;

            for (int rep = 0; rep < nrep; rep++) {
                MoveBaselineFitness(weight);
                CompMoveFitness(weight);
                if (maskmode != no_mask) { MoveMasks(weight); }
                if (gene_specific_mask_mode(maskmode)) {
                    MoveMaskHyperParameters(10 * weight);
                }  // FIXME: why move if fixed?
                if (config.withtoggle) {
                    MoveFitnessShifts(weight);
                    if (site_wise) {
                        for (int cond = 1; cond < Ncond; cond++) { move_sw_toggles(cond, 10); }
                    } else {
                        move_shift_toggles(weight);
                    }
                }
                if (config.fitness_shape.resampled() or config.fitness_center.resampled()) {
                    MoveFitnessHyperParameters(10 * weight);
                }
                if (maskepsilonmode < 2) { MoveMaskEpsilon(weight); }
            }

            if (config.nucmode != shared) { MoveNucRates(weight); }
        }

        UpdateAll();
    }

    //! Gibbs resampling of substitution histories conditional on current
    //! parameter configuration
    void ResampleSub(double frac) {
        CorruptMatrices();
        phyloprocess->Move(frac);
    }

    //! Gibbs resampling of branch lengths (based on sufficient statistics)
    void ResampleBranchLengths() {
        CollectLengthSuffStat();
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    //! MCMC move schedule on branch lengths
    void MoveBranchLengths() { ResampleBranchLengths(); }

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! Move::Profile)
    void MoveNucRates(int nrep) {
        CorruptMatrices();

        Move::Profile(nucrelrate, 0.1, 1, nrep, &DiffSelDoublySparseModel::NucRatesLogProb,
            &DiffSelDoublySparseModel::CorruptMatrices, this);
        Move::Profile(nucrelrate, 0.03, 3, nrep, &DiffSelDoublySparseModel::NucRatesLogProb,
            &DiffSelDoublySparseModel::CorruptMatrices, this);
        Move::Profile(nucrelrate, 0.01, 3, nrep, &DiffSelDoublySparseModel::NucRatesLogProb,
            &DiffSelDoublySparseModel::CorruptMatrices, this);

        Move::Profile(nucstat, 0.1, 1, nrep, &DiffSelDoublySparseModel::NucRatesLogProb,
            &DiffSelDoublySparseModel::CorruptMatrices, this);
        Move::Profile(nucstat, 0.01, 1, nrep, &DiffSelDoublySparseModel::NucRatesLogProb,
            &DiffSelDoublySparseModel::CorruptMatrices, this);

        CorruptMatrices();
    }

    //! MH compensatory move schedule on fitness parameters and hyper-parameters
    void CompMoveFitness(int nrep) { CompMoveFitness(1.0, nrep); }

    //! \brief MH compensatory move on fitness parameters and hyper-parameters
    //!
    //! for a given amino-acid.
    //! shift by a constant factor e all *active* gamma fitness parameters across
    //! all amino-acids and conditions
    double CompMoveFitness(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                const vector<int> &mask = (*sitemaskarray)[i];

                double deltalogprob = 0;

                // calculate log prior for active fitness parameters across amino-acids
                // and conditions before the move
                for (int k = 0; k < Ncond; k++) {
                    for (int a = 0; a < Naa; a++) {
                        if ((mask[k]) && ((!k) || (get_toggle(k, i, a)))) {
                            double alpha = fitnessshape * fitnesscenter[a];
                            deltalogprob -= -Random::logGamma(alpha) +
                                            (alpha - 1) * log((*fitness)(k, i)[a]) -
                                            (*fitness)(k, i)[a];
                        }
                    }
                }

                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);

                // multiply all active fitness parameters across amino-acids and
                // conditions by e
                int n = 0;
                for (int k = 0; k < Ncond; k++) {
                    for (int a = 0; a < Naa; a++) {
                        if ((mask[k]) && ((!k) || (get_toggle(k, i, a)))) {
                            (*fitness)(k, i)[a] *= e;
                            n++;
                        }
                    }
                }

                double loghastings = n * m;

                // calculate log prior for active fitness parameters across amino-acids
                // and conditions after the move
                for (int k = 0; k < Ncond; k++) {
                    for (int a = 0; a < Naa; a++) {
                        if ((mask[k]) && ((!k) || (get_toggle(k, i, a)))) {
                            double alpha = fitnessshape * fitnesscenter[a];
                            deltalogprob += -Random::logGamma(alpha) +
                                            (alpha - 1) * log((*fitness)(k, i)[a]) -
                                            (*fitness)(k, i)[a];
                        }
                    }
                }

                deltalogprob += loghastings;

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    // restore previous value for active fitness parameters
                    for (int k = 0; k < Ncond; k++) {
                        for (int a = 0; a < Naa; a++) {
                            if ((mask[k]) && ((!k) || (get_toggle(k, i, a)))) {
                                (*fitness)(k, i)[a] /= e;
                            }
                        }
                    }
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on baseline gamma fitness parameters (for condition k=0)
    void MoveBaselineFitness(int nrep) {
        // if masks are not activated (all entries equal to 1), move a random subset
        // of entries over the 20 amino-acids (2d parameter of call)
        if (maskmode == no_mask) {
            MoveAllBaselineFitness(1.0, 3, nrep);
            MoveAllBaselineFitness(1.0, 10, nrep);
            MoveAllBaselineFitness(1.0, 20, nrep);
            MoveAllBaselineFitness(0.3, 20, nrep);
        }
        // if masks are activated, move all active entries
        else {
            MoveBaselineFitness(1.0, nrep);
            MoveBaselineFitness(0.3, nrep);
        }
    }

    //! elementary MH move on baseline gamma fitness parameters (for condition
    //! k=0)
    double MoveAllBaselineFitness(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                vector<double> &x = (*fitness)(0, i);

                bk = x;

                double deltalogprob = -fitness->GetLogProb(0, i) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, n);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(0, i) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! elementary MH move on baseline fitness parameters (for condition k=0):
    //! version used when masks are activated
    double MoveBaselineFitness(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                vector<double> &fit = (*fitness)(0, i);
                const vector<int> &mask = (*sitemaskarray)[i];

                bk = fit;

                // calculate log prob before the move (for the prior factor, only those
                // entries that are not masked are counted)
                double deltalogprob = -fitness->GetLogProb(0, i, mask) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(fit, Naa, tuning, mask);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(0, i, mask) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    fit = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on gamma fitness parameters (fitness shifts) for
    //! non-baseline conditions
    void MoveFitnessShifts(int nrep) {
        for (int k = 1; k < Ncond; k++) {
            MoveFitnessShifts(k, 1, nrep);
            MoveFitnessShifts(k, 0.3, nrep);
        }
    }

    //! elementary MH move on fitness shifts for non-baseline conditions
    double MoveFitnessShifts(int k, double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                vector<double> &x = (*fitness)(k, i);
                const vector<int> &m = sitemaskarray->GetVal(i);

                // compute condition-specific mask, which is the conjunction of baseline
                // mask and condition-specific vector of toggles: s = m*t this mask
                // specifies which amino-acids are both active (across the tree) and
                // undergoing a fitness shift in current condition
                vector<int> s(Naa, 0);

                // nshift: number of amino-acids that are active in baseline and
                // undergoing a fitness shift in current condition nmask : number of
                // amino-acids that are active in baseline
                int nshift = 0;
                int nmask = 0;
                for (int a = 0; a < Naa; a++) {
                    s[a] = get_toggle(k, i, a) * m[a];
                    nmask += m[a];
                    nshift += s[a];
                }

                // fitness shifts are impacting the final fitness only if
                // (1) there are at least 2 active amino-acids
                // (2) at least one amino-acid is undergoing a fitness shift
                if ((nmask > 1) && nshift) {
                    bk = x;

                    // log prob before the move (for the prior factor, only those
                    // amino-acids that are concerned, such as specified by s, are
                    // counted)
                    double deltalogprob = -fitness->GetLogProb(k, i, s) - SiteSuffStatLogProb(i);

                    // propose move (only for the relevant amino-acids)
                    double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, s);

                    deltalogprob += loghastings;

                    UpdateSite(i);

                    // log prob after the move
                    deltalogprob += fitness->GetLogProb(k, i, s) + SiteSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        x = bk;
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! MH moves on hyperparameters of distribution of fitness factors
    void MoveFitnessHyperParameters(int nrep) {
        // collect suff stats across all active fitness parameters
        hyperfitnesssuffstat.Clear();
        hyperfitnesssuffstat.AddSuffStat(*fitness, *sitemaskarray, *toggle);

        if (config.fitness_shape.resampled()) {
            Move::Scaling(fitnessshape, 1.0, nrep, &DiffSelDoublySparseModel::FitnessHyperLogProb,
                &DiffSelDoublySparseModel::NoUpdate, this);
            Move::Scaling(fitnessshape, 0.3, nrep, &DiffSelDoublySparseModel::FitnessHyperLogProb,
                &DiffSelDoublySparseModel::NoUpdate, this);
            Move::Scaling(fitnessshape, 0.1, nrep, &DiffSelDoublySparseModel::FitnessHyperLogProb,
                &DiffSelDoublySparseModel::NoUpdate, this);
        }
        fitness->SetShape(fitnessshape);

        if (config.fitness_center.resampled()) {
            Move::Profile(fitnesscenter, 0.3, 1, nrep,
                &DiffSelDoublySparseModel::FitnessHyperLogProb, &DiffSelDoublySparseModel::NoUpdate,
                this);
            Move::Profile(fitnesscenter, 0.1, 1, nrep,
                &DiffSelDoublySparseModel::FitnessHyperLogProb, &DiffSelDoublySparseModel::NoUpdate,
                this);
            Move::Profile(fitnesscenter, 0.1, 3, nrep,
                &DiffSelDoublySparseModel::FitnessHyperLogProb, &DiffSelDoublySparseModel::NoUpdate,
                this);
        }
    }

    //! Move schedule for Gibbs resampling of shifting probabilities
    void ResampleShiftProb() {
        if (!shiftprobinvconc) {
            std::cerr << "error: in resample shift prob\n";
            exit(1);
        }
        for (int k = 1; k < Ncond; k++) { ResampleShiftProb(k); }
    }

    //! Gibbs resampling of shifting probability under condition k
    void ResampleShiftProb(int k) {
        assert(not site_wise);
        // pre-calculate parameters of the Beta distribution for non-zero case
        double alpha = shiftprobhypermean[k - 1] / shiftprobhyperinvconc[k - 1];
        double beta = (1 - shiftprobhypermean[k - 1]) / shiftprobhyperinvconc[k - 1];

        // nshift: number of amino-acids that are active in baseline and undergoing
        // a fitness shift in current condition nmask : number of amino-acids that
        // are active in baseline both are summed across all sites: sufficient
        // statistics for shiftprob
        int nshift = 0;
        int nmask = 0;
        for (int i = 0; i < Nsite; i++) {
            // const vector<int> &t = (*toggle)(k - 1, i);
            const vector<int> &m = sitemaskarray->GetVal(i);
            int ns = 0;
            int nm = 0;
            for (int a = 0; a < Naa; a++) {
                ns += m[a] * get_toggle(k, i, a);
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least
            // 2 active amino-acids
            if (nm > 1) {
                nshift += ns;
                nmask += nm;
            }
        }

        if (nshift || (pi[k - 1] == 1.0)) {
            shiftprob[k - 1] = Random::BetaSample(alpha + nshift, beta + nmask - nshift);
        } else {
            double logp0 = log(1 - pi[k - 1]);

            double logp1 = log(pi[k - 1]);
            logp1 -=
                Random::logGamma(alpha) + Random::logGamma(beta) - Random::logGamma(alpha + beta);
            logp1 += Random::logGamma(alpha + nshift) + Random::logGamma(beta + nmask - nshift) -
                     Random::logGamma(alpha + beta + nmask);

            double max = (logp0 > logp1) ? logp0 : logp1;
            double p0 = exp(logp0 - max);
            double p1 = exp(logp1 - max);
            double tot = p0 + p1;
            p0 /= tot;
            p1 /= tot;

            if (Random::Uniform() < p0) {
                shiftprob[k - 1] = 0;
            } else {
                shiftprob[k - 1] = Random::BetaSample(alpha + nshift, beta + nmask - nshift);
            }
        }
    }

    //! empirical fraction of allowed positions that undergo a shift
    double GetPropShift(int k) const {
        // nshift: number of amino-acids that are active in baseline and undergoing
        // a fitness shift in current condition nmask : number of amino-acids that
        // are active in baseline both are summed across all sites: sufficient
        // statistics for shiftprob
        int nshift = 0;
        int nmask = 0;
        for (int i = 0; i < Nsite; i++) {
            const vector<int> &m = sitemaskarray->GetVal(i);
            int ns = 0;
            int nm = 0;
            for (int a = 0; a < Naa; a++) {
                ns += m[a] * get_toggle(k, i, a);
                nm += m[a];
            }

            // fitness shifts are counted (have an effect) only if there are at least
            // 2 active amino-acids
            if (nm > 1) {
                nshift += ns;
                nmask += nm;
            }
        }

        return ((double)nshift) / nmask;
    }

    //! MH move schedule on mask hyperparameter (maskprob)
    void MoveMaskHyperParameters(int nrep) {
        Move::Sliding(maskprob, 1.0, nrep, 0.05, 0.975, &DiffSelDoublySparseModel::MaskLogProb,
            &DiffSelDoublySparseModel::UpdateMask, this);
        Move::Sliding(maskprob, 0.1, nrep, 0.05, 0.975, &DiffSelDoublySparseModel::MaskLogProb,
            &DiffSelDoublySparseModel::UpdateMask, this);
    }

    //! MH move schedule on background fitness (maskepsilon)
    void MoveMaskEpsilon(int nrep) {
        Move::Sliding(maskepsilon, 1.0, nrep, 0, 1.0, &DiffSelDoublySparseModel::MaskEpsilonLogProb,
            &DiffSelDoublySparseModel::UpdateAll, this);
        Move::Sliding(maskepsilon, 0.1, nrep, 0, 1.0, &DiffSelDoublySparseModel::MaskEpsilonLogProb,
            &DiffSelDoublySparseModel::UpdateAll, this);
        Move::Scaling(maskepsilon, 1.0, nrep, &DiffSelDoublySparseModel::MaskEpsilonLogProb,
            &DiffSelDoublySparseModel::UpdateAll, this);
        Move::Scaling(maskepsilon, 0.1, nrep, &DiffSelDoublySparseModel::MaskEpsilonLogProb,
            &DiffSelDoublySparseModel::UpdateAll, this);
    }

    //! MH move on fitness masks across sites
    double MoveMasks(int nrep) {
        double nacc = 0;
        double ntot = 0;

        for (int i = 0; i < Nsite; i++) {
            vector<int> &mask = (*sitemaskarray)[i];

            // compute number of active entries
            int naa = 0;
            for (int a = 0; a < Naa; a++) { naa += mask[a]; }

            for (int rep = 0; rep < nrep; rep++) {
                // randomly choose amino-acid
                int a = (int)(Naa * Random::Uniform());

                // don't propose move if this leads to 0 active entry in the end
                if ((!mask[a]) || (naa > 1)) {
                    // logprob before the move
                    double deltalogprob = -MaskLogPrior(i) - SiteSuffStatLogProb(i);

                    // do the move on the entry mask[k]
                    int oldnaa = naa;
                    naa -= mask[a];
                    mask[a] = 1 - mask[a];
                    naa += mask[a];

                    // if move is from inactive to active
                    if (mask[a]) {
                        // resample baseline fitness
                        resample_fitness(0, i, a);

                        // resample toggles and fitness shifts across all non-baseline
                        // conditions
                        for (int k = 1; k < Ncond; k++) {
                            if (site_wise) {
                                get_toggle(k, i) = (Random::Uniform() < sw_toggle_prob[k - 1]);
                                resample_fitness(k, i);
                            } else {
                                get_toggle(k, i, a) = (Random::Uniform() < shiftprob[k - 1]);
                                if (get_toggle(k, i, a)) { resample_fitness(k, i, a); }
                            }
                        }
                    }

                    // if move is from 1 to 2 active amino-acids,
                    // the other unmasked amino-acid should also redraw its toggles and
                    // fitness shifts (but not the baseline fitness)
                    if ((oldnaa == 1) && (naa == 2)) {
                        if (site_wise) {
                            for (int k = 1; k < Ncond; k++) {
                                get_toggle(k, i) =
                                    indicator_t(Random::Uniform() < sw_toggle_prob.at(k));
                                if (get_toggle(k, i)) { resample_fitness(k, i); }
                            }
                        } else {
                            // grep the index of the other active amino-acid
                            int b = 0;
                            while ((b < Naa) && ((!mask[b]) || (b == a))) { b++; }
                            assert(b != Naa);

                            // resample toggles and fitness shifts across all non-baseline
                            // conditions
                            for (int k = 1; k < Ncond; k++) {
                                get_toggle(k, i, b) =
                                    indicator_t(Random::Uniform() < shiftprob.at(k - 1));
                                if (get_toggle(k, i, b)) { resample_fitness(k, i, b); }
                            }
                        }
                    }

                    // logprob after the move
                    UpdateSite(i);
                    deltalogprob += MaskLogPrior(i) + SiteSuffStatLogProb(i);

                    // decision
                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        // restore state (not necessary to restore all fitness parameters
                        // and toggles that have been modified during the move, since they
                        // all go back to inactive state anyway and will be resampled upon
                        // next activation)
                        naa -= mask[a];
                        mask[a] = 1 - mask[a];
                        naa += mask[a];
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //! MH move schedule on toggles
    void move_shift_toggles(int nrep) {
        for (int k = 1; k < Ncond; k++) { move_shift_toggles(k, nrep); }
    }

    //! helper function: returns the marginal log prob of distribution of toggles
    //! for a condition, given number of toggles in active state and given
    //! hyperparameters
    double ToggleMarginalLogPrior(int nmask, int nshift, int k) const {
        assert(not site_wise);
        double ret = 0;
        if (shiftprobhyperinvconc[k - 1]) {
            // pre-calculate parameters of the Beta distribution for non-zero case
            double alpha = shiftprobhypermean[k - 1] / shiftprobhyperinvconc[k - 1];
            double beta = (1 - shiftprobhypermean[k - 1]) / shiftprobhyperinvconc[k - 1];
            double pp = pi[k - 1];

            double logp1 = log(pp) + Random::logGamma(alpha + beta) - Random::logGamma(alpha) -
                           Random::logGamma(beta) + Random::logGamma(alpha + nshift) +
                           Random::logGamma(beta + nmask - nshift) -
                           Random::logGamma(alpha + beta + nmask);
            if (nshift || (pp == 1.0)) { return logp1; }
            double logp0 = log(1 - pp);
            double max = (logp0 > logp1) ? logp0 : logp1;
            double tot = exp(logp0 - max) + exp(logp1 - max);
            ret = log(tot) + max;
        } else {
            ret = nshift * log(shiftprob[k - 1]) + (nmask - nshift) * (1 - shiftprob[k - 1]);
        }
        return ret;
    }

    //-----------------------------------------
    // Toggle-related move and utilities
    // ----------------------------------------

    // toggle move for the site-wise case
    double move_sw_toggles(int cond, int nrep) {
        assert(site_wise);
        int sw_nb_on = count_indicators(sw_toggles.at(cond));
        update_mask_counts(cond);
        AcceptanceStats acceptance_stats;  // used for nmask only (and thus not updated here)

        for (int rep = 0; rep < nrep; rep++) {
            for (int site = 0; site < Nsite; site++) {
                if (get_mask_counts(cond).nb_active(site) < 2) {
                    double log_prob_before = SiteSuffStatLogProb(site) +
                                             sw_nb_on * log(sw_toggle_prob.at(cond)) +
                                             (Nsite - sw_nb_on) * log(1 - sw_toggle_prob.at(cond));

                    auto &togref = sw_toggles.at(cond).at(site);
                    assert(togref == 1 or togref == 0);
                    togref = 1 - togref;          // change value
                    sw_nb_on += togref ? 1 : -1;  // update toggle count
                    UpdateSite(site);             // update site logprobs

                    if (togref == 1) {  // if toggle turned on then redraw fitness for all aas
                        resample_fitness(cond, site);
                    }

                    double log_prob_after = SiteSuffStatLogProb(site) +
                                            sw_nb_on * log(sw_toggle_prob.at(cond)) +
                                            (Nsite - sw_nb_on) * log(1 - sw_toggle_prob.at(cond));

                    double acceptance_prob = log_prob_after - log_prob_before;
                    if (decide(acceptance_prob)) {
                        acceptance_stats.accept();
                    } else {
                        acceptance_stats.reject();
                        togref = 1 - togref;          // change value
                        sw_nb_on += togref ? 1 : -1;  // update toggle count
                        UpdateSite(site);             // update site logprobs
                    }
                }
            }
        }
        return acceptance_stats.ratio();
    }

    //! elementary MH move on toggles
    double move_shift_toggles(int k, int nrep) {
        assert(not site_wise);
        // to achieve better MCMC mixing, shiftprob[k-1] is integrated out during
        // this MH move on toggles (and Gibbs-resampled upon leaving this MH update)
        // nshift: number of amino-acids that are active in baseline and undergoing
        // a fitness shift in current condition nmask : number of amino-acids that
        // are active in baseline both are summed across all sites: sufficient
        // statistics for shiftprob
        update_mask_counts(k);
        DEBUG("move_shift_toggles k={}; nmask={}; nshift={}", k, get_mask_counts(k).nmask(),
            get_mask_counts(k).nshift());

        AcceptanceStats acceptance_stats;
        for (int rep = 0; rep < nrep; rep++) {  // repeating move nrep times
            for (int i = 0; i < Nsite; i++) {   // for every site...
                assert(check_mask_counts(k));   // checking mask_count consistency (costly)

                const vector<int> &site_mask = sitemaskarray->GetVal(i);

                int nb_active = get_mask_counts(k).nb_active(i);
                // do move only if there are at least 2 active amino-acids
                if (nb_active > 1) {
                    // randomly choose one active amino-acid (for which mask[a] == 1)
                    int nth_active = (int)(nb_active * Random::Uniform()) + 1;
                    int index = 0;
                    while (nth_active && (index < Naa)) {
                        if (site_mask[index] == 1) { nth_active--; }
                        if (nth_active > 0) { index++; }
                    }
                    assert(index != Naa);
                    assert(site_mask[index] > 0);

                    int chosen_aa = index;
                    int &chosen_toggle_ref = get_toggle(k, i, chosen_aa);

                    // core of the MH move (logprobs and value change)
                    double logprob_before = ToggleMarginalLogPrior(get_mask_counts(k).nmask(),
                                                get_mask_counts(k).nshift(), k) +
                                            SiteSuffStatLogProb(i);
                    assert(chosen_toggle_ref == 0 or chosen_toggle_ref == 1);
                    chosen_toggle_ref = 1 - chosen_toggle_ref;  // 1->0 or 0->1

                    if (chosen_toggle_ref == 1) {  // if toggle turned on then redraw fitness
                        resample_fitness(k, i, chosen_aa);
                    }

                    UpdateSite(i);
                    get_mask_counts(k).update_toggle(i, chosen_toggle_ref);

                    double logprob_after = ToggleMarginalLogPrior(get_mask_counts(k).nmask(),
                                               get_mask_counts(k).nshift(), k) +
                                           SiteSuffStatLogProb(i);

                    // MH move accept/reject
                    double acceptance_prob = logprob_after - logprob_before;
                    int accepted = (log(Random::Uniform()) < acceptance_prob);
                    if (accepted) {
                        acceptance_stats.accept();
                    } else {
                        acceptance_stats.reject();
                        chosen_toggle_ref = 1 - chosen_toggle_ref;  // toggle back
                        get_mask_counts(k).update_toggle(i, chosen_toggle_ref);
                        UpdateSite(i);
                    }
                }
            }
        }
        if (shiftprobinvconc) { ResampleShiftProb(k); }
        return acceptance_stats.ratio();
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    const CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //-------------------
    // Traces and monitors
    // ------------------

    //! return mean width of masks across sites
    double GetMeanWidth() const { return sitemaskarray->GetMeanWidth(); }

    //! write complete current parameter configuration to stream
    void ToStream(ostream &os) { os << *this; }

    template <class Info>
    void declare_interface(Info info) {
        if (config.branch_lengths.mode() != shared) {
            model_node(info, "branchlength", *branchlength);
        }
        if (config.nucmode != shared) {
            model_node(info, "nucrelrate", nucrelrate);
            model_node(info, "nucstat", nucstat);
        }
        if (config.fitness_shape.resampled()) { model_node(info, "fitnessshape", fitnessshape); }
        if (config.fitness_center.resampled()) { model_node(info, "fitnesscenter", fitnesscenter); }
        model_node(info, "fitness", *fitness);
        if (gene_specific_mask_mode(maskmode)) { model_node(info, "maskprob", maskprob); }
        if (maskmode != no_mask) { model_node(info, "sitemaskarray", *sitemaskarray); }
        if (maskepsilonmode < 2) { model_node(info, "maskepsilon", maskepsilon); }
        if (Ncond > 1) {
            if (site_wise) {
                model_node(info, "sw_toggles", sw_toggles);
                model_node(info, "sw_toggle_prob", sw_toggle_prob);
            } else {
                model_node(info, "shiftprob", shiftprob);
                model_node(info, "toggle", *toggle);
            }
        }

        model_stat(info, "logprior", [this]() { return GetLogPrior(); });
        model_stat(info, "lnL", [this]() { return GetLogLikelihood(); });
        model_stat(info, "length", [this]() { return 3 * branchlength->GetTotalLength(); });
        model_stat(info, "maskprob", maskprob);
        model_stat(info, "meanwidth", [this]() { return GetMeanWidth(); });
        model_stat(info, "maskepsilon", maskepsilon);
        model_stat(info, "fitnessshape", fitnessshape);
        model_stat(
            info, "fitnesscenter_entropy", [&]() { return Random::GetEntropy(fitnesscenter); });
        for (int k = 1; k < Ncond; k++) {
            if (not site_wise) {
                model_stat(info, "shiftprob_" + std::to_string(k), shiftprob.at(k - 1));
            }
            model_stat(
                info, "propshift_" + std::to_string(k), [k, this]() { return GetPropShift(k); });
        }
        model_stat(info, "nucstat_entropy", [&]() { return Random::GetEntropy(nucstat); });
        model_stat(info, "nucrelrate_entropy", [&]() { return Random::GetEntropy(nucrelrate); });
        model_stat(info, "gammanullcount", gammanullcount);
    }
};

istream &operator>>(istream &is, unique_ptr<DiffSelDoublySparseModel> &m) {
    string model_name;
    string datafile;
    string treefile;
    int Ncond, Nlevel;
    int fitnessshapemode, fitnesscentermode;  // implicit cast of enum into int
    double maskepsilonmode, pihypermean, shiftprobmean, shiftprobinvconc;

    is >> model_name;
    if (model_name != "DiffselDoublySparse") {
        std::cerr << "Expected DiffselDoublySparse for model name, got " << model_name << "\n";
        exit(1);
    }
    is >> datafile;
    is >> treefile;
    is >> Ncond >> Nlevel >> maskepsilonmode >> fitnessshapemode >> pihypermean >> shiftprobmean >>
        shiftprobinvconc >> fitnesscentermode;
    m = std::make_unique<DiffSelDoublySparseModel>(datafile, treefile, Ncond, Nlevel,
        maskepsilonmode, param_mode_t(fitnessshapemode), pihypermean, shiftprobmean,
        shiftprobinvconc, param_mode_t(fitnesscentermode), true, false);  // FIXME site_wise bool!
    assert("implemented" == "false");
    Tracer tracer{*m};
    tracer.read_line(is);
    return is;
}

ostream &operator<<(ostream &os, DiffSelDoublySparseModel &m) {
    Tracer tracer{m};
    os << "DiffselDoublySparse" << '\t';
    os << m.datafile << '\t';
    os << m.treefile << '\t';
    os << m.Ncond << '\t';
    os << m.Nlevel << '\t';
    os << m.maskepsilonmode << '\t';
    os << m.fitnessshapemode << '\t';
    os << m.pihypermean << '\t';
    os << m.shiftprobmean << '\t';
    os << m.shiftprobinvconc << '\t';
    os << m.fitnesscentermode << '\t';
    tracer.write_line(os);
    return os;
}
