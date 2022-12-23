
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "Product.hpp"
#include "ConditionSpecificMeanGammaArray.hpp"
#include "dSOmegaPathSuffStatGeneBranchArray.hpp"
#include "RecursiveNewick.hpp"

class FastGeneBranchOmegaModel : public ProbModel {

  private:

    int syn_devmode;
    int om_devmode;

    Tree *tree;
    // const TaxonSet *taxonset;
    int Ntaxa;
    int Nbranch;
    int Ngene;

    double branchsyn_hypermean;
    double branchsyn_hyperinvshape;
    IIDGamma *branchsyn_array;
    GammaSuffStat branchsyn_hypersuffstat;

    double genesyn_hypermean; // constrained = 1.0
    double genesyn_hyperinvshape;
    IIDGamma *genesyn_array;
    GammaSuffStat genesyn_hypersuffstat;

    double syn_invshape;
    ProductArray *meansyn_bidimarray;
    ConditionSpecificMeanGammaBidimArray *syn_bidimarray;

    double branchom_hypermean;
    double branchom_hyperinvshape;
    IIDGamma *branchom_array;
    GammaSuffStat branchom_hypersuffstat;

    double geneom_hypermean; // constrained = 1.0
    double geneom_hyperinvshape;
    IIDGamma *geneom_array;
    GammaSuffStat geneom_hypersuffstat;

    double om_invshape;
    ProductArray *meanom_bidimarray;
    ConditionSpecificMeanGammaBidimArray *om_bidimarray;

    dSOmegaPathSuffStatGeneBranchArray *dsomss;
    vector<string> gene_names;

  public:

    FastGeneBranchOmegaModel(string datafile, string treefile, int insyn_devmode, int inom_devmode)  {

        syn_devmode = insyn_devmode;
        om_devmode = inom_devmode;

        tree = new Tree(treefile);
        tree->SetIndices();
        Ntaxa = tree->GetSize();
        Nbranch = tree->GetNbranch();

        ifstream is(datafile.c_str());
        is >> Ngene;
        gene_names.assign(Ngene,"");
        Allocate();

        // read suff stats
        for (int gene=0; gene<Ngene; gene++)    {
            is >> gene_names[gene];

            string tmp;
            is >> tmp;
            if (tmp != "counts_dS") {
                cerr << "error when reading suffstat file\n";
                exit(1);
            }
            Tree treedscount(is);
            treedscount.SetIndices();

            is >> tmp;
            if (tmp != "counts_dS_norm") {
                cerr << "error when reading suffstat file\n";
                exit(1);
            }
            Tree treedsbeta(is);
            treedsbeta.SetIndices();

            is >> tmp;
            if (tmp != "counts_dN") {
                cerr << "error when reading suffstat file\n";
                exit(1);
            }
            Tree treedncount(is);
            treedncount.SetIndices();

            is >> tmp;
            if (tmp != "counts_dN_norm") {
                cerr << "error when reading suffstat file\n";
                exit(1);
            }
            Tree treednbeta(is);
            treednbeta.SetIndices();

            (*dsomss)[gene].Add(treedscount, treedsbeta, treedncount, treednbeta);
        }
    }

    int GetNtaxa() const    {
        return Ntaxa;
    }

    int GetNbranch() const  {
        return Nbranch;
    }

    int GetNgene() const    {
        return Ngene;
    }

    string GetGeneName(int gene) const  {
        return gene_names[gene];
    }

    const Tree* GetTree() const {
        return tree;
    }

    void Allocate() {

        branchsyn_hypermean = 0.1;
        branchsyn_hyperinvshape = 1.0;
        double branchsyn_alpha = 1.0 / branchsyn_hyperinvshape;
        double branchsyn_beta = branchsyn_alpha / branchsyn_hypermean;
        branchsyn_array = new IIDGamma(Nbranch, branchsyn_alpha, branchsyn_beta);
        for (int j=0; j<Nbranch; j++) {
            (*branchsyn_array)[j] = 1.0;
        }

        genesyn_hypermean = 1.0;
        genesyn_hyperinvshape = 1.0;
        double genesyn_alpha = 1.0 / genesyn_hyperinvshape;
        double genesyn_beta = genesyn_alpha / genesyn_hypermean;
        genesyn_array = new IIDGamma(Ngene, genesyn_alpha, genesyn_beta);

        meansyn_bidimarray = new ProductArray(*branchsyn_array, *genesyn_array);
        syn_invshape = 1.0;

        syn_bidimarray = 0;
        if (syn_devmode)    {
            syn_bidimarray =
                new ConditionSpecificMeanGammaBidimArray(*meansyn_bidimarray, syn_invshape);
        }

        branchom_hypermean = 1.0;
        branchom_hyperinvshape = 1.0;
        double branchom_alpha = 1.0 / branchom_hyperinvshape;
        double branchom_beta = branchom_alpha / branchom_hypermean;
        branchom_array = new IIDGamma(Nbranch, branchom_alpha, branchom_beta);
        for (int j=0; j<Nbranch; j++) {
            (*branchom_array)[j] = 1.0;
        }

        geneom_hypermean = 1.0;
        geneom_hyperinvshape = 1.0;
        double geneom_alpha = 1.0 / geneom_hyperinvshape;
        double geneom_beta = geneom_alpha / geneom_hypermean;
        geneom_array = new IIDGamma(Ngene, geneom_alpha, geneom_beta);

        meanom_bidimarray = new ProductArray(*branchom_array, *geneom_array);
        om_invshape = 1.0;

        om_bidimarray = 0;
        if (om_devmode) {
            om_bidimarray =
                new ConditionSpecificMeanGammaBidimArray(*meanom_bidimarray, om_invshape);
        }

        dsomss = new dSOmegaPathSuffStatGeneBranchArray(*tree, Ngene);
    }

    /*
    void PrintdNdS(ostream& os) const  {
        dSOmegaPathSuffStatBranchArray globdsomss(*tree);
        for (int i=0; i<Ngene; i++) {
            globdsomss.Add((*dsomss)[i]);
        }
        SimpleBranchArray<double> dnds(*tree);
        globdsomss.GetdNdS(dnds);
        Tabulate(os, dnds, true);
        os << '\n';
        dsomss->ToStream(os);
    }
    */

    void TraceHeader(ostream &os) const override {
        os << "logprior\tlogl";

        os << "\tlength\tinvshape";
        os << "\tgeneds_mean";
        os << "\tgeneds_invshape";

        os << "\tgeneom\tgeneom_invshape";
        os << "\tbranchom_mean";
        os << "\tbranchom_invshape";

        if (syn_devmode)    {
            os << "\tomdev";
        }
        if (om_devmode) {
            os << "\tsyndev";
        }
        os << '\n';
    }


    double GetLength() const    {
        double tot = 0;
        for (int j=0; j<Nbranch; j++)   {
            tot += branchsyn_hypermean;
        }
        return tot;
    }

    double GetGeneMeanSyn() const   {
        double tot = 0;
        for (int i=0; i<Ngene; i++) {
            tot += genesyn_array->GetVal(i);
        }
        tot /= Ngene;
        return tot;
    }

    double GetBranchMeanOmega() const   {
        double tot = 0;
        for (int j=0; j<Nbranch; j++)   {
            tot += branchom_array->GetVal(j);
        }
        tot /= Nbranch;
        return tot;
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t' << GetLogLikelihood();

        os << '\t' << GetLength() << '\t' << branchsyn_hyperinvshape;
        os << '\t' << GetGeneMeanSyn();
        os << '\t' << genesyn_hyperinvshape;

        os << '\t' << geneom_hypermean << '\t' << geneom_hyperinvshape;
        os << '\t' << GetBranchMeanOmega();
        os << '\t' << branchom_hyperinvshape;

        if (syn_devmode)    {
            os << '\t' << om_invshape;
        }

        if (om_devmode) {
            os << '\t' << syn_invshape;
        }

        os << '\n';
    }

    void Monitor(ostream &os) const override {
    }

    void ToStream(ostream &os) const override {

        os << branchsyn_hypermean << '\t' << branchsyn_hyperinvshape;
        os << '\t' << *branchsyn_array;
        os << '\t' << genesyn_hyperinvshape;
        os << '\t' << *genesyn_array;


        os << '\t' << branchom_hyperinvshape;
        os << '\t' << *branchom_array;
        os << '\t' << geneom_hypermean << '\t' << geneom_hyperinvshape;
        os << '\t' << *geneom_array;

        if (syn_devmode)    {
            os << '\t' << syn_invshape;
        }
        if (om_devmode) {
            os << '\t' << om_invshape;
        }

        /*
        if (syn_devmode)    {
            os << '\t' << *syn_bidimarray;
        }
        if (om_devmode) {
            os << '\t' << *om_bidimarray;
        }
        */
    }

    void FromStream(istream &is) override {

        is >> branchsyn_hypermean >> branchsyn_hyperinvshape;
        is >> *branchsyn_array;
        is >> genesyn_hyperinvshape;
        is >> *genesyn_array;

        is >> branchom_hyperinvshape;
        is >> *branchom_array;
        is >> geneom_hypermean >> geneom_hyperinvshape;
        is >> *geneom_array;

        if (syn_devmode)    {
            is >> syn_invshape;
        }
        if (om_devmode) {
            is >> om_invshape;
        }

        /*
        if (syn_devmode)    {
            is >> *syn_bidimarray;
        }
        if (om_devmode) {
            is >> *om_bidimarray;
        }
        */
    }

    void AddGeneOmegaArrayTo(vector<double>& array) const {
        for (int i=0; i<Ngene; i++)   {
            array[i] += geneom_array->GetVal(i);
        }
    }

    void AddBranchOmegaArrayTo(vector<double>& array) const {
        for (int j=0; j<Nbranch; j++)   {
            array[j] += branchom_array->GetVal(j);
        }
    }

    void Update(IIDGamma& array, double hypermean, double hyperinvshape)    {
        double alpha = 1.0 / hyperinvshape;
        double beta = alpha / hypermean;
        array.SetShape(alpha);
        array.SetScale(beta);
    }

    void Update() override {
        meansyn_bidimarray->Update();
        meanom_bidimarray->Update();
        if (syn_devmode)    {
            syn_bidimarray->SetInvShape(syn_invshape);
        }
        if (om_devmode) {
            om_bidimarray->SetInvShape(om_invshape);
        }
        Update(*branchsyn_array, branchsyn_hypermean, branchsyn_hyperinvshape);
        Update(*genesyn_array, genesyn_hypermean, genesyn_hyperinvshape);
        Update(*branchom_array, branchom_hypermean, branchom_hyperinvshape);
        Update(*geneom_array, geneom_hypermean, geneom_hyperinvshape);
    }

    void NoUpdate() {}

    // Log Priors

    double GetLogPrior() const {
        double total = 0;

        total += BranchSynHyperLogPrior();
        total += BranchSynLogPrior();
        total += GeneSynHyperLogPrior();
        total += GeneSynLogPrior();
        if (syn_devmode)    {
            total += SynInvShapeLogPrior();
            total += SynLogPrior();
        }

        total += BranchOmegaHyperLogPrior();
        total += BranchOmegaLogPrior();
        total += GeneOmegaHyperLogPrior();
        total += GeneOmegaLogPrior();
        if (om_devmode) {
            total += OmegaInvShapeLogPrior();
            total += OmegaLogPrior();
        }

        return total;
    }

    double BranchSynHyperLogPrior() const { return -branchsyn_hypermean - branchsyn_hyperinvshape; }

    double BranchSynLogPrior() const { 
        return branchsyn_array->GetLogProb(); 
    }

    double GeneSynHyperLogPrior() const { return -genesyn_hypermean - genesyn_hyperinvshape; }

    double GeneSynLogPrior() const {
        return genesyn_array->GetLogProb(); 
    }

    double SynInvShapeLogPrior() const { return -syn_invshape; }

    double SynLogPrior() const  {
        return syn_bidimarray->GetLogProb();
    }

    double BranchOmegaHyperLogPrior() const { return -branchom_hypermean - branchom_hyperinvshape; }

    double BranchOmegaLogPrior() const { 
        return branchom_array->GetLogProb(); 
    }

    double GeneOmegaHyperLogPrior() const { return -geneom_hypermean - geneom_hyperinvshape; }

    double GeneOmegaLogPrior() const { 
        return geneom_array->GetLogProb(); 
    }

    double OmegaInvShapeLogPrior() const { return -om_invshape; }

    double OmegaLogPrior() const {
        return om_bidimarray->GetLogProb();
    }

    // likelihood
    
    double GetLogLikelihood() const {

        double tot = 0;
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                tot += GetLogLikelihood(i,j);
            }
        }
        return tot;
    }

    double GetLogLikelihood(int gene, int branch) const {

        const dSOmegaPathSuffStat &suffstat = dsomss->GetVal(gene).GetVal(branch);
        double syncount = suffstat.GetSynCount();
        double synbeta = suffstat.GetSynBeta();
        double nonsyncount = suffstat.GetNonSynCount();
        double nonsynbeta = suffstat.GetNonSynBeta();

        double syn = syn_devmode ? syn_bidimarray->GetVal(gene).GetVal(branch):
            meansyn_bidimarray->GetVal(gene).GetVal(branch);

        double om = om_devmode ? om_bidimarray->GetVal(gene).GetVal(branch):
            meanom_bidimarray->GetVal(gene).GetVal(branch);

        return syncount * log(syn) - synbeta * syn 
            + nonsyncount * log(syn*om) - nonsynbeta * (syn*om);

    }

    double GetLogProb() const override  {
        return GetLogPrior() + GetLogLikelihood();
    }

    // suff stat log probs

    double HyperSuffStatLogProb(const GammaSuffStat& suffstat, 
            double mean, double invshape) const    {
        double alpha = 1.0 / invshape;
        double beta = alpha / mean;
        return suffstat.GetLogProb(alpha, beta);
    }

    // working marginally on l_ij, conditionally on omega_ij

    double SynSuffStatLogProb() const {
        double total = 0;
        for (int i=0; i<Ngene; i++) {
            total += GeneSynSuffStatLogProb(i);
        }
        return total;
    }

    double GeneSynSuffStatLogProb(int gene) const {
        double total = 0;
        for (int j=0; j<Nbranch; j++) {
            total += GeneBranchSynSuffStatLogProb(gene, j);
        }
        return total;
    }

    double BranchSynSuffStatLogProb(int branch) const {
        double total = 0;
        for (int i=0; i<Ngene; i++) {
            total += GeneBranchSynSuffStatLogProb(i, branch);
        }
        return total;
    }

    double GeneBranchSynSuffStatLogProb(int gene, int branch) const {

        const dSOmegaPathSuffStat &suffstat = dsomss->GetVal(gene).GetVal(branch);
        double syncount = suffstat.GetSynCount();
        double synbeta = suffstat.GetSynBeta();
        double nonsyncount = suffstat.GetNonSynCount();
        double nonsynbeta = suffstat.GetNonSynBeta();

        double om = om_devmode ? om_bidimarray->GetVal(gene).GetVal(branch):
            meanom_bidimarray->GetVal(gene).GetVal(branch);

        double synmean = genesyn_array->GetVal(gene) * branchsyn_array->GetVal(branch);

        if (! syn_devmode)  {
            return (syncount + nonsyncount) * log(synmean) 
                - (synbeta + om * nonsynbeta) * synmean;
        }

        double alpha = 1.0 / syn_invshape;
        double beta = alpha / synmean;

        double postalpha = alpha + syncount + nonsyncount;
        double postbeta = beta + synbeta + om * nonsynbeta;

        return alpha * log(beta) - Random::logGamma(alpha) 
            + Random::logGamma(postalpha) - postalpha * log(postbeta);
    }

    void ResampleSyn(int gene, int branch)  {

        const dSOmegaPathSuffStat &suffstat = dsomss->GetVal(gene).GetVal(branch);
        double syncount = suffstat.GetSynCount();
        double synbeta = suffstat.GetSynBeta();
        double nonsyncount = suffstat.GetNonSynCount();
        double nonsynbeta = suffstat.GetNonSynBeta();

        double om = om_bidimarray->GetVal(gene).GetVal(branch);
        double synmean = genesyn_array->GetVal(gene) * branchsyn_array->GetVal(branch);

        double alpha = 1.0 / syn_invshape;
        double beta = alpha / synmean;

        double postalpha = alpha + syncount + nonsyncount;
        double postbeta = beta + synbeta + om * nonsynbeta;

        (*syn_bidimarray)[gene][branch] = Random::Gamma(postalpha, postbeta);
    }

    void ResampleSyn()   {
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                ResampleSyn(i,j);
            }
        }
    }

    // working marginally on omega_ij, conditionally on l_ij

    double OmegaSuffStatLogProb() const {
        double total = 0;
        for (int i=0; i<Ngene; i++) {
            total += GeneOmegaSuffStatLogProb(i);
        }
        return total;
    }

    double GeneOmegaSuffStatLogProb(int gene) const {
        double total = 0;
        for (int j=0; j<Nbranch; j++) {
            total += GeneBranchOmegaSuffStatLogProb(gene, j);
        }
        return total;
    }

    double BranchOmegaSuffStatLogProb(int branch) const {
        double total = 0;
        for (int i=0; i<Ngene; i++) {
            total += GeneBranchOmegaSuffStatLogProb(i, branch);
        }
        return total;
    }

    double GeneBranchOmegaSuffStatLogProb(int gene, int branch) const {

        const dSOmegaPathSuffStat &suffstat = dsomss->GetVal(gene).GetVal(branch);
        double nonsyncount = suffstat.GetNonSynCount();
        double nonsynbeta = suffstat.GetNonSynBeta();

        double syn = syn_devmode ? syn_bidimarray->GetVal(gene).GetVal(branch):
            meansyn_bidimarray->GetVal(gene).GetVal(branch);

        double ommean = geneom_array->GetVal(gene) * branchom_array->GetVal(branch);

        if (! om_devmode)  {
            return nonsyncount * log(ommean) -  syn * nonsynbeta * ommean;
        }

        double alpha = 1.0 / om_invshape;
        double beta = alpha / ommean;

        double postalpha = alpha + nonsyncount;
        double postbeta = beta + syn * nonsynbeta;

        return alpha * log(beta) - Random::logGamma(alpha) 
            + Random::logGamma(postalpha) - postalpha * log(postbeta);
    }

    void ResampleOmega(int gene, int branch)    {

        const dSOmegaPathSuffStat &suffstat = dsomss->GetVal(gene).GetVal(branch);
        double nonsyncount = suffstat.GetNonSynCount();
        double nonsynbeta = suffstat.GetNonSynBeta();

        double syn = syn_bidimarray->GetVal(gene).GetVal(branch);
        double ommean = geneom_array->GetVal(gene) * branchom_array->GetVal(branch);

        double alpha = 1.0 / om_invshape;
        double beta = alpha / ommean;

        double postalpha = alpha + nonsyncount;
        double postbeta = beta + syn * nonsynbeta;

        (*om_bidimarray)[gene][branch] = Random::Gamma(postalpha, postbeta);
    }

    void ResampleOmega()    {
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                ResampleOmega(i,j);
            }
        }
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // log prob for moving branch omega hyperparameters
    double BranchSynHyperLogProb() const { 
        return BranchSynHyperLogPrior() + 
            HyperSuffStatLogProb(branchsyn_hypersuffstat, branchsyn_hypermean, branchsyn_hyperinvshape); 
    }

    // log prob for moving gene omega hyperparameters
    double GeneSynHyperLogProb() const { 
        return GeneSynHyperLogPrior() + 
            HyperSuffStatLogProb(genesyn_hypersuffstat, genesyn_hypermean, genesyn_hyperinvshape); 
    }

    // log prob for moving branch omega hyperparameters
    double BranchOmegaHyperLogProb() const { 
        return BranchOmegaHyperLogPrior() + 
            HyperSuffStatLogProb(branchom_hypersuffstat, branchom_hypermean, branchom_hyperinvshape); 
    }

    // log prob for moving gene omega hyperparameters
    double GeneOmegaHyperLogProb() const { 
        return GeneOmegaHyperLogPrior() + 
            HyperSuffStatLogProb(geneom_hypersuffstat, geneom_hypermean, geneom_hyperinvshape); 
    }

    double SynInvShapeLogProb() const {
        return SynInvShapeLogPrior() + SynSuffStatLogProb();
    }

    double OmegaInvShapeLogProb() const {
        return OmegaInvShapeLogPrior() + OmegaSuffStatLogProb();
    }

    double Move() override  {
        int nrep = 10;
        int nsmallrep = 1;

        for (int rep=0; rep<nrep; rep++) {

            for (int smallrep=0; smallrep<nsmallrep; smallrep++)    {

                MoveGeneArray(*genesyn_array, 1.0, 1, 
                        [this] (int gene) {return this->GeneSynSuffStatLogProb(gene);});
                MoveBranchArray(*branchsyn_array, 1.0, 1, 
                        [this] (int branch) {return this->BranchSynSuffStatLogProb(branch);});
                CompensatoryMove(*genesyn_array, *branchsyn_array, 1.0, 3);
                meansyn_bidimarray->Update();

                if (syn_devmode)    {
                    ScalingMove(syn_invshape, 0.3, 1, 
                            &FastGeneBranchOmegaModel::SynInvShapeLogProb,
                            &FastGeneBranchOmegaModel::NoUpdate,
                            this);
                    syn_bidimarray->SetInvShape(syn_invshape);
                }
            }

            if (syn_devmode)    {
                ResampleSyn();
            }

            for (int smallrep=0; smallrep<nsmallrep; smallrep++)    {

                MoveGeneArray(*geneom_array, 1.0, 1, 
                        [this] (int gene) {return this->GeneOmegaSuffStatLogProb(gene);});
                MoveBranchArray(*branchom_array, 1.0, 1, 
                        [this] (int branch) {return this->BranchOmegaSuffStatLogProb(branch);});
                CompensatoryMove(*geneom_array, *branchom_array, 1.0, 3);
                meanom_bidimarray->Update();

                if (om_devmode)    {
                    ScalingMove(om_invshape, 0.3, 1, 
                            &FastGeneBranchOmegaModel::OmegaInvShapeLogProb,
                            &FastGeneBranchOmegaModel::NoUpdate,
                            this);
                    om_bidimarray->SetInvShape(om_invshape);
                }
            }

            if (om_devmode) {
                ResampleOmega();
            }

            MoveHyperParam(branchsyn_hypersuffstat, *branchsyn_array,
                    branchsyn_hypermean, branchsyn_hypermean, branchsyn_hyperinvshape,
                    [this] () {return this->BranchSynHyperLogProb();},
                    [] () {},
                    1.0,100);

            MoveHyperParam(branchsyn_hypersuffstat, *branchsyn_array,
                    branchsyn_hyperinvshape, branchsyn_hypermean, branchsyn_hyperinvshape,
                    [this] () {return this->BranchSynHyperLogProb();},
                    [] () {},
                    1.0,100);

            MoveHyperParam(genesyn_hypersuffstat, *genesyn_array,
                    genesyn_hyperinvshape, genesyn_hypermean, genesyn_hyperinvshape,
                    [this] () {return this->GeneSynHyperLogProb();},
                    [] () {},
                    1.0,100);

            MoveHyperParam(branchom_hypersuffstat, *branchom_array,
                    branchom_hyperinvshape, branchom_hypermean, branchom_hyperinvshape,
                    [this] () {return this->BranchOmegaHyperLogProb();},
                    [] () {},
                    1.0,100);

            MoveHyperParam(geneom_hypersuffstat, *geneom_array,
                    geneom_hypermean, geneom_hypermean, geneom_hyperinvshape,
                    [this] () {return this->GeneOmegaHyperLogProb();},
                    [] () {},
                    1.0,100);

            MoveHyperParam(geneom_hypersuffstat, *geneom_array,
                    geneom_hyperinvshape, geneom_hypermean, geneom_hyperinvshape,
                    [this] () {return this->GeneOmegaHyperLogProb();},
                    [] () {},
                    1.0,100);

        }
        Update();
        return 1.0;
    }

    template <class LogProb>
    double MoveBranchArray(IIDGamma& array, double tuning, int nrep, LogProb logprob) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++) {
            for (int j=0; j<Nbranch; j++) {
                double deltalogprob = -array.GetLogProb(j);
                deltalogprob -= logprob(j);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                array[j] *= e;
                deltalogprob += array.GetLogProb(j);
                deltalogprob += logprob(j);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    array[j] /= e;
                }
            }
        }
        double ret = nacc / GetNbranch() / nrep;
        return ret;
    }

    template <class LogProb>
    double MoveGeneArray(IIDGamma& array, double tuning, int nrep, LogProb logprob)  {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++) {
            for (int i=0; i<Ngene; i++) {
                double deltalogprob = -array.GetLogProb(i);
                deltalogprob -= logprob(i);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                array[i] *= e;
                deltalogprob += array.GetLogProb(i);
                deltalogprob += logprob(i);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    array[i] /= e;
                }
            }
        }
        double ret = nacc / GetNgene() / nrep;
        return ret;
    }

    double CompensatoryMove(IIDGamma& genearray, IIDGamma& brancharray, double tuning, int nrep)    {
        double nacc = 0;
        if (genearray.GetSize() != Ngene) {
            cerr << "error in comp move: first argument should be gene array\n";
            exit(1);
        }
        if (brancharray.GetSize() != Nbranch)  {
            cerr << "error in comp move: second argument should be branch array\n";
            exit(1);
        }
        for (int rep=0; rep<nrep; rep++)    {

            double deltalogprob = - genearray.GetLogProb() - brancharray.GetLogProb();
            
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            for (int i=0; i<Ngene; i++) {
                genearray[i] *= e;
            }
            for (int j=0; j<Nbranch; j++)   {
                brancharray[j] /= e;
            }

            deltalogprob += genearray.GetLogProb() + brancharray.GetLogProb();
            deltalogprob += (Ngene - Nbranch) * m;

            int acc = (log(Random::Uniform()) < deltalogprob);
            if (acc) {
                nacc++;
            } else {
                for (int i=0; i<Ngene; i++) {
                    genearray[i] /= e;
                }
                for (int j=0; j<Nbranch; j++)   {
                    brancharray[j] *= e;
                }
            }
        }
        double ret = nacc / GetNgene() / nrep;
        return ret;
    }

    template <class LogProb, class UpdateF>
    double MoveHyperParam(GammaSuffStat& hypersuffstat, IIDGamma& array, 
            double& hyperparam, double& hypermean, double& hyperinvshape, 
            LogProb logprob, UpdateF update, double tuning, int nrep)    {

        hypersuffstat.Clear();
        hypersuffstat.AddSuffStat(array);

        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -logprob();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            hyperparam *= e;
            update();
            deltalogprob += logprob();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                hyperparam /= e;
                update();
            }
            ntot++;
        }

        Update(array, hypermean, hyperinvshape);

        return nacc / ntot;
    }
};

