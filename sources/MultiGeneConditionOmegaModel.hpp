
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

#include "ConditionOmegaModel.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"

/**
 * \brief An Array of BranchArray's of OmegaPathSuffStat
 *
 * used in MultiGeneConditionOmegaModel, where each gene has a BranchArray of
 * OmegaPathSuffStat.
 */

class OmegaPathSuffStatBidimArray : public Array<OmegaPathSuffStatArray> {
  public:
    //! constructor, parameterized by underlying tree and size (number of genes)
    OmegaPathSuffStatBidimArray(int inncond, int insize)
        : ncond(inncond), array(insize, (OmegaPathSuffStatArray *)0) {
        for (int i = 0; i < GetSize(); i++) {
            array[i] = new OmegaPathSuffStatArray(ncond);
        }
    }

    ~OmegaPathSuffStatBidimArray() {
        for (int i = 0; i < GetSize(); i++) {
            delete[] array[i];
        }
    }

    int GetSize() const override { return array.size(); }

    const OmegaPathSuffStatArray &GetVal(int i) const override { return *array[i]; }

    OmegaPathSuffStatArray &operator[](int i) override { return *array[i]; }

    //! clear all suff stats
    void Clear() {
        for (int i = 0; i < GetSize(); i++) {
            array[i]->Clear();
        }
    }

  private:
    int ncond;
    vector<OmegaPathSuffStatArray *> array;
};

class MultiGeneConditionOmegaModel : public MultiGeneProbModel {

  private:
    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;

    string treefile;

    int Ntaxa;
    int Nbranch;
    int Ncond;
    int Nlevel;

    int blmode;
    int nucmode;
    int devmode;

    int ppredmode;

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

    BranchAllocationSystem* branchalloc;

    double condvhypermean;
    double condvhyperinvshape;
    IIDGamma *condvarray;
    GammaSuffStat hypercondvsuffstat;

    double genewhypermean;
    double genewhyperinvshape;
    IIDGamma *genewarray;
    GammaSuffStat hypergenewsuffstat;

    double omegainvshape;

    ProductArray *meanomegabidimarray;
    ConditionSpecificMeanGammaBidimArray *condomegabidimarray;

    OmegaPathSuffStatBidimArray *omegapathsuffstatbidimarray;

    // each gene defines its own ConditionOmegaModel
    std::vector<ConditionOmegaModel *> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneConditionOmegaModel(string datafile, string intreefile, int inNcond, int inNlevel, 
                                 int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {

        blmode = 1;
        nucmode = 1;
        devmode = 1;

        ppredmode = 1;

        AllocateAlignments(datafile);
        treefile = intreefile;
        Ncond = inNcond;
        Nlevel = inNlevel;

        refcodondata = new CodonSequenceAlignment(refdata, true);
        taxonset = refdata->GetTaxonSet();
        Ntaxa = refdata->GetNtaxa();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        if (! Ncond)    {
            Ncond = Nbranch;
        }

        branchalloc = new BranchAllocationSystem(*tree, Ncond);

        if (!myid) {
            std::cerr << "number of taxa : " << Ntaxa << '\n';
            std::cerr << "number of branches : " << Nbranch << '\n';
            std::cerr << "-- Tree and data fit together\n";
        }
    }

    void Allocate() {

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
        blhyperinvshape = 0.1;
        if (blmode == 2) {
            lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
            lengthhypersuffstatarray = 0;
        } else {
            branchlength->SetAllBranches(1.0 / lambda);
            branchlengtharray =
                new GammaWhiteNoiseArray(GetLocalNgene(), *tree, *branchlength, 1.0 / blhyperinvshape);
            lengthpathsuffstatarray = 0;
            lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
        }

        // Nucleotide rates

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;

        if (nucmode == 2) {
            nucrelratearray = new IIDDirichlet(1, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray = new IIDDirichlet(1, nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = new GTRSubMatrix(Nnuc, (*nucrelratearray)[0], (*nucstatarray)[0], true);
        } else {
            nucrelratearray =
                new IIDDirichlet(GetLocalNgene(), nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
            nucstatarray =
                new IIDDirichlet(GetLocalNgene(), nucstathypercenter, 1.0 / nucstathyperinvconc);
            nucmatrix = 0;
        }

        condvhypermean = 1.0;
        condvhyperinvshape = 1.0;
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condvarray = new IIDGamma(Ncond, alpha, beta);
        for (int j = 0; j < Ncond; j++) {
            (*condvarray)[j] = 1.0;
        }

        genewhypermean = 1.0;
        genewhyperinvshape = 1.0;
        double genealpha = 1.0 / genewhyperinvshape;
        double genebeta = genealpha / genewhypermean;
        genewarray = new IIDGamma(GetLocalNgene(), genealpha, genebeta);

        meanomegabidimarray = new ProductArray(*condvarray, *genewarray);
        if (devmode)    {
            omegainvshape = 1.0;
        }
        else    {
            omegainvshape = 0;
        }
        condomegabidimarray =
            new ConditionSpecificMeanGammaBidimArray(*meanomegabidimarray, omegainvshape);

        omegapathsuffstatbidimarray = new OmegaPathSuffStatBidimArray(Ncond, GetLocalNgene());

        lnL = 0;
        GeneLogPrior = 0;

        if (!GetMyid()) {
            geneprocess.assign(0, (ConditionOmegaModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (ConditionOmegaModel *)0);

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene] =
                    new ConditionOmegaModel(GetLocalGeneName(gene), treefile, Ncond, Nlevel);
                geneprocess[gene]->SetAcrossGenesModes(blmode, nucmode);
                geneprocess[gene]->Allocate();
            }
        }
    }

    // called upon constructing the model
    // mode == 2: global 
    // mode == 1: gene specific, with hyperparameters estimated across genes
    // mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode)   {
        blmode = inblmode;
        nucmode = innucmode;
    }

    void SetDeviationMode(int indevmode)    {
        devmode = indevmode;
    }

    void SetPostPredMode(int mode)  {
        ppredmode = mode;
    }

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    int GetNbranch() const { return tree->GetNbranch(); }

    const Tree *GetTree() const { return tree; }

    int GetNcond() const { return Ncond; }

    //-------------------
    // Traces and Monitors
    //-------------------

    void PrintBranchIndices(ostream& os) const {
        RecursivePrintBranchIndices(os,GetTree()->GetRoot());
    } 

    void RecursivePrintBranchIndices(ostream& os, const Link* from) const {
        if (!from->isLeaf()) {
            os << '(';
            for (const Link *link = from->Next(); link != from; link = link->Next()) {
                RecursivePrintBranchIndices(os, link->Out());
                if (link->Next() != from) {
                    os << ',';
                }
            }
            os << ')';
        }
        os << from->GetNode()->GetName();
        if (!from->isRoot()) {
            os << ':' << branchalloc->GetBranchAlloc(from->GetBranch()->GetIndex());
        }
    }

    void TraceHeader(ostream &os) const {
        os << "#logprior\tlnL";
        if (blmode == 2) {
            os << "\tlength";
        } else {
            os << "\tmeanlength\tstdev";
        }
        os << "\tgenemean\tinvshape";
        os << "\tcondinvshape";
        os << "\tomegainvshape";
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
        os << '\t' << genewhypermean << '\t' << genewhyperinvshape;
        os << '\t' << condvhyperinvshape;
        os << '\t' << omegainvshape;

        os << '\t' << nucstatarray->GetMeanEntropy();
        os << '\t' << nucrelratearray->GetMeanEntropy();
        if (nucmode != 2) {
            os << '\t' << sqrt(GetVarNucRelRate()) << '\t' << Random::GetEntropy(nucrelratehypercenter)
               << '\t' << nucrelratehyperinvconc;
            os << '\t' << sqrt(GetVarNucStat()) << '\t' << Random::GetEntropy(nucstathypercenter)
               << '\t' << nucstathyperinvconc;
        }
        os << '\n';
        os.flush();
    }

    void PrintGeneEffects(ostream &os) const {
        for (int i = 0; i < GetNgene(); i++) {
            os << genewarray->GetVal(i) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void PrintCondEffects(ostream &os) const {
        for (int j = 0; j < GetNcond(); j++) {
            os << condvarray->GetVal(j) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void PrintDeviations(ostream &os) const {
        for (int j = 0; j < GetNcond(); j++) {
            for (int i = 0; i < GetNgene(); i++) {
                os << log(GetOmega(i, j) / GetMeanOmega(i, j)) << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

    void NoDeviations() {
        meanomegabidimarray->Update();
        for (int j = 0; j < GetNcond(); j++) {
            for (int i = 0; i < GetNgene(); i++) {
                (*condomegabidimarray)[i][j] = meanomegabidimarray->GetVal(i).GetVal(j);
            }
        }
    }

    // Branch lengths

    double GetMeanTotalLength() const {
        double tot = 0;
        for (int j = 0; j < Nbranch; j++) {
            tot += branchlength->GetVal(j);
        }
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
            for (int g = 0; g < Ngene; g++) {
                double tmp = (*nucrelratearray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= Ngene;
            var /= Ngene;
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
            for (int g = 0; g < Ngene; g++) {
                double tmp = (*nucstatarray)[g][j];
                mean += tmp;
                var += tmp * tmp;
            }
            mean /= Ngene;
            var /= Ngene;
            var -= mean * mean;
            tot += var;
        }
        tot /= Nnuc;
        return tot;
    }

    void Monitor(ostream &os) const {}

    void MasterToStream(ostream &os) const {

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

        os << condvhypermean << '\t' << condvhyperinvshape << '\t';
        os << *condvarray << '\t';
        os << genewhypermean << '\t' << genewhyperinvshape << '\t';
        os << *genewarray << '\t';
        os << omegainvshape << '\t';
        os << *condomegabidimarray << '\n';
    }

    void MasterFromStream(istream &is) {

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

        is >> condvhypermean >> condvhyperinvshape;
        is >> *condvarray;
        is >> genewhypermean >> genewhyperinvshape;
        is >> *genewarray;
        is >> omegainvshape;
        is >> *condomegabidimarray;
    }

    //-------------------
    // Updates
    //-------------------

    void FastUpdate() {

        branchlength->SetScale(lambda);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);

        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condvarray->SetShape(alpha);
        condvarray->SetScale(beta);
        double genealpha = 1.0 / genewhyperinvshape;
        double genebeta = genealpha / genewhypermean;
        genewarray->SetShape(genealpha);
        genewarray->SetScale(genebeta);
        meanomegabidimarray->Update();
        condomegabidimarray->SetInvShape(omegainvshape);
    }

    void MasterUpdate() override {
        FastUpdate();

        if (nprocs > 1) {
            MasterSendBranchLengthsHyperParameters();
            MasterSendNucRatesHyperParameters();

            if (blmode == 2) {
                MasterSendGlobalBranchLengths();
            } else {
                MasterSendGeneBranchLengths();
            }

            if (nucmode == 2) {
                MasterSendGlobalNucRates();
            } else {
                MasterSendGeneNucRates();
            }

            MasterSendOmegaHyperParameters();
            MasterSendOmega();
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {

        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveNucRatesHyperParameters();

        if (blmode == 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveGeneBranchLengths();
        }
        if (nucmode == 2) {
            SlaveReceiveGlobalNucRates();
        } else {
            SlaveReceiveGeneNucRates();
        }

        SlaveReceiveOmegaHyperParameters();
        SlaveReceiveOmega();
        GeneUpdate();
        SlaveSendLogProbs();
    }

    void GeneUpdate() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->Update();
        }
    }

    void MasterPostPred(string name) override {
        FastUpdate();

        if (! ppredmode)    {
            NoDeviations();
        }

        if (nprocs > 1) {
            MasterSendBranchLengthsHyperParameters();
            MasterSendNucRatesHyperParameters();

            if (blmode == 2) {
                MasterSendGlobalBranchLengths();
            } else {
                MasterSendGeneBranchLengths();
            }

            if (nucmode == 2) {
                MasterSendGlobalNucRates();
            } else {
                MasterSendGeneNucRates();
            }

            MasterSendOmegaHyperParameters();
            MasterSendOmega();
        }
    }

    void SlavePostPred(string name) override {

        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveNucRatesHyperParameters();

        if (blmode == 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveGeneBranchLengths();
        }

        if (nucmode == 2) {
            SlaveReceiveGlobalNucRates();
        } else {
            SlaveReceiveGeneNucRates();
        }

        SlaveReceiveOmegaHyperParameters();
        SlaveReceiveOmega();

        GenePostPred(name);
    }

    void GenePostPred(string name) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->PostPred(name + GetLocalGeneName(gene));
        }
    }

    void TouchNucMatrix() {
        nucmatrix->CopyStationary((*nucstatarray)[0]);
        nucmatrix->CorruptMatrix();
    }

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetOmega(int gene, int cond) const {
        return condomegabidimarray->GetVal(gene).GetVal(cond);
    }

    double GetMeanOmega(int gene, int cond) const {
        double tmp1 = meanomegabidimarray->GetVal(gene).GetVal(cond);
        double tmp2 = condvarray->GetVal(cond) * genewarray->GetVal(gene);
        if (fabs(tmp1 - tmp2) > 1e-6)   {
            cerr << "error in GetMeanOmega\n";
            cerr << tmp1 << '\t' << tmp2 << '\n';
            exit(1);
        }
        return meanomegabidimarray->GetVal(gene).GetVal(cond);
    }

    double GetLogPrior() const {
        double total = 0;
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

        total += CondVHyperLogPrior();
        total += CondVLogPrior();
        total += GeneWHyperLogPrior();
        total += GeneWLogPrior();
        if (devmode)    {
            total += OmegaInvShapeLogPrior();
            total += OmegaLogPrior();
        }
        return total;
    }

    double OmegaLogPrior() const {
        meanomegabidimarray->Update();
        condomegabidimarray->SetInvShape(omegainvshape);
        return condomegabidimarray->GetLogProb();
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

    double CondVHyperLogPrior() const { return -condvhypermean - condvhyperinvshape; }

    double GeneWHyperLogPrior() const { return -genewhypermean - genewhyperinvshape; }

    double CondVLogPrior() const { return condvarray->GetLogProb(); }

    double GeneWLogPrior() const { return genewarray->GetLogProb(); }

    double OmegaInvShapeLogPrior() const { return omegainvshape; }

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

    double CondVHyperSuffStatLogProb() const {
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        return hypercondvsuffstat.GetLogProb(alpha, beta);
    }

    double GeneWHyperSuffStatLogProb() const {
        double alpha = 1.0 / genewhyperinvshape;
        double beta = alpha / genewhypermean;
        return hypergenewsuffstat.GetLogProb(alpha, beta);
    }

    double OmegaSuffStatLogProb() const {
        double total = 0;
        for (int i = 0; i < GetNgene(); i++) {
            total += GeneOmegaSuffStatLogProb(i);
        }
        return total;
    }

    double GeneOmegaSuffStatLogProb(int gene) const {
        double total = 0;
        for (int j = 0; j < GetNcond(); j++) {
            total += GeneConditionOmegaSuffStatLogProb(gene, j);
        }
        return total;
    }

    double ConditionOmegaSuffStatLogProb(int cond) const {
        double total = 0;
        for (int i = 0; i < GetNgene(); i++) {
            total += GeneConditionOmegaSuffStatLogProb(i, cond);
        }
        return total;
    }

    double GeneConditionOmegaSuffStatLogProb(int gene, int cond) const {
        const OmegaPathSuffStat &suffstat = omegapathsuffstatbidimarray->GetVal(gene).GetVal(cond);
        int count = suffstat.GetCount();
        double b = suffstat.GetBeta();

        double mean = genewarray->GetVal(gene) * condvarray->GetVal(cond);

        if (!devmode)    {
            return count*log(mean) - b*mean;
        }

        double alpha = 1.0 / omegainvshape;
        double beta = alpha / mean;

        return alpha * log(beta) - Random::logGamma(alpha) + Random::logGamma(alpha + count) -
               (alpha + count) * log(beta + b);
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

    // log prob for moving branch v hyperparameters
    double CondVHyperLogProb() const { return CondVHyperLogPrior() + CondVHyperSuffStatLogProb(); }

    // log prob for moving gene w hyperparameters
    double GeneWHyperLogProb() const { return GeneWHyperLogPrior() + GeneWHyperSuffStatLogProb(); }

    //-------------------
    // Moves
    //-------------------

    // all methods starting with Master are called only by master
    // for each such method, there is a corresponding method called by slave, and
    // starting with Slave
    //
    // all methods starting with Gene are called only be slaves, and do some work
    // across all genes allocated to that slave

    void MasterMove() override {
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            MasterReceiveOmegaSuffStat();
            MoveOmegaHyperParameters(3);
            MasterSendOmegaHyperParameters();

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                MasterReceiveBranchLengthsSuffStat();
                ResampleBranchLengths();
                MoveLambda();
                MasterSendGlobalBranchLengths();
            } else if (blmode == 1) {
                MasterReceiveBranchLengthsHyperSuffStat();
                MoveBranchLengthsHyperParameters();
                MasterSendBranchLengthsHyperParameters();
            }

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                MasterReceiveNucPathSuffStat();
                MoveNucRates();
                MasterSendGlobalNucRates();
            } else if (nucmode == 1) {
                MasterReceiveNucRatesHyperSuffStat();
                MoveNucRatesHyperParameters();
                MasterSendNucRatesHyperParameters();
            }
        }

        // collect current state
        if (blmode != 2) {
            MasterReceiveGeneBranchLengths();
        }
        if (nucmode != 2) {
            MasterReceiveGeneNucRates();
        }
        MasterReceiveOmega();
        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {
        GeneResampleSub(1.0);

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {

            GeneCollectPathSuffStat();
            SlaveSendOmegaSuffStat();
            SlaveReceiveOmegaHyperParameters();
            GeneResampleOmega();

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                SlaveSendBranchLengthsSuffStat();
                SlaveReceiveGlobalBranchLengths();
            } else if (blmode == 1) {
                MoveGeneBranchLengths();
                SlaveSendBranchLengthsHyperSuffStat();
                SlaveReceiveBranchLengthsHyperParameters();
            }

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                SlaveSendNucPathSuffStat();
                SlaveReceiveGlobalNucRates();
            } else if (nucmode == 1) {
                MoveGeneNucRates();
                SlaveSendNucRatesHyperSuffStat();
                SlaveReceiveNucRatesHyperParameters();
            }
        }

        // collect current state
        if (blmode != 2) {
            SlaveSendGeneBranchLengths();
        }
        if (nucmode != 2) {
            SlaveSendGeneNucRates();
        }

        SlaveSendOmega();
        SlaveSendLogProbs();
    }

    void GeneResampleSub(double frac) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void GeneCollectPathSuffStat() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectPathSuffStat();
        }
    }

    void GeneResampleOmega() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleOmega();
        }
    }

    void MoveGeneNucRates() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveNucRates();
            geneprocess[gene]->GetNucRates((*nucrelratearray)[gene],(*nucstatarray)[gene]);
        }
    }

    void MoveGeneBranchLengths() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveBranchLengths();
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
    }

    // Branch lengths

    void ResampleBranchLengths() {
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneConditionOmegaModel::LambdaHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneConditionOmegaModel::LambdaHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {

        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);

        ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneConditionOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneConditionOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);

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
        ProfileMove(nucrelratehypercenter, 1.0, 1, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10, &MultiGeneConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveNucRates() {
        vector<double> &nucrelrate = (*nucrelratearray)[0];
        ProfileMove(nucrelrate, 0.1, 1, 10, &MultiGeneConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &MultiGeneConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &MultiGeneConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneConditionOmegaModel::TouchNucMatrix, this);

        vector<double> &nucstat = (*nucstatarray)[0];
        ProfileMove(nucstat, 0.1, 1, 10, &MultiGeneConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 10, &MultiGeneConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneConditionOmegaModel::TouchNucMatrix, this);
    }

    void MoveOmegaHyperParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveGeneW(1.0, 1);
            MoveCondV(1.0, 1);
            MoveGeneW(0.3, 1);
            MoveCondV(0.3, 1);
            if (devmode)    {
                MoveOmegaInvShape(0.3, 1);
            }
        }
        // MoveCondVHyperParams(1.0,100);
        MoveGeneWHyperParams(1.0, 100);
        meanomegabidimarray->Update();
    }

    double MoveCondV(double tuning, int nrep) {
        double nacc = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int j = 1; j < GetNcond(); j++) {
                double deltalogprob = -condvarray->GetLogProb(j);
                deltalogprob -= ConditionOmegaSuffStatLogProb(j);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*condvarray)[j] *= e;
                deltalogprob += condvarray->GetLogProb(j);
                deltalogprob += ConditionOmegaSuffStatLogProb(j);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*condvarray)[j] /= e;
                }
            }
        }
        return nacc / GetNbranch() / nrep;
    }

    double MoveGeneW(double tuning, int nrep) {
        double nacc = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < GetNgene(); i++) {
                double deltalogprior = -genewarray->GetLogProb(i);
                double deltasuffstatlogprob = -GeneOmegaSuffStatLogProb(i);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*genewarray)[i] *= e;
                deltalogprior += genewarray->GetLogProb(i);
                deltasuffstatlogprob += GeneOmegaSuffStatLogProb(i);
                double deltalogprob = deltalogprior + deltasuffstatlogprob + m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*genewarray)[i] /= e;
                }
            }
        }
        return nacc / GetNgene() / nrep;
    }

    double MoveOmegaInvShape(double tuning, int nrep) {
        double nacc = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -OmegaInvShapeLogPrior() - OmegaSuffStatLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            omegainvshape *= e;
            deltalogprob += OmegaInvShapeLogPrior() + OmegaSuffStatLogProb();
            deltalogprob += m;
            int acc = (log(Random::Uniform()) < deltalogprob);
            if (acc) {
                nacc++;
            } else {
                omegainvshape /= e;
            }
        }
        condomegabidimarray->SetInvShape(omegainvshape);
        return nacc / nrep;
    }

    void MoveCondVHyperParams(double tuning, int nrep) {
        hypercondvsuffstat.Clear();
        hypercondvsuffstat.AddSuffStat(*condvarray);

        ScalingMove(condvhyperinvshape, 1.0, 10, &MultiGeneConditionOmegaModel::CondVHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(condvhyperinvshape, 0.3, 10, &MultiGeneConditionOmegaModel::CondVHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);

        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condvarray->SetShape(alpha);
        condvarray->SetScale(beta);
    }

    void MoveGeneWHyperParams(double tuning, int nrep) {
        hypergenewsuffstat.Clear();
        hypergenewsuffstat.AddSuffStat(*genewarray);

        ScalingMove(genewhypermean, 1.0, 10, &MultiGeneConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(genewhypermean, 0.3, 10, &MultiGeneConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(genewhyperinvshape, 1.0, 10, &MultiGeneConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);
        ScalingMove(genewhyperinvshape, 0.3, 10, &MultiGeneConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneConditionOmegaModel::NoUpdate, this);

        double alpha = 1.0 / genewhyperinvshape;
        double beta = alpha / genewhypermean;
        genewarray->SetShape(alpha);
        genewarray->SetScale(beta);
    }

    //-------------------
    // MPI send / receive
    //-------------------

    // Branch lengths

    void MasterSendGlobalBranchLengths() { MasterSendGlobal(*branchlength); }

    void SlaveReceiveGlobalBranchLengths() {
        SlaveReceiveGlobal(*branchlength);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengths(*branchlength);
        }
    }

    void MasterSendBranchLengthsHyperParameters() {
        MasterSendGlobal(*branchlength, blhyperinvshape);
    }

    void SlaveReceiveBranchLengthsHyperParameters() {
        SlaveReceiveGlobal(*branchlength, blhyperinvshape);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
        }
    }

    void MasterSendGeneBranchLengths() {
        MasterSendGeneArray(*branchlengtharray);
    }

    void SlaveReceiveGeneBranchLengths() {
        SlaveReceiveGeneArray(*branchlengtharray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengths(branchlengtharray->GetVal(gene));
        }
    }

    void SlaveSendGeneBranchLengths() {
        // in principle, redundant..
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
        SlaveSendGeneArray(*branchlengtharray);
    }

    void MasterReceiveGeneBranchLengths() {
        MasterReceiveGeneArray(*branchlengtharray);
    }

    void SlaveSendBranchLengthsSuffStat() {
        lengthpathsuffstatarray->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectLengthSuffStat();
            lengthpathsuffstatarray->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
        }
        SlaveSendAdditive(*lengthpathsuffstatarray);
    }

    void MasterReceiveBranchLengthsSuffStat() {
        lengthpathsuffstatarray->Clear();
        MasterReceiveAdditive(*lengthpathsuffstatarray);
    }

    void SlaveSendBranchLengthsHyperSuffStat() {
        lengthhypersuffstatarray->Clear();
        lengthhypersuffstatarray->AddSuffStat(*branchlengtharray);
        SlaveSendAdditive(*lengthhypersuffstatarray);
    }

    void MasterReceiveBranchLengthsHyperSuffStat() {
        lengthhypersuffstatarray->Clear();
        MasterReceiveAdditive(*lengthhypersuffstatarray);
    }

    // Nucleotide Rates

    void MasterSendGlobalNucRates() {
        MasterSendGlobal(nucrelratearray->GetVal(0), nucstatarray->GetVal(0));
    }

    void SlaveReceiveGlobalNucRates() {
        SlaveReceiveGlobal((*nucrelratearray)[0], (*nucstatarray)[0]);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[0], (*nucstatarray)[0]);
        }
    }

    void MasterSendGeneNucRates() {
        MasterSendGeneArray(*nucrelratearray, *nucstatarray);
    }

    void SlaveReceiveGeneNucRates() {
        SlaveReceiveGeneArray(*nucrelratearray, *nucstatarray);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
        }
    }

    void SlaveSendGeneNucRates() {
        SlaveSendGeneArray(*nucrelratearray, *nucstatarray);
    }

    void MasterReceiveGeneNucRates() {
        MasterReceiveGeneArray(*nucrelratearray, *nucstatarray);
    }

    void MasterSendNucRatesHyperParameters() {
        MasterSendGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        MasterSendGlobal(nucstathypercenter, nucstathyperinvconc);
    }

    void SlaveReceiveNucRatesHyperParameters() {
        SlaveReceiveGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        SlaveReceiveGlobal(nucstathypercenter, nucstathyperinvconc);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter, nucrelratehyperinvconc,
                                                          nucstathypercenter, nucstathyperinvconc);
        }
    }

    void SlaveSendNucRatesHyperSuffStat() {
        nucrelratesuffstat.Clear();
        nucrelratearray->AddSuffStat(nucrelratesuffstat);
        SlaveSendAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        nucstatarray->AddSuffStat(nucstatsuffstat);
        SlaveSendAdditive(nucstatsuffstat);
    }

    void MasterReceiveNucRatesHyperSuffStat() {
        nucrelratesuffstat.Clear();
        MasterReceiveAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        MasterReceiveAdditive(nucstatsuffstat);
    }

    void SlaveSendNucPathSuffStat() {
        nucpathsuffstat.Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
        }

        SlaveSendAdditive(nucpathsuffstat);
    }

    void MasterReceiveNucPathSuffStat() {
        nucpathsuffstat.Clear();
        MasterReceiveAdditive(nucpathsuffstat);
    }

    // omega arrays

    void SlaveSendOmega() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            (*condomegabidimarray)[gene].Copy(*geneprocess[gene]->GetOmegaArray());
        }
        SlaveSendGeneArray(*condomegabidimarray);
    }

    void MasterReceiveOmega() { MasterReceiveGeneArray(*condomegabidimarray); }

    void MasterSendOmega() { MasterSendGeneArray(*condomegabidimarray); }

    void SlaveReceiveOmega() {
        SlaveReceiveGeneArray(*condomegabidimarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmegaTree(condomegabidimarray->GetVal(gene));
        }
    }

    // omega hyperparameters

    void MasterSendOmegaHyperParameters() {
        MasterSendGlobal(condvhypermean, condvhyperinvshape);
        MasterSendGlobal(genewhypermean, genewhyperinvshape);
        MasterSendGlobal(omegainvshape);
        MasterSendGlobal(*condvarray);
        MasterSendGeneArray(*genewarray);
    }

    void SlaveReceiveOmegaHyperParameters() {
        SlaveReceiveGlobal(condvhypermean, condvhyperinvshape);
        SlaveReceiveGlobal(genewhypermean, genewhyperinvshape);
        SlaveReceiveGlobal(omegainvshape);
        SlaveReceiveGlobal(*condvarray);
        SlaveReceiveGeneArray(*genewarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetCondVHyperParams(condvhypermean, condvhyperinvshape);
            geneprocess[gene]->SetCondV(*condvarray);
            geneprocess[gene]->SetGeneW(genewarray->GetVal(gene));
            geneprocess[gene]->SetOmegaHyperInvShape(omegainvshape);
        }
    }

    // omega path suff stat

    void SlaveSendOmegaSuffStat() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectOmegaSuffStat();
            (*omegapathsuffstatbidimarray)[gene].Array<OmegaPathSuffStat>::Copy(
                *geneprocess[gene]->GetOmegaPathSuffStatArray());
        }
        SlaveSendGeneArray(*omegapathsuffstatbidimarray);
    }

    void MasterReceiveOmegaSuffStat() { MasterReceiveGeneArray(*omegapathsuffstatbidimarray); }

    // log probs

    void SlaveSendLogProbs() {
        GeneLogPrior = 0;
        lnL = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
    }

    void MasterReceiveLogProbs() {
        GeneLogPrior = 0;
        MasterReceiveAdditive(GeneLogPrior);
        lnL = 0;
        MasterReceiveAdditive(lnL);
    }
};
