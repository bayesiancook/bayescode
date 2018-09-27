
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

#include "SparseConditionOmegaModel.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"
#include "GeneIIDMultiGamma.hpp"
#include "GeneIIDMultiDiscrete.hpp"

/**
 * \brief An Array of BranchArray's of OmegaPathSuffStat
 *
 * used in MultiGeneSparseConditionOmegaModel, where each gene has a BranchArray of
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

class MultiGeneSparseConditionOmegaModel : public MultiGeneProbModel {

  private:
    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;

    string treefile;

    int Ntaxa;
    int Nbranch;
    int Ncond;
    int Nlevel;

    int burnin;

    int blmode;
    int nucmode;

    double pineg;
    double pipos;

    int ppredmode;

    int modalprior;

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
    double condvhypervar;
    IIDNormal *condvarray;
    NormalSuffStat condvhypersuffstat;

    double genewhypermean;
    double genewhypervar;
    IIDNormal *genewarray;
    NormalSuffStat genewhypersuffstat;

    SumArray *meanlogomegabidimarray;

    vector<double> picenter;
    double piconcentration;
    IIDDirichlet* pi;
    IIDGamma* meanpos;
    IIDGamma* invshapepos;
    IIDGamma* meanneg;
    IIDGamma* invshapeneg;

    GeneIIDMultiDiscrete* alloc;
    GeneIIDMultiGamma* devpos;
    GeneIIDMultiGamma* devneg;

    MultiGammaSuffStat* devposhypersuffstat;
    MultiGammaSuffStat* devneghypersuffstat;

    GeneIIDLogNormalMixArray *condomegabidimarray;

    OmegaPathSuffStatBidimArray *omegapathsuffstatbidimarray;

    // each gene defines its own SparseConditionOmegaModel
    std::vector<SparseConditionOmegaModel *> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneSparseConditionOmegaModel(string datafile, string intreefile, int inNcond, int inNlevel, double inpipos, double inpineg,
                                 int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {

        blmode = 1;
        nucmode = 1;

        burnin = 0;

        ppredmode = 1;

        pipos = inpipos;
        pineg = inpineg;

        modalprior = 1;

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

        condvhypermean = 0;
        condvhypervar = 1.0;
        condvarray = new IIDNormal(Ncond, condvhypermean, condvhypervar);
        for (int j = 0; j < Ncond; j++) {
            (*condvarray)[j] = 0;
        }

        genewhypermean = -1;
        genewhypervar = 0.1;
        genewarray = new IIDNormal(GetLocalNgene(), genewhypermean, genewhypervar);
        /*
        for (int gene=0; gene<GetNgene(); gene++)   {
            (*genewarray)[gene] = 0;
        }
        */

        meanlogomegabidimarray = new SumArray(*condvarray, *genewarray);

        picenter.assign(3,0);
        if (pipos + pineg >= 1.0)   {
            cerr << "error: pi's sum up to more than 1.0\n";
            exit(1);
        }
        picenter[0] = pineg;
        picenter[2] = pipos;
        picenter[1] = 1.0 - pipos - pineg;
        piconcentration = 20;
        pi = new IIDDirichlet(Ncond,picenter,piconcentration);
        for (int cond=0; cond<Ncond; cond++)    {
            for (int k=0; k<3; k++) {
                (*pi)[cond][k] = picenter[k];
            }
        }
        /*
        if (! pineg)  {
            for (int cond=0; cond<Ncond; cond++)    {
                (*pi)[cond][0] = 0;
            }
        }
        if (! pipos)  {
            for (int cond=0; cond<Ncond; cond++)    {
                (*pi)[cond][2] = 0;
            }
        }
        for (int cond=0; cond<Ncond; cond++)    {
            (*pi)[cond][1] = 1.0 - (*pi)[cond][2] - (*pi)[cond][0];
        }
        */

        alloc = new GeneIIDMultiDiscrete(GetLocalNgene(),*pi);
        alloc->SetVal(1);

        meanpos = new IIDGamma(Ncond,1.0,1.0);
        invshapepos = new IIDGamma(Ncond,1.0,1.0);
        for (int k=0; k<Ncond; k++) {
            (*meanpos)[k] = 0.01;
            (*invshapepos)[k] = 0.5;
        }
        devpos = new GeneIIDMultiGamma(GetLocalNgene(),*meanpos,*invshapepos);
        devposhypersuffstat = new MultiGammaSuffStat(Ncond);

        meanneg = new IIDGamma(Ncond,1.0,1.0);
        invshapeneg = new IIDGamma(Ncond,1.0,1.0);
        for (int k=0; k<Ncond; k++) {
            (*meanneg)[k] = 0.01;
            (*invshapeneg)[k] = 0.5;
        }
        devneg = new GeneIIDMultiGamma(GetLocalNgene(),*meanneg,*invshapeneg);
        devneghypersuffstat = new MultiGammaSuffStat(Ncond);

        condomegabidimarray = new GeneIIDLogNormalMixArray(*meanlogomegabidimarray, *alloc, *devpos, *devneg);

        omegapathsuffstatbidimarray = new OmegaPathSuffStatBidimArray(Ncond, GetLocalNgene());

        lnL = 0;
        GeneLogPrior = 0;

        if (!GetMyid()) {
            geneprocess.assign(0, (SparseConditionOmegaModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (SparseConditionOmegaModel *)0);

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene] =
                    new SparseConditionOmegaModel(GetLocalGeneName(gene), treefile, Ncond, Nlevel);
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
        os << "\tgenemean\tvar";
        os << "\tcondmean\tcondvar";
        for (int cond=0; cond<Ncond; cond++)    {
            os << "\tpipos" << cond;
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << "\tpineg" << cond;
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << "\tmeanpos" << cond;
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << "\tmeanneg" << cond;
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << "\tinvshapepos" << cond;
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << "\tinvshapeneg" << cond;
        }
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
        os << '\t' << genewarray->GetEmpiricalMean() << '\t' << genewarray->GetEmpiricalVar();
        os << '\t' << condvarray->GetEmpiricalMean() << '\t' << condvarray->GetEmpiricalVar();
        for (int cond=0; cond<Ncond; cond++)    {
            os << '\t' << pi->GetVal(cond)[2];
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << '\t' << pi->GetVal(cond)[0];
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << '\t' << meanpos->GetVal(cond);
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << '\t' << meanneg->GetVal(cond);
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << '\t' << invshapepos->GetVal(cond);
        }
        for (int cond=0; cond<Ncond; cond++)    {
            os << '\t' << invshapeneg->GetVal(cond);
        }

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
            os << exp(genewarray->GetVal(i)) << '\t';
        }
        // os << *genewarray;
        os << '\n';
        os.flush();
    }

    void PrintCondEffects(ostream &os) const {
        for (int j = 0; j < GetNcond(); j++) {
            os << exp(condvarray->GetVal(j)) << '\t';
        }
        // os << *condvarray;
        os << '\n';
        os.flush();
    }

    void PrintDeviations(ostream &os) const {
        for (int j = 0; j < GetNcond(); j++) {
            for (int i = 0; i < GetNgene(); i++) {
                os << log(GetOmega(i, j)) - GetMeanLogOmega(i, j) << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

    void NoDeviations() {
        meanlogomegabidimarray->Update();
        for (int j = 0; j < GetNcond(); j++) {
            for (int i = 0; i < GetNgene(); i++) {
                (*condomegabidimarray)[i][j] = exp(meanlogomegabidimarray->GetVal(i).GetVal(j));
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

        os << condvhypermean << '\t' << condvhypervar << '\t';
        os << *condvarray << '\t';
        os << genewhypermean << '\t' << genewhypervar << '\t';
        os << *genewarray << '\t';
        os << *pi << '\t';
        os << *alloc << '\t';
        os << *meanpos << '\t' << *invshapepos << '\t';
        os << *devpos << '\t';
        os << *meanneg << '\t' << *invshapeneg << '\t';
        os << *devneg;
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

        is >> condvhypermean >> condvhypervar;
        is >> *condvarray;
        is >> genewhypermean >> genewhypervar;
        is >> *genewarray;
        is >> *pi >> *alloc;
        is >> *meanpos >> *invshapepos >> *devpos;
        is >> *meanneg >> *invshapeneg >> *devneg;
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

        condvarray->SetMean(condvhypermean);
        condvarray->SetVar(condvhypervar);
        genewarray->SetMean(genewhypermean);
        genewarray->SetVar(genewhypervar);

        meanlogomegabidimarray->Update();
        condomegabidimarray->Update();
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

            // MasterSendOmegaHyperParameters();
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

        // SlaveReceiveOmegaHyperParameters();
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

            // MasterSendOmegaHyperParameters();
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

        // SlaveReceiveOmegaHyperParameters();
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

    int GetAlloc(int gene, int cond) const  {
        return alloc->GetVal(gene).GetVal(cond);
    }

    double GetOmega(int gene, int cond) const {
        return condomegabidimarray->GetVal(gene).GetVal(cond);
    }

    double GetMeanLogOmega(int gene, int cond) const {
        /*
        double tmp1 = meanlogomegabidimarray->GetVal(gene).GetVal(cond);
        double tmp2 = genewarray->GetVal(gene) + condvarray->GetVal(cond);
        if (fabs(tmp1 - tmp2) > 1e-6)   {
            cerr << "error in get mean log omega\n";
            cerr << tmp1 << '\t' << tmp2 << '\n';
            exit(1);
        }
        return meanlogomegabidimarray->GetVal(gene).GetVal(cond);
        */
        return genewarray->GetVal(gene) + condvarray->GetVal(cond);
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

        if (pipos)  {
            total += DevPosHyperLogPrior();
            total += DevPosLogPrior();
        }
        if (pineg)  {
            total += DevNegHyperLogPrior();
            total += DevNegLogPrior();
        }

        // total += PiLogPrior();
        total += AllocLogPrior();

        if (std::isnan(total))   {
            cerr << "in GetLogPrior: nan\n";
            exit(1);
        }

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

    double CondVHyperLogPrior() const { return -condvhypermean*condvhypermean - condvhypervar; }
    double GeneWHyperLogPrior() const { return -genewhypermean*genewhypermean - genewhypervar; }
    double CondVLogPrior() const { return condvarray->GetLogProb(); }
    double GeneWLogPrior() const { return genewarray->GetLogProb(); }

    double DevPosHyperLogPrior() const {
        double total = 0;
        total += meanpos->GetLogProb();
        total += invshapepos->GetLogProb();
        if (modalprior) {
            for (int k=0; k<Ncond; k++) {
                if (invshapepos->GetVal(k) > 1.0)   {
                    return log(0);
                }
            }
        }
        return total;
    }

    double DevNegHyperLogPrior() const {
        double total = 0;
        total += meanneg->GetLogProb();
        total += invshapeneg->GetLogProb();
        if (modalprior) {
            for (int k=0; k<Ncond; k++) {
                if (invshapeneg->GetVal(k) > 1.0)   {
                    return log(0);
                }
            }
        }
        return total;
    }

    double DevPosLogPrior() const { return devpos->GetLogProb(); }
    double DevNegLogPrior() const { return devneg->GetLogProb(); }

    // double PiLogPrior() const { return pi->GetLogProb(); }
    double AllocLogPrior() const { return alloc->GetLogProb(); }

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
        return condvhypersuffstat.GetLogProb(condvhypermean,condvhypervar);
    }

    double GeneWHyperSuffStatLogProb() const {
        return genewhypersuffstat.GetLogProb(genewhypermean,genewhypervar);
    }

    double DevPosHyperSuffStatLogProb() const   {
        return devposhypersuffstat->GetLogProb(*meanpos,*invshapepos);
    }

    double DevNegHyperSuffStatLogProb() const   {
        return devneghypersuffstat->GetLogProb(*meanneg,*invshapeneg);
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

        double logom = GetMeanLogOmega(gene,cond);

        double om0 = exp(logom);
        double logp0 = count*log(om0) - b*om0;

        if ((!pipos) && (!pineg))  {
            return logp0;
        }

        double ompos = exp(logom + devpos->GetVal(gene).GetVal(cond));
        double logppos = count*log(ompos) - b*ompos;

        double omneg = exp(logom - devneg->GetVal(gene).GetVal(cond));
        double logpneg = count*log(omneg) - b*omneg;

        double max = logp0;
        if (max < logppos)  {
            max = logppos;
        }
        if (max < logpneg)  {
            max = logpneg;
        }
        const vector<double>& w = (*pi)[cond];
        double pneg = w[0]*exp(logpneg-max);
        double p0 = w[1]*exp(logp0-max);
        double ppos = w[2]*exp(logppos-max);
        double tot = pneg + p0 + ppos;
        if (burnin >= 30)   {
            double u = Random::Uniform() * tot;
            if (u < pneg)   {
                (*alloc)[gene][cond] = 0;
            }
            else if (u < (pneg + p0))   {
                (*alloc)[gene][cond] = 1;
            }
            else    {
                (*alloc)[gene][cond] = 2;
            }
        }
        return log(tot) + max;
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

    // log prob for moving meanpos and invshapepos
    double DevPosHyperLogProb() const { return DevPosHyperLogPrior() + DevPosHyperSuffStatLogProb(); }

    // log prob for moving meaneg and invshapeneg
    double DevNegHyperLogProb() const { return DevNegHyperLogPrior() + DevNegHyperSuffStatLogProb(); }

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
            MasterSendOmega();

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
        MasterReceiveLogProbs();
        burnin++;
    }

    // slave move
    void SlaveMove() override {
        GeneResampleSub(1.0);

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            GeneCollectPathSuffStat();
            SlaveSendOmegaSuffStat();
            SlaveReceiveOmega();

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
        SlaveSendLogProbs();
        burnin++;
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
        ScalingMove(lambda, 1.0, 10, &MultiGeneSparseConditionOmegaModel::LambdaHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneSparseConditionOmegaModel::LambdaHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {

        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);

        ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneSparseConditionOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneSparseConditionOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);

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
        ProfileMove(nucrelratehypercenter, 1.0, 1, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10, &MultiGeneSparseConditionOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveNucRates() {
        vector<double> &nucrelrate = (*nucrelratearray)[0];
        ProfileMove(nucrelrate, 0.1, 1, 10, &MultiGeneSparseConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneSparseConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &MultiGeneSparseConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneSparseConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &MultiGeneSparseConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneSparseConditionOmegaModel::TouchNucMatrix, this);

        vector<double> &nucstat = (*nucstatarray)[0];
        ProfileMove(nucstat, 0.1, 1, 10, &MultiGeneSparseConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneSparseConditionOmegaModel::TouchNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 10, &MultiGeneSparseConditionOmegaModel::NucRatesLogProb,
                    &MultiGeneSparseConditionOmegaModel::TouchNucMatrix, this);
    }

    void MultipleTryMove(int nrep, int ntry)  {

        for (int rep=0; rep<nrep; rep++)    {
            for (int cond=0; cond<Ncond; cond++)    {

                vector<int> count(3,0);
                for (int k=0; k<3; k++) {
                    count[k] = 0;
                }
                for (int gene=0; gene<GetNgene(); gene++)   {
                    count[alloc->GetVal(gene).GetVal(cond)]++;
                }

                for (int gene=0; gene<GetNgene(); gene++)    {

                    int bkalloc = alloc->GetVal(gene).GetVal(cond);
                    count[bkalloc]--;
                    if (count[bkalloc] < 0) {
                        cerr << "error in mtry move: count overflow\n";
                        exit(1);
                    }

                    const OmegaPathSuffStat &suffstat = omegapathsuffstatbidimarray->GetVal(gene).GetVal(cond);
                    int c = suffstat.GetCount();
                    double b = suffstat.GetBeta();

                    double logom0 = GetMeanLogOmega(gene,cond);
                    double logp[3];
                    logp[1] = c*logom0 - b*exp(logom0);

                    vector<double> trydevneg(ntry,0);
                    vector<double> trydevpos(ntry,0);

                    vector<double> logpneg(ntry,0);
                    vector<double> logppos(ntry,0);

                    vector<double> cumulpneg(ntry,0);
                    vector<double> cumulppos(ntry,0);

                    int minneg = 0;
                    int minpos = 0;

                    if (bkalloc == 0)   {
                        minneg = 1;
                        trydevneg[0] = devneg->GetVal(gene).GetVal(cond);
                    }
                    else if (bkalloc == 2)  {
                        minpos = 1;
                        trydevpos[0] = devpos->GetVal(gene).GetVal(cond);
                    }

                    double shapeneg = 1.0 / invshapeneg->GetVal(cond);
                    double scaleneg = shapeneg / meanneg->GetVal(cond);
                    for (int k=minneg; k<ntry; k++) {
                        trydevneg[k] = Random::GammaSample(shapeneg,scaleneg);
                    }

                    double maxneg = 0;
                    for (int k=0; k<ntry; k++)  {
                        double logom = logom0 - trydevneg[k];
                        double om = exp(logom);
                        logpneg[k] = c*logom - b*om;
                        if ((!k) || (maxneg < logpneg[k]))  {
                            maxneg = logpneg[k];
                        }
                    }

                    double totneg = 0;
                    for (int k=0; k<ntry; k++)  {
                        double p = exp(logpneg[k] - maxneg);
                        totneg += p;
                        cumulpneg[k] = totneg;
                    }
                    logp[0] = log(totneg / ntry) + maxneg;

                    double uneg = totneg * Random::Uniform();
                    int kneg = 0;
                    while ((kneg < ntry) && (uneg > cumulpneg[kneg]))   {
                        kneg++;
                    }
                    if (kneg == ntry)   {
                        cerr << "error in multiple try: neg overflow\n";
                        exit(1);
                    }

                    double shapepos = 1.0 / invshapepos->GetVal(cond);
                    double scalepos = shapepos / meanpos->GetVal(cond);
                    for (int k=minpos; k<ntry; k++) {
                        trydevpos[k] = Random::GammaSample(shapepos,scalepos);
                    }
                    double maxpos = 0;
                    for (int k=0; k<ntry; k++)  {
                        double logom = logom0 + trydevpos[k];
                        double om = exp(logom);
                        logppos[k] = c*logom - b*om;
                        if ((!k) || (maxpos < logppos[k]))  {
                            maxpos = logppos[k];
                        }
                    }

                    double totpos = 0;
                    for (int k=0; k<ntry; k++)  {
                        double p = exp(logppos[k] - maxpos);
                        totpos += p;
                        cumulppos[k] = totpos;
                    }
                    logp[2] = log(totpos / ntry) + maxpos;

                    double upos = totpos * Random::Uniform();
                    int kpos = 0;
                    while ((kpos < ntry) && (upos > cumulpneg[kpos]))   {
                        kpos++;
                    }
                    if (kpos == ntry)   {
                        cerr << "error in multiple try: pos overflow\n";
                        exit(1);
                    }

                    double max = logp[0];
                    for (int k=1; k<3; k++) {
                        if (max < logp[k])  {
                            max = logp[k];
                        }
                    }

                    double cumul[3];
                    double tot = 0;
                    for (int k=0; k<3; k++) {
                        tot += (piconcentration*picenter[k] + count[k]) * exp(logp[k] - max);
                        cumul[k] = max;
                    }
                    double u = tot*Random::Uniform();
                    int k = 0;
                    while ((k<3) && (cumul[k] < u)) {
                        k++;
                    }
                    if (k == 3) {
                        cerr << "error in multiple try: overflow\n";
                    }

                    (*alloc)[gene][cond] = k;
                    count[k]++;

                    if (k == 0) {
                        (*devneg)[gene][cond] = trydevneg[kneg];
                    }
                    else if (k == 1)    {
                    }
                    else if (k == 2)    {
                        (*devpos)[gene][cond] = trydevpos[kneg];
                    }

                }
            }
        }

        ResamplePi();
    }

    void GammaResampleGeneW()   {

        for (int gene=0; gene<GetNgene(); gene++)   {
            OmegaPathSuffStat suffstat;
            suffstat.Clear();
            for (int cond=0; cond<GetNcond(); cond++)   {
                suffstat += omegapathsuffstatbidimarray->GetVal(gene).GetVal(cond);
            }
            double tmp = Random::GammaSample(1.0 + suffstat.GetCount(), 1.0 + suffstat.GetBeta());
            (*genewarray)[gene] = log(tmp);
        }
        meanlogomegabidimarray->Update();
        condomegabidimarray->Update();
    }

    void MoveOmegaHyperParameters(int nrep) {

        if (burnin <= 10)   {
            GammaResampleGeneW();
        }
        else    {
        for (int rep = 0; rep < nrep; rep++) {
            MoveGeneW(1.0, 1);
            MoveCondV(1.0, 1);
            if (burnin >= 30) {
                if (pipos)    {
                    MoveDevPos(1.0,1);
                }
                if (pineg)  {
                    MoveDevNeg(1.0,1);
                }
            }

            MoveGeneW(0.3, 1);
            MoveCondV(0.3, 1);
            if (burnin >= 30)   {
                if (pipos)    {
                    MoveDevPos(0.3,1);
                }
                if (pineg)  {
                    MoveDevNeg(0.3,1);
                }
            }
        }

        OmegaSuffStatLogProb();
        meanlogomegabidimarray->Update();
        condomegabidimarray->Update();

        if (burnin >= 30)   {
            if (pipos || pineg) {
                ResamplePi();
            }
            if (pipos)  {
                MoveDevPosHyperParams(1.0,100);
            }
            if (pineg)  {
                MoveDevNegHyperParams(1.0,100);
            }
        }
        }
        MoveGeneWHyperParams(1.0, 100);
    }

    void ResamplePi()   {
        for (int cond=0; cond<Ncond; cond++)    {
            int count[3];
            count[0] = count[1] = count[2] = 0;
            for (int gene=0; gene<GetNgene(); gene++)    {
                count[(*alloc)[gene][cond]]++;
            }
            double p[3];
            if (pipos)  {
                p[2] = Random::sGamma(piconcentration*picenter[2] + count[2]);
            }
            else    {
                p[2] = 0;
            }
            if (pineg)  {
                p[0] = Random::sGamma(piconcentration*picenter[0] + count[0]);
            }
            else    {
                p[0] = 0;
            }
            p[1] = Random::sGamma(piconcentration*picenter[1] + count[1]);

            double tot = 0;
            for (int k=0; k<3; k++) {
                tot += p[k];
            }

            for (int k=0; k<3; k++) {
                p[k] /= tot;
                (*pi)[cond][k] = p[k];
            }
        }
    }

    double MoveDevPos(double tuning, int nrep)  {
        double nacc = 0;
        double ntot = 0;
        for (int gene=0; gene<GetNgene(); gene++)    {
            for (int cond=0; cond<Ncond; cond++)    {
                double deltalogprob = -devpos->GetLogProb(gene,cond) - GeneConditionOmegaSuffStatLogProb(gene,cond);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*devpos)[gene][cond] *= e;
                deltalogprob += devpos->GetLogProb(gene,cond) + GeneConditionOmegaSuffStatLogProb(gene,cond);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*devpos)[gene][cond] /= e;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveDevNeg(double tuning, int nrep)  {
        double nacc = 0;
        double ntot = 0;
        for (int gene=0; gene<GetNgene(); gene++)    {
            for (int cond=0; cond<Ncond; cond++)    {
                double deltalogprob = -devneg->GetLogProb(gene,cond) - GeneConditionOmegaSuffStatLogProb(gene,cond);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*devneg)[gene][cond] *= e;
                deltalogprob += devneg->GetLogProb(gene,cond) + GeneConditionOmegaSuffStatLogProb(gene,cond);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*devneg)[gene][cond] /= e;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveCondV(double tuning, int nrep) {
        double nacc = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int j = 1; j < GetNcond(); j++) {
                double deltalogprob = -condvarray->GetLogProb(j);
                deltalogprob -= ConditionOmegaSuffStatLogProb(j);
                double m = tuning * (Random::Uniform() - 0.5);
                (*condvarray)[j] += m;
                deltalogprob += condvarray->GetLogProb(j);
                deltalogprob += ConditionOmegaSuffStatLogProb(j);
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*condvarray)[j] -= m;
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
                (*genewarray)[i] += m;
                deltalogprior += genewarray->GetLogProb(i);
                deltasuffstatlogprob += GeneOmegaSuffStatLogProb(i);
                double deltalogprob = deltalogprior + deltasuffstatlogprob;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*genewarray)[i] -= m;
                }
            }
        }
        return nacc / GetNgene() / nrep;
    }

    void MoveDevPosHyperParams(double tuning, int nrep) {
        devposhypersuffstat->Clear();
        devposhypersuffstat->AddSuffStat(*devpos);

        for (int cond=0; cond<Ncond; cond++)    {
            ScalingMove((*meanpos)[cond], 1.0, 10, &MultiGeneSparseConditionOmegaModel::DevPosHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
            ScalingMove((*meanpos)[cond], 0.3, 10, &MultiGeneSparseConditionOmegaModel::DevPosHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
            ScalingMove((*invshapepos)[cond], 1.0, 10, &MultiGeneSparseConditionOmegaModel::DevPosHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
            ScalingMove((*invshapepos)[cond], 0.3, 10, &MultiGeneSparseConditionOmegaModel::DevPosHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        }
    }

    void MoveDevNegHyperParams(double tuning, int nrep) {
        devneghypersuffstat->Clear();
        devneghypersuffstat->AddSuffStat(*devneg);

        for (int cond=0; cond<Ncond; cond++)    {
            ScalingMove((*meanneg)[cond], 1.0, 10, &MultiGeneSparseConditionOmegaModel::DevNegHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
            ScalingMove((*meanneg)[cond], 0.3, 10, &MultiGeneSparseConditionOmegaModel::DevNegHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
            ScalingMove((*invshapeneg)[cond], 1.0, 10, &MultiGeneSparseConditionOmegaModel::DevNegHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
            ScalingMove((*invshapeneg)[cond], 0.3, 10, &MultiGeneSparseConditionOmegaModel::DevNegHyperLogProb,
                        &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        }
    }

    void MoveCondVHyperParams(double tuning, int nrep) {
        condvhypersuffstat.Clear();
        condvhypersuffstat.AddSuffStat(*condvarray);

        ScalingMove(condvhypervar, 1.0, 10, &MultiGeneSparseConditionOmegaModel::CondVHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(condvhypervar, 0.3, 10, &MultiGeneSparseConditionOmegaModel::CondVHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);

        condvarray->SetMean(condvhypermean);
        condvarray->SetVar(condvhypervar);
    }

    void MoveGeneWHyperParams(double tuning, int nrep) {
        genewhypersuffstat.Clear();
        genewhypersuffstat.AddSuffStat(*genewarray);

        SlidingMove(genewhypermean, 1.0, 10, 0, 0, &MultiGeneSparseConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        SlidingMove(genewhypermean, 0.3, 10, 0, 0, &MultiGeneSparseConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(genewhypervar, 1.0, 10, &MultiGeneSparseConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);
        ScalingMove(genewhypervar, 0.3, 10, &MultiGeneSparseConditionOmegaModel::GeneWHyperLogProb,
                    &MultiGeneSparseConditionOmegaModel::NoUpdate, this);

        genewarray->SetMean(genewhypermean);
        genewarray->SetVar(genewhypervar);
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

    void MasterSendOmega() { 
        condomegabidimarray->Update();
        MasterSendGeneArray(*condomegabidimarray); 
    }

    void SlaveReceiveOmega() {
        SlaveReceiveGeneArray(*condomegabidimarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmegaTree(condomegabidimarray->GetVal(gene));
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
