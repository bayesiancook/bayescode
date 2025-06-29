
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "NHT92SingleOmegaModel.hpp"
#include "IIDDirichlet.hpp"
#include "IIDGamma.hpp"

class MultiGeneSingleOmegaModel : public MultiGeneProbModel {

  private:

    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;
    std::vector<CodonSequenceAlignment*> alivector;

    string datafile;
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
    double rootkappa;
    double kappahypermean;
    double kappahyperinvshape;
    BranchIIDGamma* branchkappa;

    GammaSuffStat kappahypersuffstat;

    double rootgamma;
    double gammahypermean;
    double gammahyperinvconc;
    BranchIIDBeta* branchgamma;

    BetaSuffStat gammahypersuffstat;

    T92SubMatrix* rootnucmatrix;
    T92SubMatrixBranchArray* nucmatrixarray;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    // includes root suffstat in addition to branch suff stats
    NucPathSuffStatBranchArray* nucpathsuffstat;

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

    // each gene defines its own SingleOmegaModel
    std::vector<SingleOmegaModel *> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneSingleOmegaModel(string indatafile, string intreefile, int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs)	{

        blmode = 1;
        nucmode = 1;
        omegamode = 1;

        datafile = indatafile;
        treefile = intreefile;
        AllocateAlignments(datafile);

        refcodondata = new CodonSequenceAlignment(refdata, true);
        taxonset = refdata->GetTaxonSet();
        Ntaxa = refdata->GetNtaxa();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        if (!myid) {
            cerr << "number of taxa : " << Ntaxa << '\n';
            cerr << "number of branches : " << Nbranch << '\n';
            cerr << "tree and data fit together\n";
        }
    }

    const Tree& GetTree() const {return *tree;}

    int GetNbranch() const  {
        return Nbranch;
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
            branchlengtharray =
                new GammaWhiteNoiseArray(GetLocalNgene(), *tree, *branchlength, 1.0 / blhyperinvshape);
            lengthpathsuffstatarray = 0;
            lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
        }

        // Nucleotide rates

        if (nucmode == 2) {
            rootkappa = 2.0;
            kappahypermean = 2.0;
            kappahyperinvshape = 0.1;
            double kappaalpha = 1.0 / kappahyperinvshape;
            double kappabeta = kappaalpha / kappahypermean;
            branchkappa = new BranchIIDGamma(*tree, kappaalpha, kappabeta);
            branchkappa->SetAllBranches(2.0);

            rootgamma = 0.5;
            gammahypermean = 0.5;
            gammahyperinvconc = 0.1;
            double gammaalpha = gammahypermean / gammahyperinvconc;
            double gammabeta = (1-gammahypermean) / gammahyperinvconc;
            branchgamma = new BranchIIDBeta(*tree, gammaalpha, gammabeta);
            branchgamma->SetAllBranches(0.5);

            rootnucmatrix = new T92SubMatrix(rootkappa, rootgamma, true);
            nucmatrixarray = new T92SubMatrixBranchArray(branchkappa, branchgamma, true);
            nucpathsuffstat = new NucPathSuffStatBranchArray(*tree);
        }
        else    {
            cerr << "gene specific nuc rates not yet implemented\n";
            exit(1);
        }

        // Omega

        omegaarray = new IIDGamma(GetLocalNgene(), omegahypermean, omegahyperinvshape);

        // Gene processes 

        lnL = 0;
        GeneLogPrior = 0;

        if (!GetMyid()) {
            geneprocess.assign(0, (SingleOmegaModel*)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (SingleOmegaModel*)0);

            ifstream is(datafile.c_str());
            string tmp;
            is >> tmp;
            if (tmp == "ALI")   {
                int ngene;
                is >> ngene;
                if (ngene != GetNgene())    {
                    cerr << "error when reading alignments from cat file: non matching number of genes\n";
                    exit(1);
                }
                alivector.assign(GetLocalNgene(), (CodonSequenceAlignment*) 0);
                int index = 0;
                for (int gene=0; gene<GetNgene(); gene++)   {
                    string name;
                    is >> name;
                    FileSequenceAlignment tmp(is);
                    if ((index < GetLocalNgene()) && (name == GeneName[index]))    {
                        if (GetLocalGeneName(index) != name)    {
                            cerr << "error: non matching gene name\n";
                            exit(1);
                        }
                        if (alivector[index]) {
                            cerr << "error: alignment already allocated\n";
                            exit(1);
                        }
                        alivector[index] = new CodonSequenceAlignment(&tmp, true);
                        index++;
                    }
                }
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    if (! alivector[gene])  {
                        cerr << "error: alignment not allocated\n";
                        exit(1);
                    }
		    geneprocess[gene] = new SingleOmegaModel(alivector[gene], tree);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
		    geneprocess[gene] = new SingleOmegaModel(GetLocalGeneName(gene), treefile);
                }
            }

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetAcrossGenesModes(blmode, nucmode);
                geneprocess[gene]->Allocate();
            }
        }
    }

    // called upon constructing the model
    // mode == 2: global 
    // mode == 1: gene specific, with hyperparameters estimated across genes
    // mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode, int inomegamode)   {
        blmode = inblmode;
        nucmode = innucmode;
        omegamode = inomegamode;
    }

    int GetBLMode() const   {
        return blmode;
    }

    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape)  {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    void FastUpdate() {

        branchlength->SetScale(lambda);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }
        if (nucmode == 2)   {
            double kappaalpha = 1.0 / kappahyperinvshape;
            double kappabeta = kappaalpha / kappahypermean;
            branchkappa->SetShape(kappaalpha);
            branchkappa->SetScale(kappabeta);
            branchgamma->SetMeanInvConc(gammahypermean, gammahyperinvconc);
            TouchNucMatrices();
        }

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
    }

    void MasterUpdate() override {

        FastUpdate();

        if (nprocs > 1) {
            MasterSendBranchLengthsHyperParameters();
            // MasterSendNucRatesHyperParameters();

            if (blmode == 2) {
                MasterSendGlobalBranchLengths();
            } else {
                MasterSendGeneBranchLengths();
            }

            if (nucmode == 2) {
                MasterSendGlobalNucRates();
            } else {
                // MasterSendGeneNucRates();
            }

            MasterSendOmegaHyperParameters();
            MasterSendOmega();
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {

        SlaveReceiveBranchLengthsHyperParameters();
        // SlaveReceiveNucRatesHyperParameters();

        if (blmode == 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveGeneBranchLengths();
        }
        if (nucmode == 2) {
            SlaveReceiveGlobalNucRates();
        } else {
            // SlaveReceiveGeneNucRates();
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
        if (nprocs > 1) {
            MasterSendBranchLengthsHyperParameters();
            // MasterSendNucRatesHyperParameters();

            if (blmode == 2) {
                MasterSendGlobalBranchLengths();
            } else {
                MasterSendGeneBranchLengths();
            }

            if (nucmode == 2) {
                MasterSendGlobalNucRates();
            } else {
                // MasterSendGeneNucRates();
            }

            MasterSendOmegaHyperParameters();
            MasterSendOmega();
        }
    }

    void SlavePostPred(string name) override {
        SlaveReceiveBranchLengthsHyperParameters();
        // SlaveReceiveNucRatesHyperParameters();

        if (blmode == 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveGeneBranchLengths();
        }
        if (nucmode == 2) {
            SlaveReceiveGlobalNucRates();
        } else {
            // SlaveReceiveGeneNucRates();
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

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    const vector<double> &GetOmegaArray() const { return omegaarray->GetArray(); }

    //-------------------
    // Traces and Monitors
    //-------------------

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL";
        if (blmode == 2) {
            os << "\tlength";
        } else {
            os << "\tmeanlength\tstdev";
        }
        os << "\tmeanomega";
        os << "\tvaromega";
        os << "\tomegahypermean\tinvshape";
        os << "\tmeantstv\tinvshape\trootgc\tmeangc\tinvconc";
        os << '\n';
    }

    void Trace(ostream &os) const override {
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

        os << '\t' << kappahypermean;
        os << '\t' << kappahyperinvshape;
        os << '\t' << rootgamma;
        os << '\t' << gammahypermean << '\t' << gammahyperinvconc;
        os << '\n';
        os.flush();
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

    void AddLength(SimpleBranchArray<double>& in)   {
        for (int j = 0; j < Nbranch; j++) {
            in[j] += branchlength->GetVal(j);
        }
    }

    // Nucleotide rates

    void Monitor(ostream &os) const override {}

    void MasterFromStream(istream &is) override {

        if (blmode == 2) {
            is >> lambda;
            is >> *branchlength;
        } else {
            is >> lambda;
            is >> *branchlength;
            is >> blhyperinvshape;
            is >> *branchlengtharray;
        }

        if (nucmode == 2)   {
            is >> kappahypermean >> kappahyperinvshape >> *branchkappa;
            is >> rootgamma >> gammahypermean >> gammahyperinvconc >> *branchgamma;
        }

        is >> omegahypermean;
        is >> omegahyperinvshape;
        is >> *omegaarray;
    }

    void MasterToStream(ostream &os) const override {

        if (blmode == 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        } else {
            os << lambda << '\t';
            os << *branchlength << '\t';
            os << blhyperinvshape << '\t';
            os << *branchlengtharray << '\t';
        }

        if (nucmode == 2)   {
            os << kappahypermean << '\t';
            os << kappahyperinvshape << '\t';
            os << *branchkappa << '\t';
            os << rootgamma << '\t';
            os << gammahypermean << '\t';
            os << gammahyperinvconc << '\t';
            os << *branchgamma << '\t';
        }

        os << omegahypermean << '\t';
        os << omegahyperinvshape << '\t';
        os << *omegaarray << '\n';
    }

    void TraceOmega(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << omegaarray->GetVal(gene) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TraceGeneTreeLength(ostream &os) const {
        if (blmode == 2)    {
            cerr << "error: in trace gene length tree yet branch lengths are shared\n";
            exit(1);
        }
        for (int gene = 0; gene < Ngene; gene++) {
            os << branchlengtharray->GetVal(gene).GetTotalLength() << '\t';
        }
        os << '\n';
        os.flush();
    }

    //-------------------
    // Updates
    //-------------------

    void NoUpdate() {}

    void TouchNucMatrices() {
        TouchRootNucMatrix();
        nucmatrixarray->UpdateMatrices();
    }

    void TouchRootNucMatrix()	{
        rootnucmatrix->SetKappa(rootkappa);
        rootnucmatrix->SetGC(rootgamma);
        rootnucmatrix->CorruptMatrix();
    }

    void TouchBranchNucMatrix(int branchindex)  {
        (*nucmatrixarray)[branchindex].SetKappa(branchkappa->GetVal(branchindex));
        (*nucmatrixarray)[branchindex].SetGC(branchgamma->GetVal(branchindex));
        (*nucmatrixarray)[branchindex].CorruptMatrix();
    }

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
            total += NucRatesHyperLogPrior();
            total += NucRatesLogPrior();
        }

        if (omegamode == 1) {
            total += OmegaHyperLogPrior();
        }
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

    double NucRatesHyperLogPrior() const    {
        double total = 0;
        total -= gammahyperinvconc;
        total -= kappahypermean / 5.0;
        total -= kappahyperinvshape;
        return total;
    }

    double KappaHyperLogPrior() const   {
        return -kappahypermean / 5.0 - kappahyperinvshape;
    }

    double GammaHyperLogPrior() const    {
        return -gammahyperinvconc;
    }

    double NucRatesLogPrior() const {
        double total = 0;

        // uniform prior over rootgamma
        // root kappa is fixed

        total += branchkappa->GetLogProb();
        total += branchgamma->GetLogProb();
        return total;
    }

    double RootNucRatesLogPrior() const {
        return 0;
    }

    double NucRatesLogPrior(int branchindex) const  {
        return branchkappa->GetLogProb(branchindex) + branchgamma->GetLogProb(branchindex);
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

    double GammaHyperSuffStatLogProb() const {
        return gammahypersuffstat.GetMeanInvConcLogProb(gammahypermean, gammahyperinvconc);
    }

    double KappaHyperSuffStatLogProb() const    {
        double kappaalpha = 1.0 / kappahyperinvshape;
        double kappabeta = kappaalpha / kappahypermean;
        return kappahypersuffstat.GetLogProb(kappaalpha, kappabeta);
    }

    void CollectNucPathSuffStat() {
        // collect MPI
    }

    double NucRatesSuffStatLogProb(int branchindex) const   {
        return nucpathsuffstat->GetVal(branchindex).GetLogProb(nucmatrixarray->GetVal(branchindex), *GetCodonStateSpace());
    }

    double RootNucRatesSuffStatLogProb() const  {
        return nucpathsuffstat->GetRootVal().GetLogProb(*rootnucmatrix, *GetCodonStateSpace());
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

    double GammaHyperLogProb() const    {
        return GammaHyperLogPrior() + GammaHyperSuffStatLogProb();
    }

    double KappaHyperLogProb() const    {
        return KappaHyperLogPrior() + KappaHyperSuffStatLogProb();
    }


    // log prob for moving nuc rates hyper params
    /*
    double NucRatesHyperLogProb() const {
        return GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
    }
    */

    // log prob for moving nuc rates
    double NucRatesLogProb(int branchindex) const   {
        return NucRatesLogPrior(branchindex) + NucRatesSuffStatLogProb(branchindex);
    }

    double RootNucRatesLogProb() const  {
        return RootNucRatesLogPrior() + RootNucRatesSuffStatLogProb();
    }

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

    void MasterMove() override {
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            if (omegamode == 1) {
                MasterReceiveOmega();
                MoveOmegaHyperParameters();
                MasterSendOmegaHyperParameters();
            }

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
                MoveNucRatesHyperParameters();
                MasterSendGlobalNucRates();
            } else if (nucmode == 1) {
                /*
                MasterReceiveNucRatesHyperSuffStat();
                MoveNucRatesHyperParameters();
                MasterSendNucRatesHyperParameters();
                */
            }
        }

        // collect current state
        if (blmode != 2) {
            MasterReceiveGeneBranchLengths();
        }
        if (nucmode != 2) {
            // MasterReceiveGeneNucRates();
        }
        MasterReceiveOmega();
        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {
        GeneResampleSub(1.0);

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {

            MoveGeneParameters(1.0);

            if (omegamode == 1) {
                SlaveSendOmega();
                SlaveReceiveOmegaHyperParameters();
            }

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                SlaveSendBranchLengthsSuffStat();
                SlaveReceiveGlobalBranchLengths();
            } else if (blmode == 1) {
                SlaveSendBranchLengthsHyperSuffStat();
                SlaveReceiveBranchLengthsHyperParameters();
            }

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                SlaveSendNucPathSuffStat();
                SlaveReceiveGlobalNucRates();
            } else if (nucmode == 1) {
                /*
                SlaveSendNucRatesHyperSuffStat();
                SlaveReceiveNucRatesHyperParameters();
                */
            }
        }

        // collect current state
        if (blmode != 2) {
            SlaveSendGeneBranchLengths();
        }
        if (nucmode != 2) {
            // SlaveSendGeneNucRates();
        }
        SlaveSendOmega();
        SlaveSendLogProbs();
    }

    void GeneResampleSub(double frac) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void MoveGeneParameters(int nrep)   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveParameters(nrep);

            (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
            if (blmode != 2) {
                geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
            }
            if (nucmode != 2) {
                // geneprocess[gene]->GetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
            }
        }
    }

    // Branch lengths

    void ResampleBranchLengths() {
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneSingleOmegaModel::LambdaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneSingleOmegaModel::LambdaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {

        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);

        ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneSingleOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneSingleOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);

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

    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        CollectNucPathSuffStat();

        SlidingMove(rootgamma, 0.3, 3, 0, 1, &MultiGeneSingleOmegaModel::RootNucRatesLogProb,
                & MultiGeneSingleOmegaModel::TouchRootNucMatrix, this);
        SlidingMove(rootgamma, 0.1, 3, 0, 1, &MultiGeneSingleOmegaModel::RootNucRatesLogProb,
                & MultiGeneSingleOmegaModel::TouchRootNucMatrix, this);

        for (int j=0; j<GetNbranch(); j++)  {
            MoveBranchGamma(j, 0.3, 3);
            MoveBranchKappa(j, 0.3, 3);
            MoveBranchGamma(j, 0.1, 3);
            MoveBranchKappa(j, 0.1, 3);
        }

        // important for codon matrices (nucl. matrices are already updated)
        TouchNucMatrices();
    }

    double MoveBranchGamma(int branch, double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        double& x = (*branchgamma)[branch];
        double min = 0;
        double max = 1.0;

        for (int rep = 0; rep < nrep; rep++) {
            double bk = x;
            double deltalogprob = - NucRatesLogProb(branch);
            double m = tuning * (Random::Uniform() - 0.5);
            x += m;
            if (max > min) {
                while ((x < min) || (x > max)) {
                    if (x < min) {
                        x = 2 * min - x;
                    }
                    if (x > max) {
                        x = 2 * max - x;
                    }
                }
            }
            TouchBranchNucMatrix(branch);
            deltalogprob += NucRatesLogProb(branch);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x = bk;
                TouchBranchNucMatrix(branch);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveBranchKappa(int branch, double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        double& x = (*branchkappa)[branch];

        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -NucRatesLogProb(branch);
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            x *= e;
            TouchBranchNucMatrix(branch);
            deltalogprob += NucRatesLogProb(branch);
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x /= e;
                TouchBranchNucMatrix(branch);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    void MoveNucRatesHyperParameters()  {
        MoveKappaHyperParameters();
        MoveGammaHyperParameters();
    }

    void MoveKappaHyperParameters() {
        kappahypersuffstat.Clear();
        kappahypersuffstat.AddSuffStat(*branchkappa);
        ScalingMove(kappahypermean, 1.0, 10, &MultiGeneSingleOmegaModel::KappaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        ScalingMove(kappahyperinvshape, 1.0, 10, &MultiGeneSingleOmegaModel::KappaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        double alpha = 1.0 / kappahyperinvshape;
        double beta = alpha / kappahypermean;
        branchkappa->SetShape(alpha);
        branchkappa->SetScale(beta);
    }

    void MoveGammaHyperParameters() {
        gammahypersuffstat.Clear();
        branchgamma->AddSuffStat(gammahypersuffstat);
        ScalingMove(gammahypermean, 1.0, 10, &MultiGeneSingleOmegaModel::GammaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        ScalingMove(gammahyperinvconc, 1.0, 10, &MultiGeneSingleOmegaModel::GammaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        branchgamma->SetMeanInvConc(gammahypermean, gammahyperinvconc);
    }

    // Omega

    void MoveOmegaHyperParameters() {
        omegahypersuffstat.Clear();
        omegahypersuffstat.AddSuffStat(*omegaarray);

        ScalingMove(omegahypermean, 1.0, 10, &MultiGeneSingleOmegaModel::OmegaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        ScalingMove(omegahypermean, 0.3, 10, &MultiGeneSingleOmegaModel::OmegaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        ScalingMove(omegahyperinvshape, 1.0, 10, &MultiGeneSingleOmegaModel::OmegaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);
        ScalingMove(omegahyperinvshape, 0.3, 10, &MultiGeneSingleOmegaModel::OmegaHyperLogProb,
                    &MultiGeneSingleOmegaModel::NoUpdate, this);

        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
        omegaarray->SetShape(alpha);
        omegaarray->SetScale(beta);
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
        MasterSendGlobal(rootkappa, rootgamma);
        MasterSendGlobal(*branchkappa, *branchgamma);
    }

    void SlaveReceiveGlobalNucRates() {
        SlaveReceiveGlobal(rootkappa, rootgamma);
        SlaveReceiveGlobal(*branchkappa, *branchgamma);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates(rootkappa, rootgamma, *branchkappa, *branchgamma);
        }
    }

    /*
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
    */

    void SlaveSendNucPathSuffStat() {
        nucpathsuffstat->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat->Add(geneprocess[gene]->GetNucPathSuffStat());
        }
        SlaveSendAdditive(*nucpathsuffstat);
    }

    void MasterReceiveNucPathSuffStat() {
        nucpathsuffstat->Clear();
        MasterReceiveAdditive(*nucpathsuffstat);
    }

    // omega (and hyperparameters)

    void SlaveSendOmega() { SlaveSendGeneArray(*omegaarray); }

    void MasterReceiveOmega() { MasterReceiveGeneArray(*omegaarray); }

    void MasterSendOmega() { MasterSendGeneArray(*omegaarray); }

    void SlaveReceiveOmega() {
        SlaveReceiveGeneArray(*omegaarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmega((*omegaarray)[gene]);
        }
    }

    // omega hyperparameters

    void MasterSendOmegaHyperParameters() { MasterSendGlobal(omegahypermean, omegahyperinvshape); }

    void SlaveReceiveOmegaHyperParameters() {
        SlaveReceiveGlobal(omegahypermean, omegahyperinvshape);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape);
        }
    }

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

    /*
    void SlaveAddGeneNodePathSuffStat(vector<RelativePathSuffStatNodeArray>& array)   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->AddNodePathSuffStat(array[gene]);
        }
    }

    void SlaveAdddSOmegaPathSuffStat(vector<dSOmegaPathSuffStatBranchArray>& array)   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->AdddSOmegaPathSuffStat(array[gene]);
        }
    }

    void SlaveAddGCConsdSOmegaPathSuffStat(vector<GCConsdSOmegaPathSuffStatBranchArray>& array)   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->AddGCConsdSOmegaPathSuffStat(array[gene]);
        }
    }

    void SlaveAddGeneDoubleCounts(vector<vector<double>>& counts) const {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->AddDoubleCounts(counts[gene]);
        }
    }
    */
};
