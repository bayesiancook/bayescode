
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

#include "Chrono.hpp"
#include "DiffSelDoublySparseModel.hpp"
#include "IIDMultiBernBeta.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"

class MultiGeneDiffSelDoublySparseModel : public MultiGeneProbModel {
  private:
    double epsilon;
    double fitnessshape;
    int fitnesscentermode;
    int blmode;
    int nucmode;
    int shiftmode;

    const double minshiftprobhypermean = 0.01;
    const double maxshiftprobhypermean = 0.3;
    const double maxshiftprobhyperinvconc = 0.2;
    const double minpi = 0.01;
    const double maxpi = 0.50;

    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;

    string treefile;

    int Ntaxa;
    int Nbranch;

    int Ncond;
    int Nlevel;
    int codonmodel;

    // branch lengths
    double lambda;
    BranchIIDGamma *branchlength;
    GammaSuffStat hyperlengthsuffstat;

    double blhyperinvshape;
    GammaWhiteNoiseArray *branchlengtharray;
    GammaSuffStatBranchArray *lengthhypersuffstatarray;

    // gene-specific nuc rates
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    IIDDirichlet *nucrelratearray;
    DirichletSuffStat nucrelratesuffstat;

    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    IIDDirichlet *nucstatarray;
    DirichletSuffStat nucstatsuffstat;

    // shiftprob arrays across genes
    double shiftprobmean;
    double shiftprobinvconc;
    vector<double> shiftprobhypermean;
    vector<double> shiftprobhyperinvconc;
    vector<double> pi;
    double pihypermean;
    double pihyperinvconc;

    IIDMultiBernBeta *shiftprobarray;
    SimpleArray<int> *totcount;
    IIDMultiCount *shiftcountarray;

    // each gene defines its own DiffSelDoublySparseModel
    std::vector<DiffSelDoublySparseModel *> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;
    double meanWidth;
    double meanEps;

    SimpleArray<double> *genednds;

    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    int withtoggle;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneDiffSelDoublySparseModel(string datafile, string intreefile, int inNcond, int inNlevel,
                                      int incodonmodel, double inepsilon, double infitnessshape,
                                      int inblmode, int innucmode, int inshiftmode,
                                      double inpihypermean, double inpihyperinvconc,
                                      double inshiftprobmean, double inshiftprobinvconc, int inmyid,
                                      int innprocs)
        : MultiGeneProbModel(inmyid, innprocs), nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {
        withtoggle = 0;

        blmode = inblmode;
        nucmode = innucmode;
        shiftmode = inshiftmode;
        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;
        shiftprobmean = inshiftprobmean;
        shiftprobinvconc = inshiftprobinvconc;

        epsilon = inepsilon;
        fitnessshape = infitnessshape;
        fitnesscentermode = 3;

        codonmodel = incodonmodel;
        Ncond = inNcond;
        Nlevel = inNlevel;

        AllocateAlignments(datafile);
        treefile = intreefile;

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
            std::cerr << "number of taxa : " << Ntaxa << '\n';
            std::cerr << "number of branches : " << Nbranch << '\n';
            std::cerr << "-- Tree and data fit together\n";
        }

        GeneLogPrior = 0;
        lnL = 0;
        geneprocess.assign(0, (DiffSelDoublySparseModel *)0);
    }

    void Allocate() {
        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
        blhyperinvshape = 0.1;
        branchlength->SetAllBranches(1.0 / lambda);
        branchlengtharray =
            new GammaWhiteNoiseArray(GetLocalNgene(), *tree, *branchlength, 1.0 / blhyperinvshape);
        lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;
        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;
        nucrelratearray =
            new IIDDirichlet(GetLocalNgene(), nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        nucstatarray =
            new IIDDirichlet(GetLocalNgene(), nucstathypercenter, 1.0 / nucstathyperinvconc);

        shiftprobhypermean.assign(Ncond - 1, shiftprobmean);
        shiftprobhyperinvconc.assign(Ncond - 1, shiftprobinvconc);
        pi.assign(Ncond - 1, pihypermean);
        shiftprobarray =
            new IIDMultiBernBeta(GetLocalNgene(), pi, shiftprobhypermean, shiftprobhyperinvconc);
        totcount = new SimpleArray<int>(GetLocalNgene(), 0);
        shiftcountarray =
            new IIDMultiCount(totcount->GetArray(), pi, shiftprobhypermean, shiftprobhyperinvconc);

        lnL = 0;
        GeneLogPrior = 0;

        if (Ncond == 1) {
            genednds = new SimpleArray<double>(GetLocalNgene());
        } else {
            genednds = 0;
        }

        if (!GetMyid()) {
            geneprocess.assign(0, (DiffSelDoublySparseModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (DiffSelDoublySparseModel *)0);

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene] = new DiffSelDoublySparseModel(
                    GetLocalGeneName(gene), treefile, Ncond, Nlevel, codonmodel, epsilon,
                    fitnessshape, pihypermean, shiftprobmean, shiftprobinvconc);
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->SetFitnessCenterMode(fitnesscentermode);
                geneprocess[gene]->SetWithToggles(withtoggle);
                geneprocess[gene]->Allocate();
            }
        }
    }

    void FastUpdate() {
        branchlength->SetScale(lambda);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MasterUpdate() override {
        FastUpdate();
        if (nprocs > 1) {
            MasterSendBranchLengthsHyperParameters();
            MasterSendGeneBranchLengths();
            MasterSendNucRatesHyperParameters();
            MasterSendGeneNucRates();
            MasterSendShiftProbHyperParameters();
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {
        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveGeneBranchLengths();
        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveGeneNucRates();
        SlaveReceiveShiftProbHyperParameters();
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
            MasterSendGeneBranchLengths();
            MasterSendNucRatesHyperParameters();
            MasterSendGeneNucRates();
            MasterSendShiftProbHyperParameters();
        }
    }

    void SlavePostPred(string name) override {
        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveGeneBranchLengths();
        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveGeneNucRates();
        SlaveReceiveShiftProbHyperParameters();
        GenePostPred(name);
    }

    void GenePostPred(string name) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->PostPred(name + GetLocalGeneName(gene));
        }
    }

    void SetWithToggles(int in) {
        withtoggle = in;
        if (myid && geneprocess.size()) {
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetWithToggles(in);
            }
        }
    }

    //! \brief set estimation method for fitness hyperparameter (center of
    //! Dirichlet distribution)
    //!
    //! - mode == 3: fixed (uniform)
    //! - mode == 2: shared across genes, estimated
    //! - mode == 1: gene specific, with hyperparameters estimated across genes
    //! - mode == 0: gene-specific, with fixed hyperparameters
    void SetFitnessCenterMode(int in) { fitnesscentermode = in; }

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    int GetNbranch() const { return tree->GetNbranch(); }

    const Tree *GetTree() const { return tree; }

    //-------------------
    // Traces and Monitors
    //-------------------

    void TracePredictedDNDS(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << genednds->GetVal(gene) << '\t';
        }
        os << '\n';
    }

    void MasterReceivePredictedDNDS() { MasterReceiveGeneArray(*genednds); }

    void SlaveSendPredictedDNDS() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            (*genednds)[gene] = geneprocess[gene]->GetPredictedDNDS(0);
        }
        SlaveSendGeneArray(*genednds);
    }

    double GetMeanLength() const { return branchlengtharray->GetMeanLength(); }

    double GetVarLength() const { return branchlengtharray->GetVarLength(); }

    void TraceHeader(ostream &os) const {
        os << "#logprior\tlnL\tlength\t";
        os << "invshape\t";
        os << "meanwidth\t";
        os << "meaneps\t";
        for (int k = 1; k < Ncond; k++) {
            os << "pi" << k << '\t';
            os << "mean\t";
            os << "invconc\t";
        }
        os << "nucstatcenter\tinvconc\trelratecenter\tinvconc\t";
        os << '\n';
    }

    void Trace(ostream &os) const {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t' << blhyperinvshape << '\t';
        // os << GetMeanLength() << '\t' << GetVarLength() << '\t';
        os << meanWidth << '\t';
        os << meanEps << '\t';
        for (int k = 1; k < Ncond; k++) {
            os << pi[k - 1] << '\t';
            os << shiftprobhypermean[k - 1] << '\t';
            os << shiftprobhyperinvconc[k - 1] << '\t';
        }
        os << Random::GetEntropy(nucstathypercenter) << '\t';
        os << nucstathyperinvconc << '\t';
        os << Random::GetEntropy(nucrelratehypercenter) << '\t';
        os << nucrelratehyperinvconc << '\t';
        os << '\n';
        os.flush();
    }

    void Monitor(ostream &os) const {}

    void MasterToStream(ostream &os) const {
        os << lambda << '\t';
        os << *branchlength << '\t';
        os << blhyperinvshape << '\t';
        os << *branchlengtharray << '\t';

        os << nucrelratehypercenter << '\t';
        os << nucrelratehyperinvconc << '\t';
        os << nucstathypercenter << '\t';
        os << nucstathyperinvconc << '\t';
        os << *nucrelratearray << '\t';
        os << *nucstatarray << '\t';

        os << shiftprobhypermean << '\t';
        os << shiftprobhyperinvconc << '\t';
        os << pi << '\t';

        for (int proc = 1; proc < GetNprocs(); proc++) {
            MPI_Status stat;
            int size;
            MPI_Recv(&size, 1, MPI_INT, proc, TAG1, MPI_COMM_WORLD, &stat);
            MPIBuffer buffer(size);
            MPI_Recv(buffer.GetBuffer(), buffer.GetSize(), MPI_DOUBLE, proc, TAG1, MPI_COMM_WORLD,
                     &stat);
            os << size << '\n';
            buffer.ToStream(os);
        }
    }

    void SlaveToStream() const override {
        int size = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            size += geneprocess[gene]->GetMPISize();
        }
        MPIBuffer buffer(size);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            buffer << *geneprocess[gene];
        }
        MPI_Send(&size, 1, MPI_INT, 0, TAG1, MPI_COMM_WORLD);
        MPI_Send(buffer.GetBuffer(), buffer.GetSize(), MPI_DOUBLE, 0, TAG1, MPI_COMM_WORLD);
    }

    void MasterFromStream(istream &is) {
        is >> lambda;
        is >> *branchlength;
        is >> blhyperinvshape;
        is >> *branchlengtharray;

        is >> nucrelratehypercenter;
        is >> nucrelratehyperinvconc;
        is >> nucstathypercenter;
        is >> nucstathyperinvconc;
        is >> *nucrelratearray;
        is >> *nucstatarray;

        is >> shiftprobhypermean;
        is >> shiftprobhyperinvconc;
        is >> pi;

        for (int proc = 1; proc < GetNprocs(); proc++) {
            int size;
            is >> size;
            MPI_Send(&size, 1, MPI_INT, proc, TAG1, MPI_COMM_WORLD);
            MPIBuffer buffer(size);
            buffer.FromStream(is);
            MPI_Send(buffer.GetBuffer(), buffer.GetSize(), MPI_DOUBLE, proc, TAG1, MPI_COMM_WORLD);
        }
    }

    void SlaveFromStream() override {
        int checksize = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            checksize += geneprocess[gene]->GetMPISize();
        }
        MPI_Status stat;
        int size;
        MPI_Recv(&size, 1, MPI_INT, 0, TAG1, MPI_COMM_WORLD, &stat);
        if (size != checksize) {
            cerr << "error in SlaveFromStream: non matching buffer size\n";
            exit(1);
        }
        MPIBuffer buffer(size);
        MPI_Recv(buffer.GetBuffer(), buffer.GetSize(), MPI_DOUBLE, 0, TAG1, MPI_COMM_WORLD, &stat);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            buffer >> *geneprocess[gene];
        }
    }

    double GetSlaveMoveTime() const { return moveTime; }

    double GetSlaveMapTime() const { return mapTime; }

    double GetMasterMoveTime() const { return movechrono.GetTime(); }

    //-------------------
    // Updates
    //-------------------

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
        double total = GeneLogPrior;
        if (blmode == 1) {
            total += GeneBranchLengthsHyperLogPrior();
        }
        if (nucmode == 1) {
            total += GeneNucRatesHyperLogPrior();
        }
        total += ShiftProbHyperLogPrior();
        if (std::isnan(total)) {
            cerr << "GetLogPrior is nan\n";
            exit(1);
        }
        return total;
    }

    double LambdaHyperLogPrior() const { return -lambda / 10; }

    // exponential of mean 1 for blhyperinvshape
    double BranchLengthsHyperInvShapeLogPrior() const { return -blhyperinvshape; }

    double GeneBranchLengthsHyperLogPrior() const {
        return LambdaHyperLogPrior() + BranchLengthsHyperInvShapeLogPrior() +
               branchlength->GetLogProb();
    }

    double GeneNucRatesHyperLogPrior() const {
        double total = 0;
        total -= nucrelratehyperinvconc;
        total -= nucstathyperinvconc;
        return total;
    }

    double ShiftProbHyperLogPrior() const {
        double total = 0;
        if (pihyperinvconc) {
            for (int k = 1; k < Ncond; k++) {
                total += GetPiLogPrior(k);
            }
        }
        if (shiftmode) {
            for (int k = 1; k < Ncond; k++) {
                total += GetShiftProbHyperMeanLogPrior(k);
                total += GetShiftProbHyperInvConcLogPrior(k);
            }
        }
        return total;
    }

    double GetPiLogPrior(int k) const {
        if ((pi[k - 1] < minpi) || (pi[k - 1] > maxpi)) {
            return log(0);
        }
        double alpha = pihypermean / pihyperinvconc;
        double beta = (1 - pihypermean) / pihyperinvconc;
        return Random::logBetaDensity(pi[k - 1], alpha, beta);
    }

    double GetShiftProbHyperMeanLogPrior(int k) const {
        if (shiftprobhypermean[k - 1] < minshiftprobhypermean) {
            return log(0);
            // return Random::INFPROB;
        }
        return 0;
    }

    double GetShiftProbHyperInvConcLogPrior(int k) const {
        if (shiftprobhyperinvconc[k - 1] > maxshiftprobhyperinvconc) {
            return log(0);
        }
        return 0;
        /*
        double alpha = shiftprobhypermean[k - 1] / shiftprobhyperinvconc[k - 1];
        double beta = (1 - shiftprobhypermean[k - 1]) / shiftprobhyperinvconc[k - 1];
        if ((alpha < 1.0) || (beta < 1.0)) {
            return log(0);
            return Random::INFPROB;
        }
        return -10*shiftprobhyperinvconc[k - 1];
        */
    }

    double GetCountLogProb(int k) const { return shiftcountarray->GetMarginalLogProb(k - 1); }

    double GetCountLogProb() const {
        double total = 0;
        for (int k = 1; k < Ncond; k++) {
            total += GetCountLogProb(k);
        }
        return total;
    }

    double GetLogLikelihood() const { return lnL; }

    //-------------------
    // Suff Stat Log Probs
    //-------------------

    // suff stat for global branch lengths, as a function of lambda
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    // suff stat for gene-specific branch lengths, as a function of bl
    // hyperparameters
    double BranchLengthsHyperSuffStatLogProb() const {
        return lengthhypersuffstatarray->GetLogProb(*branchlength, blhyperinvshape);
    }

    // suff stat for gene-specific nuc rates, as a function of nucrate
    // hyperparameters
    double NucRatesHyperSuffStatLogProb() const {
        double total = 0;
        total += nucrelratesuffstat.GetLogProb(nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total += nucstatsuffstat.GetLogProb(nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // logprob for moving lambda
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    // logprob for moving hyperparameters of gene-specific branchlengths
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperInvShapeLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates hyper params
    double NucRatesHyperLogProb() const {
        return GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
    }

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
        int nrep = 10;

        for (int rep = 0; rep < nrep; rep++) {
            if (withtoggle) {
                MasterReceiveShiftCounts();
                movechrono.Start();
                MoveShiftProbHyperParameters(100);
                movechrono.Stop();
                MasterSendShiftProbHyperParameters();
            }

            if (blmode == 1) {
                MasterReceiveBranchLengthsHyperSuffStat();
                movechrono.Start();
                MoveBranchLengthsHyperParameters(10);
                movechrono.Stop();
                MasterSendBranchLengthsHyperParameters();
            }

            if (nucmode == 1) {
                MasterReceiveNucRatesHyperSuffStat();
                movechrono.Start();
                MoveNucRatesHyperParameters(10);
                movechrono.Stop();
                MasterSendNucRatesHyperParameters();
            }
        }

        MasterReceiveGeneBranchLengths();
        MasterReceiveGeneNucRates();
        if (Ncond == 1) {
            MasterReceivePredictedDNDS();
        }
        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {
        movechrono.Start();
        mapchrono.Start();
        GeneResampleSub(1.0);
        mapchrono.Stop();
        movechrono.Stop();

        int nrep = 10;

        for (int rep = 0; rep < nrep; rep++) {
            movechrono.Start();
            GeneMove(1, 3);
            movechrono.Stop();

            if (withtoggle) {
                SlaveSendShiftCounts();
                SlaveReceiveShiftProbHyperParameters();
                if (shiftprobinvconc) {
                    movechrono.Start();
                    GeneResampleShiftProbs();
                    movechrono.Stop();
                }
            }

            if (blmode == 1) {
                SlaveSendBranchLengthsHyperSuffStat();
                SlaveReceiveBranchLengthsHyperParameters();
            }

            if (nucmode == 1) {
                SlaveSendNucRatesHyperSuffStat();
                SlaveReceiveNucRatesHyperParameters();
            }
        }

        SlaveSendGeneBranchLengths();
        SlaveSendGeneNucRates();
        if (Ncond == 1) {
            SlaveSendPredictedDNDS();
        }
        SlaveSendLogProbs();
    }

    void GeneCollectShiftCounts() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            for (int k = 1; k < Ncond; k++) {
                (*shiftcountarray)[gene][k - 1] = geneprocess[gene]->GetNShift(k);
                (*totcount)[gene] = geneprocess[gene]->GetNTarget();
            }
        }
    }

    void GeneResampleShiftProbs() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleShiftProb();
        }
    }

    void GeneMove(int nrep0, int nrep) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveParameters(nrep0, nrep);
        }
    }

    void GeneResampleSub(double frac) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneDiffSelDoublySparseModel::LambdaHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneDiffSelDoublySparseModel::LambdaHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters(int nrep) {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, nrep);
            BranchLengthsHyperScalingMove(0.3, nrep);
        }

        ScalingMove(blhyperinvshape, 1.0, nrep,
                    &MultiGeneDiffSelDoublySparseModel::BranchLengthsHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, nrep,
                    &MultiGeneDiffSelDoublySparseModel::BranchLengthsHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.1, nrep,
                    &MultiGeneDiffSelDoublySparseModel::BranchLengthsHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);

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

    void MoveNucRatesHyperParameters(int nrep) {
        ProfileMove(nucrelratehypercenter, 1.0, 1, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, nrep,
                    &MultiGeneDiffSelDoublySparseModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelDoublySparseModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveShiftProbHyperParameters(int nrep) {
        if (pihyperinvconc) {
            MovePi(0.3, nrep);
        }
        if (shiftmode) {
            MoveShiftProbHyperMean(0.3, nrep);
            MoveShiftProbHyperInvConc(0.3, nrep);
        }
    }

    double MovePi(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int k = 1; k < Ncond; k++) {
            for (int rep = 0; rep < nrep; rep++) {
                double deltalogprob = -GetPiLogPrior(k) - GetCountLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double bk = pi[k - 1];
                pi[k - 1] += m;
                while ((pi[k - 1] < minpi) || (pi[k - 1] > maxpi)) {
                    if (pi[k - 1] < minpi) {
                        pi[k - 1] = 2 * minpi - pi[k - 1];
                    }
                    if (pi[k - 1] > maxpi) {
                        pi[k - 1] = 2 * maxpi - pi[k - 1];
                    }
                }
                deltalogprob += GetPiLogPrior(k) + GetCountLogProb(k);

                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    pi[k - 1] = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveShiftProbHyperMean(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int k = 1; k < Ncond; k++) {
            for (int rep = 0; rep < nrep; rep++) {
                double deltalogprob = -GetShiftProbHyperMeanLogPrior(k) - GetCountLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double bk = shiftprobhypermean[k - 1];
                shiftprobhypermean[k - 1] += m;
                while ((shiftprobhypermean[k - 1] < minshiftprobhypermean) ||
                       (shiftprobhypermean[k - 1] > maxshiftprobhypermean)) {
                    if (shiftprobhypermean[k - 1] < minshiftprobhypermean) {
                        shiftprobhypermean[k - 1] =
                            2 * minshiftprobhypermean - shiftprobhypermean[k - 1];
                    }
                    if (shiftprobhypermean[k - 1] > maxshiftprobhypermean) {
                        shiftprobhypermean[k - 1] =
                            2 * maxshiftprobhypermean - shiftprobhypermean[k - 1];
                    }
                }
                deltalogprob += GetShiftProbHyperMeanLogPrior(k) + GetCountLogProb(k);

                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    shiftprobhypermean[k - 1] = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveShiftProbHyperInvConc(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int k = 1; k < Ncond; k++) {
            for (int rep = 0; rep < nrep; rep++) {
                double deltalogprob = -GetShiftProbHyperInvConcLogPrior(k) - GetCountLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                double bk = shiftprobhyperinvconc[k - 1];
                shiftprobhyperinvconc[k - 1] *= e;
                deltalogprob += GetShiftProbHyperInvConcLogPrior(k) + GetCountLogProb(k);
                deltalogprob += m;

                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    shiftprobhyperinvconc[k - 1] = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    //-------------------
    // MPI send / receive
    //-------------------

    void MasterSendBranchLengthsHyperParameters() {
        MasterSendGlobal(*branchlength, blhyperinvshape);
    }

    void SlaveReceiveBranchLengthsHyperParameters() {
        SlaveReceiveGlobal(*branchlength, blhyperinvshape);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
        }
    }

    void MasterSendGeneBranchLengths() { MasterSendGeneArray(*branchlengtharray); }

    void SlaveReceiveGeneBranchLengths() {
        SlaveReceiveGeneArray(*branchlengtharray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengths(branchlengtharray->GetVal(gene));
        }
    }

    void SlaveSendGeneBranchLengths() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
        SlaveSendGeneArray(*branchlengtharray);
    }

    void MasterReceiveGeneBranchLengths() { MasterReceiveGeneArray(*branchlengtharray); }

    void MasterSendGeneNucRates() { MasterSendGeneArray(*nucrelratearray, *nucstatarray); }

    void SlaveReceiveGeneNucRates() {
        SlaveReceiveGeneArray(*nucrelratearray, *nucstatarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates(nucrelratearray->GetVal(gene),
                                           nucstatarray->GetVal(gene));
        }
    }

    void SlaveSendGeneNucRates() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->GetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
        }
        SlaveSendGeneArray(*nucrelratearray, *nucstatarray);
    }

    void MasterReceiveGeneNucRates() { MasterReceiveGeneArray(*nucrelratearray, *nucstatarray); }

    void MasterSendNucRatesHyperParameters() {
        MasterSendGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        MasterSendGlobal(nucstathypercenter, nucstathyperinvconc);
    }

    void SlaveReceiveNucRatesHyperParameters() {
        SlaveReceiveGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        SlaveReceiveGlobal(nucstathypercenter, nucstathyperinvconc);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter,
                                                          nucrelratehyperinvconc,
                                                          nucstathypercenter, nucstathyperinvconc);
        }
    }

    // log probs

    void SlaveSendLogProbs() {
        GeneLogPrior = 0;
        lnL = 0;
        meanWidth = 0;
        meanEps = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
            meanWidth += geneprocess[gene]->GetMeanWidth();
            meanEps += geneprocess[gene]->GetMaskEpsilon();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
        SlaveSendAdditive(meanWidth);
        SlaveSendAdditive(meanEps);

        moveTime = movechrono.GetTime();
        mapTime = mapchrono.GetTime();
        SlaveSendAdditive(moveTime);
        SlaveSendAdditive(mapTime);
    }

    void MasterReceiveLogProbs() {
        GeneLogPrior = 0;
        MasterReceiveAdditive(GeneLogPrior);
        lnL = 0;
        MasterReceiveAdditive(lnL);
        meanWidth = 0;
        MasterReceiveAdditive(meanWidth);
        meanWidth /= GetLocalNgene();
        meanEps = 0;
        MasterReceiveAdditive(meanEps);
        meanEps /= GetLocalNgene();

        moveTime = 0;
        mapTime = 0;
        MasterReceiveAdditive(moveTime);
        MasterReceiveAdditive(mapTime);
        moveTime /= (GetNprocs() - 1);
        mapTime /= (GetNprocs() - 1);
    }

    void MasterReceiveShiftCounts() {
        MasterReceiveGeneArray(*shiftcountarray);
        MasterReceiveGeneArray(*totcount);
    }

    void SlaveSendShiftCounts() {
        GeneCollectShiftCounts();
        SlaveSendGeneArray(*shiftcountarray);
        SlaveSendGeneArray(*totcount);
    }

    void MasterSendShiftProbHyperParameters() {
        MasterSendGlobal(shiftprobhypermean, shiftprobhyperinvconc);
        MasterSendGlobal(pi);
    }

    void SlaveReceiveShiftProbHyperParameters() {
        SlaveReceiveGlobal(shiftprobhypermean, shiftprobhyperinvconc);
        SlaveReceiveGlobal(pi);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetShiftProbHyperParameters(pi, shiftprobhypermean,
                                                           shiftprobhyperinvconc);
        }
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

    void MasterTraceSiteStats(string name, int mode) {
        MasterReceiveGeneArray(*shiftprobarray);
        for (int k = 1; k < Ncond; k++) {
            ostringstream s;
            s << name << "_" << k << ".geneshiftprob";
            ofstream os(s.str().c_str(), ios_base::app);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                os << shiftprobarray->GetVal(gene)[k - 1] << '\t';
            }
            os << '\n';
        }

        MasterReceiveShiftCounts();
        for (int k = 1; k < Ncond; k++) {
            ostringstream s;
            s << name << "_" << k << ".geneshiftcounts";
            ofstream os(s.str().c_str(), ios_base::app);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                os << shiftcountarray->GetVal(gene)[k - 1] << '/' << totcount->GetVal(gene) << '\t';
            }
            os << '\n';
        }
        ofstream os((name + ".genemaskcounts").c_str(), ios_base::app);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            os << totcount->GetVal(gene) << '\t';
        }
        os << '\n';

        if (mode == 2) {
            for (int k = 0; k < Ncond; k++) {
                ostringstream s;
                s << name << "_" << k << ".fitness";
                ofstream os(s.str().c_str(), ios_base::app);
                for (int proc = 1; proc < GetNprocs(); proc++) {
                    int totnsite = GetSlaveTotNsite(proc);
                    double *array = new double[totnsite * Naa];
                    MPI_Status stat;
                    MPI_Recv(array, totnsite * Naa, MPI_DOUBLE, proc, TAG1, MPI_COMM_WORLD, &stat);

                    int i = 0;
                    for (int gene = 0; gene < Ngene; gene++) {
                        if (GeneAlloc[gene] == proc) {
                            os << GeneName[gene] << '\t';
                            int nsite = GeneNsite[gene];
                            for (int j = 0; j < nsite; j++) {
                                for (int a = 0; a < Naa; a++) {
                                    os << array[i++] << '\t';
                                }
                            }
                        }
                    }
                    if (i != totnsite * Naa) {
                        cerr << "error in "
                                "MultiGeneDiffSelDoublySparseModel::MasterTraceSiteStats: "
                                "non matching number of sites\n";
                        exit(1);
                    }
                    delete[] array;
                }
                os << '\n';
                os.flush();
            }

            for (int k = 1; k < Ncond; k++) {
                ostringstream s;
                s << name << "_" << k << ".shifttoggle";
                ofstream os(s.str().c_str(), ios_base::app);
                for (int proc = 1; proc < GetNprocs(); proc++) {
                    int totnsite = GetSlaveTotNsite(proc);
                    int *array = new int[totnsite * Naa];
                    MPI_Status stat;
                    MPI_Recv(array, totnsite * Naa, MPI_INT, proc, TAG1, MPI_COMM_WORLD, &stat);

                    int i = 0;
                    for (int gene = 0; gene < Ngene; gene++) {
                        if (GeneAlloc[gene] == proc) {
                            os << GeneName[gene] << '\t';
                            int nsite = GeneNsite[gene];
                            for (int j = 0; j < nsite; j++) {
                                for (int a = 0; a < Naa; a++) {
                                    os << array[i++] << '\t';
                                }
                            }
                        }
                    }
                    if (i != totnsite * Naa) {
                        cerr << "error in "
                                "MultiGeneDiffSelDoublySparseModel::MasterTraceSiteStats: "
                                "non matching number of sites\n";
                        exit(1);
                    }
                    delete[] array;
                }
                os << '\n';
                os.flush();
            }
        }
    }

    void SlaveTraceSiteStats(int mode) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            (*shiftprobarray)[gene] = geneprocess[gene]->GetShiftProbVector();
        }
        SlaveSendGeneArray(*shiftprobarray);
        SlaveSendShiftCounts();

        if (mode == 2) {
            for (int k = 0; k < Ncond; k++) {
                int ngene = GetLocalNgene();
                int totnsite = GetLocalTotNsite();
                double *array = new double[totnsite * Naa];
                int i = 0;
                for (int gene = 0; gene < ngene; gene++) {
                    geneprocess[gene]->GetFitnessArray(k, array + i);
                    i += GetLocalGeneNsite(gene) * Naa;
                }
                if (i != totnsite * Naa) {
                    cerr << "error in "
                            "MultiGeneDiffSelDoublySparseModel::SlaveTraceSiteStats: non "
                            "matching number of sites\n";
                    exit(1);
                }

                MPI_Send(array, totnsite * Naa, MPI_DOUBLE, 0, TAG1, MPI_COMM_WORLD);
                delete[] array;
            }
            for (int k = 1; k < Ncond; k++) {
                int ngene = GetLocalNgene();
                int totnsite = GetLocalTotNsite();
                int *array = new int[totnsite * Naa];
                int i = 0;
                for (int gene = 0; gene < ngene; gene++) {
                    geneprocess[gene]->GetShiftToggleArray(k, array + i);
                    i += GetLocalGeneNsite(gene) * Naa;
                }
                if (i != totnsite * Naa) {
                    cerr << "error in MultiGeneCodonM2aModel::SlaveTraceSitesPostProb: "
                            "non matching number of sites\n";
                    exit(1);
                }

                MPI_Send(array, totnsite * Naa, MPI_INT, 0, TAG1, MPI_COMM_WORLD);
                delete[] array;
            }
        }
    }
};
