
#include "Chrono.hpp"
#include "DiffSelModel.hpp"
#include "IIDMultiBernBeta.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"

class MultiGeneDiffSelModel : public MultiGeneProbModel {
  private:
    // const double minshiftprobhypermean = 0.01;

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

    // each gene defines its own DiffSelModel
    std::vector<DiffSelModel *> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    int blmode;
    int nucmode;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneDiffSelModel(string datafile, string intreefile, int inNcond, int inNlevel,
                                int incodonmodel, int inblmode, int innucmode, int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs), nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {

        codonmodel = incodonmodel;
        blmode = inblmode;
        nucmode = innucmode;

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

        lnL = 0;
        GeneLogPrior = 0;

        if (!GetMyid()) {
            geneprocess.assign(0, (DiffSelModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (DiffSelModel *)0);

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                int fixglob = 1;
                int fixvar = 1;
                geneprocess[gene] = new DiffSelModel(GetLocalGeneName(gene), treefile, Ncond,
                                                           Nlevel, fixglob, fixvar, codonmodel);
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
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
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {
        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveGeneBranchLengths();
        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveGeneNucRates();
        GeneUpdate();
        SlaveSendLogProbs();
    }

    void GeneUpdate() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->Update();
        }
    }

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    int GetNbranch() const { return tree->GetNbranch(); }

    const Tree *GetTree() const { return tree; }

    //-------------------
    // Traces and Monitors
    //-------------------

    double GetMeanLength() const { return branchlengtharray->GetMeanLength(); }

    double GetVarLength() const { return branchlengtharray->GetVarLength(); }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "stdev\t";
        os << '\n';
    }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << GetMeanLength() << '\t' << GetVarLength() << '\t';
        os << '\n';
        os.flush();
    }

    void Monitor(ostream &os) const override {}

    void MasterToStream(ostream &os) const override {
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

    void MasterFromStream(istream &is) override {
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
        int nrep = 1;
        // int nrep = 3;

        for (int rep = 0; rep < nrep; rep++) {

            if (blmode == 1) {
                MasterReceiveBranchLengthsHyperSuffStat();
                movechrono.Start();
                MoveBranchLengthsHyperParameters();
                movechrono.Stop();
                MasterSendBranchLengthsHyperParameters();
            }

            if (nucmode == 1) {
                MasterReceiveNucRatesHyperSuffStat();
                movechrono.Start();
                MoveNucRatesHyperParameters();
                movechrono.Stop();
                MasterSendNucRatesHyperParameters();
            }
        }

        // collect current state
        MasterReceiveGeneBranchLengths();
        MasterReceiveGeneNucRates();
        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {
        movechrono.Start();
        mapchrono.Start();
        GeneResampleSub(1.0);
        mapchrono.Stop();
        movechrono.Stop();

        int nrep = 1;
        // int nrep = 3;

        for (int rep = 0; rep < nrep; rep++) {
            movechrono.Start();
            GeneMove();
            movechrono.Stop();

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
        SlaveSendLogProbs();
    }

    void GeneMove() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveParameters(1, 10);
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
            geneprocess[gene]->GetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
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
        ScalingMove(lambda, 1.0, 10, &MultiGeneDiffSelModel::LambdaHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneDiffSelModel::LambdaHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, 10);
            BranchLengthsHyperScalingMove(0.3, 10);
        }

        ScalingMove(blhyperinvshape, 1.0, 10,
                    &MultiGeneDiffSelModel::BranchLengthsHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10,
                    &MultiGeneDiffSelModel::BranchLengthsHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);

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

    void MoveNucRatesHyperParameters() {
        ProfileMove(nucrelratehypercenter, 1.0, 1, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10,
                    &MultiGeneDiffSelModel::NucRatesHyperLogProb,
                    &MultiGeneDiffSelModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
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
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);

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
        moveTime = 0;
        mapTime = 0;
        MasterReceiveAdditive(moveTime);
        MasterReceiveAdditive(mapTime);
        moveTime /= (GetNprocs() - 1);
        mapTime /= (GetNprocs() - 1);
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
        for (int k = 0; k < Ncond; k++) {
            ostringstream s;
            if (! k)    {
                s << name << ".fitness";
            }
            else    {
                s << name << "_" << k << ".delta";
            }
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
                            "MultiGeneDiffSelModel::MasterTraceSiteStats: non "
                            "matching number of sites\n";
                    exit(1);
                }
                delete[] array;
            }
            os << '\n';
            os.flush();
        }
    }

    void SlaveTraceSiteStats(int mode) {
        for (int k = 0; k < Ncond; k++) {
            int ngene = GetLocalNgene();
            int totnsite = GetLocalTotNsite();
            double *array = new double[totnsite * Naa];
            int i = 0;
            for (int gene = 0; gene < ngene; gene++) {
                if (! k)    {
                    geneprocess[gene]->GetBaselineArray(array + i);
                }
                else    {
                    geneprocess[gene]->GetDeltaArray(k, array + i);
                }
                i += GetLocalGeneNsite(gene) * Naa;
            }
            if (i != totnsite * Naa) {
                cerr << "error in MultiGeneDiffSelModel::SlaveTraceSiteStats: "
                        "non matching number of sites\n";
                exit(1);
            }

            MPI_Send(array, totnsite * Naa, MPI_DOUBLE, 0, TAG1, MPI_COMM_WORLD);
            delete[] array;
        }
    }
};
