
#include "AAMutSelDSBDPOmegaModel.hpp"
#include "IIDBernoulliGamma.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "Permutation.hpp"

class MultiGeneAAMutSelDSBDPOmegaModel : public MultiGeneProbModel {
  private:
    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;

    string treefile;

    int Ntaxa;
    int Nbranch;

    int basemin;
    int baseNcat;
    int Ncat;

    // branch lengths
    double lambda;
    BranchIIDGamma *branchlength;
    GammaSuffStat hyperlengthsuffstat;

    double blhyperinvshape;
    GammaWhiteNoiseArray *branchlengtharray;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
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

    // omega*: iid gamma across genes
    double omegahypermean;
    double omegahyperinvshape;
    IIDGamma *omegaarray;
    GammaSuffStat omegahypersuffstat;
    SimpleArray<double>* genednds;

    double dposompihypermean;
    double dposompihyperinvconc;
    double dposompi;
    double dposomhypermean;
    double dposomhyperinvshape;
    IIDBernoulliGamma *dposomarray;
    BernoulliGammaSuffStat dposomhypersuffstat;

    // truncated stick-breaking mixture of baseNcat Dirichlet densities
    // centers
    vector<double> basecenterhypercenter;
    double basecenterhyperinvconc;
    IIDDirichlet *basecenterarray;
    // concentrations
    double baseconchypermean;
    double baseconchyperinvshape;
    IIDGamma *baseconcentrationarray;
    // and associated suffstatarray
    DirichletSuffStatArray *basesuffstatarray;
    // weights:
    double basekappa;
    StickBreakingProcess *baseweight;
    // and associated occupancy suffstat
    OccupancySuffStat *baseoccupancy;
    // for keeping track of label-switching moves done by master
    // and sending the resulting permutation to slaves
    Permutation *permutocc;

    std::vector<AAMutSelDSBDPOmegaModel *> geneprocess;

    double lnL;
    double GeneLogPrior;
    double MeanNcluster;
    double MeanBaseNcluster;
    double MeanStatEnt;
    double MeanAAConc;
    double MeanAACenterEnt;
    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    int blmode, nucmode, basemode, omegamode, omegaprior, modalprior;

    Chrono totchrono;
    Chrono paramchrono;
    Chrono basechrono;
    Chrono blchrono;
    Chrono aachrono;

    int burnin;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneAAMutSelDSBDPOmegaModel(string datafile, string intreefile, int inNcat, int inbaseNcat,
                                     int inblmode, int innucmode, int inbasemode, int inomegamode,
                                     int inomegaprior, int inmodalprior, double indposompihypermean,
                                     double indposompihyperinvconc, int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs), nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {
        AllocateAlignments(datafile);
        treefile = intreefile;
        Ncat = inNcat;

        burnin = 0;

        basemin = 0;
        if (inbaseNcat < 0) {
            basemin = 1;
            baseNcat = -inbaseNcat;
            if (baseNcat != 2) {
                cerr << "error in basencat\n";
                exit(1);
            }
        } else {
            baseNcat = inbaseNcat;
        }

        blmode = inblmode;
        nucmode = innucmode;
        basemode = inbasemode;
        omegamode = inomegamode;
        omegaprior = inomegaprior;
        modalprior = inmodalprior;
        dposompihypermean = indposompihypermean;
        dposompihyperinvconc = indposompihyperinvconc;

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
    }

    void Allocate() {
        lambda = 10;
        branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
        blhyperinvshape = 0.1;
        if (blmode == 2) {
            lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
            lengthhypersuffstatarray = 0;
            branchlengtharray = 0;
        } else {
            // when blmode == 1:
            // blarray ~ GammaWN(bl,invshape)
            // but then bl's are iid ~ Gamma(1.0,lambda) (thus, not necessarily all
            // equal)
            //
            // when blmode == 0
            // each gene has its own lambda
            branchlength->SetAllBranches(1.0 / lambda);
            branchlengtharray = new GammaWhiteNoiseArray(GetLocalNgene(), *tree, *branchlength,
                                                         1.0 / blhyperinvshape);
            lengthpathsuffstatarray = 0;
            lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
        }

        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 0.1 / Nrr;

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 0.1 / Nnuc;

        nucrelratearray =
            new IIDDirichlet(GetLocalNgene(), nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        nucstatarray =
            new IIDDirichlet(GetLocalNgene(), nucstathypercenter, 1.0 / nucstathyperinvconc);

        omegahypermean = 1.0;
        omegahyperinvshape = 0.5;
        dposompi = dposompihypermean;
        dposomhypermean = 1.0;
        dposomhyperinvshape = 0.5;
        if (omegaprior == 0) {
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            omegaarray = new IIDGamma(GetLocalNgene(), alpha, beta);
            if (omegamode == 3) {
                for (int i = 0; i < GetLocalNgene(); i++) {
                    (*omegaarray)[i] = 1.0;
                }
            }
        } else if (omegaprior == 1) {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            dposomarray = new IIDBernoulliGamma(GetLocalNgene(), dposompi, alpha, beta);
            for (int i = 0; i < GetLocalNgene(); i++) {
                (*dposomarray)[i] = 0;
            }
        } else {
            cerr << "error: unrecognized omega prior\n";
            exit(1);
        }

        genednds = new SimpleArray<double>(GetLocalNgene());

        lnL = 0;
        GeneLogPrior = 0;
        MeanStatEnt = 0;
        MeanAAConc = 0;
        MeanAACenterEnt = 0;

        basekappa = 1.0;
        baseweight = new StickBreakingProcess(baseNcat, basekappa);
        baseoccupancy = new OccupancySuffStat(baseNcat);
        permutocc = new Permutation(baseNcat);

        basecenterhypercenter.assign(Naa, 1.0 / Naa);
        basecenterhyperinvconc = 1.0 / Naa;

        basecenterarray =
            new IIDDirichlet(baseNcat, basecenterhypercenter, 1.0 / basecenterhyperinvconc);
        basecenterarray->SetUniform();

        baseconchypermean = Naa;
        baseconchyperinvshape = 1.0;
        double alpha = 1.0 / baseconchyperinvshape;
        double beta = alpha / baseconchypermean;

        baseconcentrationarray = new IIDGamma(baseNcat, alpha, beta);
        for (int k = 0; k < baseNcat; k++) {
            (*baseconcentrationarray)[k] = 20.0;
        }
        if (basemin == 1) {
            (*baseconcentrationarray)[0] = 1.0;
        }

        // suff stats for component aa fitness arrays
        basesuffstatarray = new DirichletSuffStatArray(baseNcat, Naa);

        if (!GetMyid()) {
            geneprocess.assign(0, (AAMutSelDSBDPOmegaModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (AAMutSelDSBDPOmegaModel *)0);

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene] = new AAMutSelDSBDPOmegaModel(
                    GetLocalGeneName(gene), treefile, omegamode, omegaprior, Ncat, baseNcat);
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->SetBaseMode(basemode);
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

        if (omegaprior == 0) {
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            omegaarray->SetShape(alpha);
            omegaarray->SetScale(beta);
        } else {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            dposomarray->SetPi(dposompi);
            dposomarray->SetShape(alpha);
            dposomarray->SetScale(beta);
        }
    }

    void MasterUpdate() override {
        FastUpdate();

        if (nprocs > 1) {
            if (blmode >= 2) {
                MasterSendGlobalBranchLengths();
            } else {
                MasterSendBranchLengthsHyperParameters();
                MasterSendGeneBranchLengths();
            }

            MasterSendNucRatesHyperParameters();
            MasterSendGeneNucRates();

            if (basemode >= 2) {
                MasterSendBaseMixture();
            }

            if (omegamode == 1) {
                MasterSendOmegaHyperParameters();
                MasterSendOmega();
            }
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {
        if (blmode >= 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveBranchLengthsHyperParameters();
            SlaveReceiveGeneBranchLengths();
        }

        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveGeneNucRates();

        if (basemode >= 2) {
            SlaveReceiveBaseMixture();
        }
        if (omegamode == 1) {
            SlaveReceiveOmegaHyperParameters();
            SlaveReceiveOmega();
        }
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
            if (blmode >= 2) {
                MasterSendGlobalBranchLengths();
            } else {
                MasterSendBranchLengthsHyperParameters();
                MasterSendGeneBranchLengths();
            }

            MasterSendNucRatesHyperParameters();
            MasterSendGeneNucRates();

            if (basemode >= 2) {
                MasterSendBaseMixture();
            }

            if (omegamode == 1) {
                MasterSendOmegaHyperParameters();
                MasterSendOmega();
            }
            // MasterReceiveLogProbs();
        }
    }

    void SlavePostPred(string name) override {
        if (blmode >= 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveBranchLengthsHyperParameters();
            SlaveReceiveGeneBranchLengths();
        }

        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveGeneNucRates();

        if (basemode >= 2) {
            SlaveReceiveBaseMixture();
        }
        if (omegamode == 1) {
            SlaveReceiveOmegaHyperParameters();
            SlaveReceiveOmega();
        }
        GenePostPred(name);
        // SlaveSendLogProbs();
    }

    void GenePostPred(string name) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->PostPred(name + GetLocalGeneName(gene));
        }
    }

    void TracePredictedDNDS(ostream& os) const   {
        for (int gene = 0; gene < Ngene; gene++) {
            os << genednds->GetVal(gene) << '\t';
        }
        os << '\n';
    }

    void MasterReceivePredictedDNDS() {
        MasterReceiveGeneArray(*genednds);
    }

    void SlaveSendPredictedDNDS() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            (*genednds)[gene] = geneprocess[gene]->GetPredictedDNDS();
        }
        SlaveSendGeneArray(*genednds);
    }

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    //-------------------
    // Traces and Monitors
    //-------------------

    double GetSlaveMoveTime() const { return moveTime; }

    double GetSlaveMapTime() const { return mapTime; }

    double GetMasterMoveTime() const { return movechrono.GetTime(); }

    int GetBaseNcluster() const {
        int n = 0;
        for (int i = 0; i < baseNcat; i++) {
            if (baseoccupancy->GetVal(i)) {
                n++;
            }
        }
        return n;
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\t";
        os << "length\t";
        if (blmode < 2) {
            os << "stdev\t";
        }
        if (omegamode != 3) {
            if (omegaprior == 0) {
                os << "meanomega\t";
                os << "varomega\t";
            } else {
                os << "npos\tposmean\tdposom_pi\tmeandposom\tinvshape\t";
            }
        }
        if (Ncat > 1) {
            os << "ncluster\t";
        }
        if (basemode >= 2) {
            os << "basencluster\t";
            os << "basekappa\t";
        }
        os << "aastatent\t";
        os << "baseconc\t";
        os << "baseent\t";
        os << "nucrr\tinvconc\t";
        os << "nucstat\tinvconc\n";
    }

    double GetMeanTotalLength() const {
        double tot = 0;
        for (int j = 0; j < Nbranch; j++) {
            tot += branchlength->GetVal(j);
        }
        return tot;
    }

    double GetMeanLength() const { return branchlengtharray->GetMeanLength(); }

    double GetVarLength() const { return branchlengtharray->GetVarLength(); }

    void Trace(ostream &os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        if (blmode == 2) {
            os << GetMeanTotalLength() << '\t';
        } else {
            os << GetMeanLength() << '\t';
            os << sqrt(GetVarLength()) << '\t';
        }
        if (omegamode != 3) {
            if (omegaprior == 0) {
                os << omegaarray->GetMean() << '\t';
                os << omegaarray->GetVar() << '\t';
            } else {
                os << dposomarray->GetNpos() << '\t';
                os << dposomarray->GetPosMean() << '\t';
                os << dposompi << '\t' << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
            }
        }
        if (Ncat > 1) {
            os << MeanNcluster << '\t';
        }
        if (basemode >= 2) {
            os << GetBaseNcluster() << '\t';
            os << basekappa << '\t';
        }
        os << MeanStatEnt << '\t';
        os << MeanAAConc << '\t';
        os << MeanAACenterEnt << '\t';
        os << Random::GetEntropy(nucrelratehypercenter) << '\t' << nucrelratehyperinvconc << '\t';
        os << Random::GetEntropy(nucstathypercenter) << '\t' << nucstathyperinvconc << '\n';
        os.flush();
    }

    void TraceOmega(ostream &os) const {
        if (omegaprior == 0) {
            for (int gene = 0; gene < Ngene; gene++) {
                os << omegaarray->GetVal(gene) << '\t';
            }
        } else {
            for (int gene = 0; gene < Ngene; gene++) {
                os << 1.0 + dposomarray->GetVal(gene) << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

    void TraceMixture(ostream &os) const {
        for (int k = 0; k < baseNcat; k++) {
            os << baseweight->GetVal(k) << '\t' << baseconcentrationarray->GetVal(k) << '\t';
            for (int a = 0; a < Naa; a++) {
                os << (int)(100 * basecenterarray->GetVal(k)[a]) << '\t';
            }
            os << '\n';
        }
        os.flush();
    }

    void Monitor(ostream &os) const override {
        os << totchrono.GetTime() << '\t' << paramchrono.GetTime() << '\t' << basechrono.GetTime()
           << '\t' << blchrono.GetTime() << '\t' << aachrono.GetTime() << '\n';
        os << "prop time in param moves: " << paramchrono.GetTime() / totchrono.GetTime() << '\n';
        os << "sub prop time in base moves : " << basechrono.GetTime() / paramchrono.GetTime()
           << '\n';
        os << "sub prop time in bl moves   : " << blchrono.GetTime() / paramchrono.GetTime()
           << '\n';
        os << "sub prop time in aa moves   : " << aachrono.GetTime() / paramchrono.GetTime()
           << '\n';
    }

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

        is >> nucrelratehypercenter;
        is >> nucrelratehyperinvconc;
        is >> nucstathypercenter;
        is >> nucstathyperinvconc;
        is >> *nucrelratearray;
        is >> *nucstatarray;

        if (omegamode != 3) {
            if (omegaprior == 0) {
                is >> omegahypermean;
                is >> omegahyperinvshape;
                is >> *omegaarray;
            } else {
                is >> dposompi;
                is >> dposomhypermean;
                is >> dposomhyperinvshape;
                is >> *dposomarray;
            }
        }
        if (basemode == 2) {
            is >> *basecenterarray;
            is >> *baseconcentrationarray;
            is >> basekappa;
            baseweight->FromStreamSB(is);
        }

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

        os << nucrelratehypercenter << '\t';
        os << nucrelratehyperinvconc << '\t';
        os << nucstathypercenter << '\t';
        os << nucstathyperinvconc << '\t';
        os << *nucrelratearray << '\t';
        os << *nucstatarray << '\t';

        if (omegamode != 3) {
            if (omegaprior == 0) {
                os << omegahypermean << '\t';
                os << omegahyperinvshape << '\t';
                os << *omegaarray << '\t';
            } else {
                os << dposompi << '\t';
                os << dposomhypermean << '\t';
                os << dposomhyperinvshape << '\t';
                os << *dposomarray << '\t';
            }
        }
        if (basemode == 2) {
            os << *basecenterarray << '\t';
            os << *baseconcentrationarray << '\t';
            os << basekappa << '\t';
            baseweight->ToStreamSB(os);
        }

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

    void PrintBaseMixtureLogo(ostream &os) const {
        os << baseNcat << '\t' << Naa << '\n';
        for (int k = 0; k < baseNcat; k++) {
            os << baseweight->GetVal(k) << '\t';
            for (int a = 0; a < Naa; a++) {
                os << basecenterarray->GetVal(k)[a] << '\t';
            }
            os << '\n';
        }
    }

    void PrintBaseSampleLogo(ostream &os) const {
        int nrep = 8;
        vector<double> x(Naa);
        os << baseNcat * nrep << '\t' << Naa << '\n';
        for (int k = 0; k < baseNcat; k++) {
            for (int rep = 0; rep < nrep; rep++) {
                Random::DirichletSample(x, basecenterarray->GetVal(k),
                                        baseconcentrationarray->GetVal(k));
                os << 1 << '\t';
                for (int a = 0; a < Naa; a++) {
                    os << x[a] << '\t';
                }
                os << '\n';
            }
        }
    }

    //-------------------
    // Updates
    //-------------------

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------

    double GetLogPrior() const {
        double total = GeneLogPrior;

        if (blmode == 2) {
            total += GlobalBranchLengthsLogPrior();
        } else if (blmode == 1) {
            total += GeneBranchLengthsHyperLogPrior();
        }

        if (nucmode == 1) {
            total += GeneNucRatesHyperLogPrior();
        }

        if (omegamode != 3) {
            total += OmegaHyperLogPrior();
            total += OmegaLogPrior();
        }
        total += BaseStickBreakingHyperLogPrior();
        // total += BaseStickBreakingLogPrior();
        total += BaseLogPrior();
        return total;
    }

    double LambdaHyperLogPrior() const { return -lambda / 10; }

    double GlobalBranchLengthsLogPrior() const {
        return LambdaHyperLogPrior() + branchlength->GetLogProb();
    }

    // exponential of mean 1 for blhyperinvshape
    double BranchLengthsHyperInvShapeLogPrior() const { return -blhyperinvshape; }

    double GeneBranchLengthsHyperLogPrior() const {
        return LambdaHyperLogPrior() + BranchLengthsHyperInvShapeLogPrior() +
               branchlength->GetLogProb();
    }

    double GeneNucRatesHyperLogPrior() const {
        double total = 0;
        if (nucmode == 1) {
            total -= nucrelratehyperinvconc;
            total -= nucstathyperinvconc;
        }
        return total;
    }

    double OmegaHyperLogPrior() const {
        double total = 0;
        if (omegaprior == 0) {
            total -= omegahypermean;
            total -= omegahyperinvshape;
        } else if (omegaprior == 1) {
            double pialpha = dposompihypermean / dposompihyperinvconc;
            double pibeta = (1 - dposompihypermean) / dposompihyperinvconc;
            total += (pialpha - 1) * log(1.0 - dposompi) + (pibeta - 1) * log(dposompi);
            total -= dposomhypermean;
            total -= dposomhyperinvshape;
            if (modalprior && (dposomhyperinvshape > 1.0))  {
                total += Random::INFPROB;
            }
        }
        return total;
    }

    double OmegaLogPrior() const {
        double ret = 0;
        if (omegaprior == 0) {
            ret += omegaarray->GetLogProb();
        } else {
            ret += dposomarray->GetLogProb();
        }
        return ret;
    }

    double BaseStickBreakingHyperLogPrior() const { return -basekappa / 10; }

    double BaseStickBreakingSuffStatLogProb() const {
        double ret = baseweight->GetMarginalLogProb(*baseoccupancy);
        if (std::isinf(ret)) {
            cerr << "in base stick breaking suff stat log prob: inf\n";
            exit(1);
        }
        return ret;
    }

    double BaseLogPrior() const {
        double total = 0;
        total += basecenterarray->GetLogProb();
        total += baseconcentrationarray->GetLogProb();
        if (std::isinf(total)) {
            cerr << "in BaseLogPrior: inf\n";
            exit(1);
        }
        return total;
    }

    double BaseLogPrior(int k) const {
        double total = 0;
        total += basecenterarray->GetLogProb(k);
        total += baseconcentrationarray->GetLogProb(k);
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

    // suff stats for moving omega hyper parameters
    double OmegaHyperSuffStatLogProb() const {
        double ret = 0;
        if (omegaprior == 0) {
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            ret = omegahypersuffstat.GetLogProb(alpha, beta);
        } else {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            ret = dposomhypersuffstat.GetLogProb(dposompi, alpha, beta);
        }
        return ret;
    }

    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(basecenterarray->GetVal(k),
                                                       baseconcentrationarray->GetVal(k));
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

    // log prob for moving omega hyperparameters
    double OmegaHyperLogProb() const { return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb(); }

    // for moving aa hyper params (aacenter and aainvconc)
    // for component k of the mixture
    double BaseLogProb(int k) const { return BaseLogPrior(k) + BaseSuffStatLogProb(k); }

    // for moving kappa
    double BaseStickBreakingHyperLogProb() const {
        return BaseStickBreakingHyperLogPrior() + BaseStickBreakingSuffStatLogProb();
        // return BaseStickBreakingHyperLogPrior() + BaseStickBreakingLogPrior();
    }

    //-------------------
    // Moves
    //-------------------

    void MasterMove() override {
        totchrono.Start();
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            paramchrono.Start();
            aachrono.Start();

            aachrono.Stop();

            if (basemode >= 2) {
                basechrono.Start();
                for (int r = 0; r < 5; r++) {
                    MasterReceiveBaseSuffStat();
                    movechrono.Start();
                    MoveBaseMixture(1);
                    movechrono.Stop();
                    MasterSendBaseMixture();
                }
                basechrono.Stop();
            }

            if ((burnin > 10) && (omegamode != 3)) {
                MasterReceiveOmega();
                movechrono.Start();
                MoveOmegaHyperParameters();
                movechrono.Stop();
                MasterSendOmegaHyperParameters();
            }

            blchrono.Start();

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                MasterReceiveBranchLengthsSuffStat();
                movechrono.Start();
                ResampleBranchLengths();
                MoveLambda();
                movechrono.Stop();
                MasterSendGlobalBranchLengths();
            } else if (blmode == 1) {
                MasterReceiveBranchLengthsHyperSuffStat();
                movechrono.Start();
                MoveBranchLengthsHyperParameters();
                movechrono.Stop();
                MasterSendBranchLengthsHyperParameters();
            }
            blchrono.Stop();

            if (nucmode == 1) {
                MasterReceiveNucRatesHyperSuffStat();
                movechrono.Start();
                MoveNucRatesHyperParameters();
                movechrono.Stop();
                MasterSendNucRatesHyperParameters();
            }

            paramchrono.Stop();
        }

        if (blmode != 2) {
            MasterReceiveGeneBranchLengths();
        }
        MasterReceiveGeneNucRates();
        MasterReceiveOmega();
        MasterReceiveLogProbs();
        MasterReceivePredictedDNDS();
        
        totchrono.Stop();

        burnin++;
    }

    // slave move
    void SlaveMove() override {
        movechrono.Start();
        mapchrono.Start();
        GeneResampleSub(1.0);
        mapchrono.Stop();
        movechrono.Stop();

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            movechrono.Start();
            GeneCollectPathSuffStat();
            MoveGeneAA();
            movechrono.Stop();

            if (basemode >= 2) {
                for (int r = 0; r < 5; r++) {
                    movechrono.Start();
                    MoveGeneBase();
                    movechrono.Stop();
                    SlaveSendBaseSuffStat();
                    SlaveReceiveBaseMixture();
                }
            } else {
                movechrono.Start();
                for (int r = 0; r < 5; r++) {
                    MoveGeneBase();
                }
                movechrono.Stop();
            }

            if ((burnin > 10) && (omegamode != 3)) {
                movechrono.Start();
                MoveGeneOmegas();
                movechrono.Stop();
                SlaveSendOmega();
                SlaveReceiveOmegaHyperParameters();
            }

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                SlaveSendBranchLengthsSuffStat();
                SlaveReceiveGlobalBranchLengths();
            } else {
                MoveGeneBranchLengths();
                if (blmode == 1) {
                    SlaveSendBranchLengthsHyperSuffStat();
                    SlaveReceiveBranchLengthsHyperParameters();
                }
            }

            // global nucrates, or gene nucrates hyperparameters
            MoveGeneNucRates();
            if (nucmode == 1) {
                SlaveSendNucRatesHyperSuffStat();
                SlaveReceiveNucRatesHyperParameters();
            }

            movechrono.Start();
            movechrono.Stop();
        }

        if (blmode != 2) {
            SlaveSendGeneBranchLengths();
        }
        SlaveSendGeneNucRates();
        SlaveSendOmega();
        SlaveSendLogProbs();
        SlaveSendPredictedDNDS();

        burnin++;
    }

    void GeneResampleSub(double frac) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void GeneCollectPathSuffStat() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectSitePathSuffStat();
            geneprocess[gene]->CollectComponentPathSuffStat();
        }
    }

    void MoveGeneOmegas() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveOmega();
            if (omegaprior == 0) {
                (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
            } else {
                (*dposomarray)[gene] = geneprocess[gene]->GetOmega() - 1;
            }
        }
    }

    void MoveGeneAA() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveAAMixture(3);
        }
    }

    void MoveGeneBase() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveBase(3);
            geneprocess[gene]->CollectBaseSuffStat();
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

    void MoveBaseMixture(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseComponents(10);
            if (baseNcat > 1) {
                ResampleBaseEmptyComponents();
                MoveBaseKappa();
                if (!basemin) {
                    BaseLabelSwitchingMove();
                }
            }
        }
    }

    void MoveBaseComponents(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            MoveBaseCenters(1.0, 1);
            MoveBaseCenters(1.0, 3);
            MoveBaseCenters(0.3, 3);
            MoveBaseConcentrations(1.0);
            MoveBaseConcentrations(0.3);
        }
    }

    double MoveBaseCenters(double tuning, int n) {
        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa, 0);
        for (int k = 0; k < baseNcat; k++) {
            if (baseoccupancy->GetVal(k)) {
                vector<double> &aa = (*basecenterarray)[k];
                bk = aa;
                double deltalogprob = -BaseLogProb(k);
                double loghastings = Random::ProfileProposeMove(aa, Naa, tuning, n);
                deltalogprob += loghastings;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    aa = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveBaseConcentrations(double tuning) {
        double nacc = 0;
        double ntot = 0;
        for (int k = basemin; k < baseNcat; k++) {
            if (baseoccupancy->GetVal(k)) {
                double &c = (*baseconcentrationarray)[k];
                double bk = c;
                double deltalogprob = -BaseLogProb(k);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                c *= e;
                deltalogprob += m;
                deltalogprob += BaseLogProb(k);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    c = bk;
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    void ResampleBaseEmptyComponents() {
        basecenterarray->PriorResample(*baseoccupancy);
        baseconcentrationarray->PriorResample(*baseoccupancy);
        if (basemin == 1) {
            (*baseconcentrationarray)[0] = 1.0;
        }
    }

    void BaseLabelSwitchingMove() {
        permutocc->Reset();
        baseweight->LabelSwitchingMove(5, *baseoccupancy, *permutocc);
        basecenterarray->Permute(*permutocc);
        baseconcentrationarray->Permute(*permutocc);
        basesuffstatarray->Permute(*permutocc);
    }

    void ResampleBaseWeights() { baseweight->GibbsResample(*baseoccupancy); }

    void MoveBaseKappa() {
        ScalingMove(basekappa, 1.0, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::BaseStickBreakingHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(basekappa, 0.3, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::BaseStickBreakingHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        baseweight->SetKappa(basekappa);
        baseweight->GibbsResample(*baseoccupancy);
    }

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneAAMutSelDSBDPOmegaModel::LambdaHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneAAMutSelDSBDPOmegaModel::LambdaHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, 10);
            BranchLengthsHyperScalingMove(0.3, 10);
        }

        ScalingMove(blhyperinvshape, 1.0, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);

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
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveOmegaHyperParameters() {
        if (omegaprior == 0) {
            omegahypersuffstat.Clear();
            omegahypersuffstat.AddSuffStat(*omegaarray);

            ScalingMove(omegahypermean, 1.0, 10,
                        &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                        &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
            ScalingMove(omegahypermean, 0.3, 10,
                        &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                        &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
            ScalingMove(omegahyperinvshape, 1.0, 10,
                        &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                        &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
            ScalingMove(omegahyperinvshape, 0.3, 10,
                        &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                        &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);

            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            omegaarray->SetShape(alpha);
            omegaarray->SetScale(beta);
        } else {
            dposomhypersuffstat.Clear();
            dposomhypersuffstat.AddSuffStat(*dposomarray);

            ScalingMove(dposomhypermean, 1.0, 10,
                        &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                        &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
            ScalingMove(dposomhypermean, 0.3, 10,
                        &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                        &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);

            if (modalprior) {
                SlidingMove(dposomhyperinvshape,1.0,10,0,1.0,&MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate,this);
                SlidingMove(dposomhyperinvshape,0.3,10,0,1.0,&MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate,this);
                SlidingMove(dposomhyperinvshape,0.1,10,0,1.0,&MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,&MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate,this);
            } else  {
                ScalingMove(dposomhyperinvshape, 1.0, 10,
                            &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                            &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
                ScalingMove(dposomhyperinvshape, 0.3, 10,
                            &MultiGeneAAMutSelDSBDPOmegaModel::OmegaHyperLogProb,
                            &MultiGeneAAMutSelDSBDPOmegaModel::NoUpdate, this);
            }

            if (burnin > 10) {
                if (dposompihyperinvconc) {
                    ResampleDPosOmPi();
                }
            }
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            dposomarray->SetPi(dposompi);
            dposomarray->SetShape(alpha);
            dposomarray->SetScale(beta);
        }
    }

    void ResampleDPosOmPi() {
        int n0 = dposomhypersuffstat.GetN0();
        int n1 = dposomhypersuffstat.GetN1();
        if ((n0 + n1) != Ngene) {
            cerr << "error in resample pi\n";
            exit(1);
        }
        double pialpha = dposompihypermean / dposompihyperinvconc;
        double pibeta = (1 - dposompihypermean) / dposompihyperinvconc;
        double a0 = Random::sGamma(pialpha + n1);
        double a1 = Random::sGamma(pibeta + n0);
        dposompi = a1 / (a0 + a1);
    }

    //-------------------
    // MPI send / receive
    //-------------------

    // global branch lengths

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

    // omega (and hyperparameters)

    void SlaveSendOmega() {
        if (omegaprior == 0) {
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
            }
            SlaveSendGeneArray(*omegaarray);
        } else {
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                (*dposomarray)[gene] = geneprocess[gene]->GetOmega() - 1.0;
            }
            SlaveSendGeneArray(*dposomarray);
        }
    }

    void MasterReceiveOmega() {
        if (omegaprior == 0) {
            MasterReceiveGeneArray(*omegaarray);
        } else {
            MasterReceiveGeneArray(*dposomarray);
        }
    }

    void MasterSendOmega() {
        if (omegaprior == 0) {
            MasterSendGeneArray(*omegaarray);
        } else {
            MasterSendGeneArray(*dposomarray);
        }
    }

    void SlaveReceiveOmega() {
        if (omegaprior == 0) {
            SlaveReceiveGeneArray(*omegaarray);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmega((*omegaarray)[gene]);
            }
        } else {
            SlaveReceiveGeneArray(*dposomarray);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmega((*dposomarray)[gene] + 1);
            }
        }
    }

    // omega hyperparameters

    void MasterSendOmegaHyperParameters() {
        if (omegaprior == 0) {
            MasterSendGlobal(omegahypermean, omegahyperinvshape);
        } else {
            MasterSendGlobal(dposompi);
            MasterSendGlobal(dposomhypermean, dposomhyperinvshape);
        }
    }

    void SlaveReceiveOmegaHyperParameters() {
        if (omegaprior == 0) {
            SlaveReceiveGlobal(omegahypermean, omegahyperinvshape);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape);
            }
        } else {
            SlaveReceiveGlobal(dposompi);
            SlaveReceiveGlobal(dposomhypermean, dposomhyperinvshape);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetDPosOmHyperParameters(dposompi, dposomhypermean,
                                                            dposomhyperinvshape);
            }
        }
    }

    // aa base mixture

    void MasterSendBaseMixture() {
        MasterSendGlobal(*basecenterarray, *baseconcentrationarray);
        MasterSendGlobal(*baseweight, *permutocc);
    }

    void SlaveReceiveBaseMixture() {
        SlaveReceiveGlobal(*basecenterarray, *baseconcentrationarray);
        SlaveReceiveGlobal(*baseweight, *permutocc);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBaseMixture(*basecenterarray, *baseconcentrationarray,
                                              *baseweight, *permutocc);
        }
    }

    // aa hyper suff stat

    void SlaveSendBaseSuffStat() {
        basesuffstatarray->Clear();
        baseoccupancy->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectBaseSuffStat();
            basesuffstatarray->Add(*geneprocess[gene]->GetBaseSuffStatArray());
            geneprocess[gene]->UpdateBaseOccupancies();
            baseoccupancy->Add(*geneprocess[gene]->GetBaseOccupancies());
        }
        SlaveSendAdditive(*basesuffstatarray);
        SlaveSendAdditive(*baseoccupancy);
    }

    void MasterReceiveBaseSuffStat() {
        basesuffstatarray->Clear();
        baseoccupancy->Clear();
        MasterReceiveAdditive(*basesuffstatarray);
        MasterReceiveAdditive(*baseoccupancy);
    }

    // branch length suff stat
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

    // log probs

    void SlaveSendLogProbs() {
        GeneLogPrior = 0;
        lnL = 0;
        MeanNcluster = 0;
        MeanBaseNcluster = 0;
        MeanStatEnt = 0;
        MeanAAConc = 0;
        MeanAACenterEnt = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
            MeanNcluster += geneprocess[gene]->GetNcluster();
            MeanBaseNcluster += geneprocess[gene]->GetBaseNcluster();
            MeanStatEnt += geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanAAEntropy();
            MeanAAConc += geneprocess[gene]->GetNsite() *
                          geneprocess[gene]->GetMeanComponentAAConcentration();
            MeanAACenterEnt +=
                geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanComponentAAEntropy();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
        SlaveSendAdditive(MeanNcluster);
        SlaveSendAdditive(MeanBaseNcluster);
        SlaveSendAdditive(MeanStatEnt);
        SlaveSendAdditive(MeanAAConc);
        SlaveSendAdditive(MeanAACenterEnt);

        moveTime = movechrono.GetTime();
        mapTime = mapchrono.GetTime();
        SlaveSendAdditive(moveTime);
        SlaveSendAdditive(mapTime);
    }

    void MasterReceiveLogProbs() {
        GeneLogPrior = 0;
        lnL = 0;
        MasterReceiveAdditive(GeneLogPrior);
        MasterReceiveAdditive(lnL);
        MeanNcluster = 0;
        MeanBaseNcluster = 0;
        MeanStatEnt = 0;
        MeanAAConc = 0;
        MeanAACenterEnt = 0;
        MasterReceiveAdditive(MeanNcluster);
        MasterReceiveAdditive(MeanBaseNcluster);
        MasterReceiveAdditive(MeanStatEnt);
        MasterReceiveAdditive(MeanAAConc);
        MasterReceiveAdditive(MeanAACenterEnt);
        MeanNcluster /= GetLocalNgene();
        MeanBaseNcluster /= GetLocalNgene();
        MeanStatEnt /= GetTotNsite();
        MeanAAConc /= GetTotNsite();
        MeanAACenterEnt /= GetTotNsite();

        moveTime = 0;
        mapTime = 0;
        MasterReceiveAdditive(moveTime);
        MasterReceiveAdditive(mapTime);
        moveTime /= (GetNprocs() - 1);
        mapTime /= (GetNprocs() - 1);
    }
};
