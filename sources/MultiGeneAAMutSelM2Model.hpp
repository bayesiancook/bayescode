
#include "AAMutSelM2Model.hpp"
#include "IIDBernoulliBeta.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "Permutation.hpp"

class MultiGeneAAMutSelM2Model : public MultiGeneProbModel {
  private:
    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;
    std::vector<CodonSequenceAlignment*> alivector;

    string datapath;
    string datafile;
    string treefile;

    int Ntaxa;
    int Nbranch;

    int basemin;
    int baseNcat;
    int Ncat;

    // omega mixture
    vector<double> mixhyperparam;

    double &dposomhypermean;
    double &dposomhyperinvshape;
    IIDGamma *dposomarray;
    GammaSuffStat dposomsuffstat;

    double &poswhypermean;
    double &poswhyperinvconc;
    IIDBernoulliBeta *poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double &pi;

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

    std::vector<AAMutSelM2Model *> geneprocess;

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

    int blmode, nucmode, basemode, poswmode, dposommode, modalprior;

    Chrono totchrono;
    Chrono paramchrono;
    Chrono basechrono;
    Chrono blchrono;
    Chrono aachrono;

    int burnin;
    int chainsize;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneAAMutSelM2Model(string indatafile, string intreefile, int inNcat, int inbaseNcat,
                                     int inblmode, int innucmode, int inbasemode,
                                     int inposwmode, int indposommode,
                                     int inmodalprior, double inpihypermean, double inpihyperinvconc,
                                     int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs), 
        
          mixhyperparam(5, 0),
          dposomhypermean(mixhyperparam[0]),
          dposomhyperinvshape(mixhyperparam[1]),
          poswhypermean(mixhyperparam[2]),
          poswhyperinvconc(mixhyperparam[3]),
          pi(mixhyperparam[4]),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {

        datafile = indatafile;
        treefile = intreefile;
        AllocateAlignments(datafile);

        Ncat = inNcat;

        burnin = 20;
        chainsize = 0;

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
        poswmode = inposwmode;
        dposommode = indposommode;

        modalprior = inmodalprior;

        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;
        pi = pihypermean;

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

    void SetChainSize(double insize)	{
	    chainsize = insize;
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

        // omega mixture
        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray = new IIDGamma(GetLocalNgene(), dposomalpha, dposombeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        poswarray = new IIDBernoulliBeta(GetLocalNgene(), pi, poswalpha, poswbeta);

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
            geneprocess.assign(0, (AAMutSelM2Model *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (AAMutSelM2Model *)0);

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
                    geneprocess[gene] = new AAMutSelM2Model(alivector[gene], tree, Ncat, baseNcat);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    geneprocess[gene] = new AAMutSelM2Model(GetLocalGeneName(gene), treefile, Ncat, baseNcat);
                }
            }

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->SetBaseMode(basemode);
                geneprocess[gene]->SetOmegaMixtureHyperParameters(pi,poswhypermean, poswhyperinvconc,
                        dposomhypermean, dposomhyperinvshape);
                geneprocess[gene]->Allocate();
            }
        }
    }

    void SetOmegaMixtureHyperParameters(double inposwhypermean, double inposwhyperinvconc,
            double indposomhypermean, double indposomhyperinvshape) {
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
    }

    void SetOmegaMixtureArrays() {

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray->SetShape(dposomalpha);
        dposomarray->SetScale(dposombeta);
        if (myid) {
            dposomarray->PriorResample(*poswarray);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmegaMixtureParameters((*poswarray)[gene], (*dposomarray)[gene]);
            }
        }

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        poswarray->SetPi(pi);
        poswarray->SetAlpha(poswalpha);
        poswarray->SetBeta(poswbeta);
    }

    void FastUpdate() {
        branchlength->SetScale(lambda);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);

        SetOmegaMixtureArrays();
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

            MasterSendOmegaMixtureHyperParameters();
            MasterSendOmegaMixtureParameters();

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

        SlaveReceiveOmegaMixtureHyperParameters();
        SlaveReceiveOmegaMixtureParameters();

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

            MasterSendOmegaMixtureHyperParameters();
            MasterSendOmegaMixtureParameters();

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

        SlaveReceiveOmegaMixtureHyperParameters();
        SlaveReceiveOmegaMixtureParameters();

        GenePostPred(name);
        // SlaveSendLogProbs();
    }

    void GenePostPred(string name) {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->PostPred(name + GetLocalGeneName(gene));
        }
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
        os << "pi\t";
        os << "nposfrac\t";
        os << "dposommean\tinvshape\t";
        os << "poswmean\tinvconc\t";
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

        // omega
        os << pi << '\t';
        os << GetNpos() << '\t';
        os << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
        os << poswhypermean << '\t' << poswhyperinvconc << '\t';

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

    void TracePosWeight(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << poswarray->GetVal(gene) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TracePosOm(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            if (poswarray->GetVal(gene))    {
                os << 1 + dposomarray->GetVal(gene) << '\t';
            }
            else    {
                os << 1.0 << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

    int GetNpos() const { return GetNgene() - poswarray->GetNullSet(); }

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

        // omega
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> poswhypermean >> poswhyperinvconc;
        is >> pi;
        is >> *dposomarray;
        is >> *poswarray;

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

        // omega
        os << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
        os << poswhypermean << '\t' << poswhyperinvconc << '\t';
        os << pi << '\t';
        os << *dposomarray << '\t';
        os << *poswarray << '\t';

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
            os << size << '\t';
            buffer.ToStream(os);
        }

        os << '\n';
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

        total += OmegaMixtureHyperLogPrior();

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


    double OmegaMixtureHyperLogPrior() const {
        double total = 0;
        if (pi) {
            // beta distribution for pi, if not 0
            double pialpha = pihypermean / pihyperinvconc;
            double pibeta = (1 - pihypermean) / pihyperinvconc;
            total += (pialpha - 1) * log(1.0 - pi) + (pibeta - 1) * log(pi);
        }

        // uniform prior for poswhypermean

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
        if (modalprior == 2)    {
            if (dposomhypermean < 0.5)  {
                total += log(0);
                // total += Random::INFPROB;
            }
        }

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

    // suff stat for gene-specific mixture parameters, as a function of mixture
    // hyperparameters
    double OmegaMixtureHyperSuffStatLogProb() const {
        double total = 0;

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        total += dposomsuffstat.GetLogProb(dposomalpha, dposombeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        total += poswsuffstat.GetLogProb(pi, poswalpha, poswbeta);
        return total;
    }

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

    double BaseSuffStatLogProb(int k) const {
        return basesuffstatarray->GetVal(k).GetLogProb(basecenterarray->GetVal(k),
                                                       baseconcentrationarray->GetVal(k));
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // logprob for omega mixture hyperparameters
    double OmegaMixtureHyperLogProb() const {
        return OmegaMixtureHyperLogPrior() + OmegaMixtureHyperSuffStatLogProb();
    }

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

            // mixture hyperparameters

            MasterReceiveOmegaMixtureHyperSuffStat();
            MoveOmegaMixtureHyperParameters();
            MasterSendOmegaMixtureHyperParameters();

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
        MasterReceiveOmegaMixtureParameters();
        MasterReceiveLogProbs();
        
        totchrono.Stop();

        chainsize++;
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

            MoveGeneOmegaMixtureParameters();
            SlaveSendOmegaMixtureHyperSuffStat();
            SlaveReceiveOmegaMixtureHyperParameters();
            SetOmegaMixtureArrays();

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
        SlaveSendOmegaMixtureParameters();
        SlaveSendLogProbs();

        chainsize++;
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

    void MoveGeneOmegaMixtureParameters() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveOmega();
            geneprocess[gene]->GetOmegaMixtureParameters((*poswarray)[gene], (*dposomarray)[gene]);
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
                    &MultiGeneAAMutSelM2Model::BaseStickBreakingHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(basekappa, 0.3, 10,
                    &MultiGeneAAMutSelM2Model::BaseStickBreakingHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        baseweight->SetKappa(basekappa);
        baseweight->GibbsResample(*baseoccupancy);
    }

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneAAMutSelM2Model::LambdaHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneAAMutSelM2Model::LambdaHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, 10);
            BranchLengthsHyperScalingMove(0.3, 10);
        }

        ScalingMove(blhyperinvshape, 1.0, 10,
                    &MultiGeneAAMutSelM2Model::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10,
                    &MultiGeneAAMutSelM2Model::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);

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
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelM2Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelM2Model::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveOmegaMixtureHyperParameters() {

        if (dposommode == 1) {
            ScalingMove(dposomhypermean, 1.0, 10, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            ScalingMove(dposomhypermean, 0.3, 10, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            ScalingMove(dposomhyperinvshape, 1.0, 10, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            ScalingMove(dposomhyperinvshape, 0.3, 10, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
        }

        if (poswmode == 1) {
            MovePoswHyper();
        }

        if (chainsize >= burnin)    {
            if (pihyperinvconc) {
                ResamplePi();
            }
        }
    }


    void MovePoswHyper()	{

            SlidingMove(poswhypermean, 1.0, 10, 0, 1, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            SlidingMove(poswhypermean, 0.3, 10, 0, 1, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            SlidingMove(poswhypermean, 0.1, 10, 0, 1, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 1.0, 10, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 0.3, 10, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 0.1, 10, &MultiGeneAAMutSelM2Model::OmegaMixtureHyperLogProb,
                        &MultiGeneAAMutSelM2Model::NoUpdate, this);

        for (int rep=0; rep<10; rep++)	{
            PoswCompMove(1.0);
            PoswCompMove(0.3);
            PoswCompMove(0.1);
        }
    }

    int PoswCompMove(double tuning)	{

        double bkposwhypermean = poswhypermean;
        double bkposwhyperinvconc = poswhyperinvconc;
        double deltalogprob = - OmegaMixtureHyperLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        if (poswhypermean * e > 1.0)	{
            return 0;
        }
        poswhypermean *= e;
        poswhyperinvconc *= e;
        deltalogprob += OmegaMixtureHyperLogProb();
        deltalogprob += 2*m;

        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted) {
            poswhypermean = bkposwhypermean;
            poswhyperinvconc = bkposwhyperinvconc;
        }
        return accepted;
    }


    void ResamplePi() {
        int n0 = poswsuffstat.GetN0();
        int n1 = poswsuffstat.GetN1();
        if ((n0 + n1) != Ngene) {
            cerr << "error in resample pi\n";
            exit(1);
        }
        double pialpha = pihypermean / pihyperinvconc;
        double pibeta = (1 - pihypermean) / pihyperinvconc;
        double postalpha = Random::sGamma(pialpha + n1);
        double postbeta = Random::sGamma(pibeta + n0);
        pi = postalpha / (postalpha + postbeta);
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

    void MasterSendOmegaMixtureHyperParameters() { MasterSendGlobal(mixhyperparam); }

    void SlaveReceiveOmegaMixtureHyperParameters() {
        SlaveReceiveGlobal(mixhyperparam);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmegaMixtureHyperParameters(pi, poswhypermean, poswhyperinvconc,
                    dposomhypermean, dposomhyperinvshape);
        }
    }


    void SlaveSendOmegaMixtureHyperSuffStat() {
        dposomsuffstat.Clear();
        dposomsuffstat.AddSuffStat(*dposomarray, *poswarray);
        SlaveSendAdditive(dposomsuffstat);
        poswsuffstat.Clear();
        poswarray->AddSuffStat(poswsuffstat);
        SlaveSendAdditive(poswsuffstat);
    }

    void MasterReceiveOmegaMixtureHyperSuffStat() {
        dposomsuffstat.Clear();
        MasterReceiveAdditive(dposomsuffstat);
        poswsuffstat.Clear();
        MasterReceiveAdditive(poswsuffstat);
    }

    void SlaveSendOmegaMixtureParameters() {
        SlaveSendGeneArray(*poswarray, *dposomarray);
    }

    void MasterReceiveOmegaMixtureParameters() {
        MasterReceiveGeneArray(*poswarray, *dposomarray);
    }

    void MasterSendOmegaMixtureParameters() {
        MasterSendGeneArray(*poswarray, *dposomarray);
    }

    void SlaveReceiveOmegaMixtureParameters() {
        SlaveReceiveGeneArray(*poswarray, *dposomarray);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetOmegaMixtureParameters((*poswarray)[gene], (*dposomarray)[gene]);
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

    void MasterTraceSitesPostProb(ostream &os) {
        for (int proc = 1; proc < GetNprocs(); proc++) {
            int totnsite = GetSlaveTotNsite(proc);
            double *array = new double[totnsite];
            MPI_Status stat;
            MPI_Recv(array, totnsite, MPI_DOUBLE, proc, TAG1, MPI_COMM_WORLD, &stat);

            int i = 0;
            for (int gene = 0; gene < Ngene; gene++) {
                if (GeneAlloc[gene] == proc) {
                    os << GeneName[gene] << '\t';
                    int nsite = GeneNsite[gene];
                    for (int k = 0; k < nsite; k++) {
                        if (array[i] < 0) {
                            cerr << "error: negative post prob\n";
                            cerr << GeneName[gene] << '\n';
                            cerr << GeneNsite[gene] << '\n';
                            cerr << i << '\n';
                            exit(1);
                        }
                        os << array[i++] << '\t';
                    }
                }
            }
            if (i != totnsite) {
                cerr << "error in MultiGeneCodonM2aModel::MasterTraceSitesPostProb: non "
                        "matching number of sites\n";
                cerr << i << '\t' << totnsite << '\n';
                exit(1);
            }
            delete[] array;
        }
        os << '\n';
        os.flush();
    }

    void SlaveTraceSitesPostProb() {
        int ngene = GetLocalNgene();
        int totnsite = GetLocalTotNsite();
        double *array = new double[totnsite];
        int i = 0;
        for (int gene = 0; gene < ngene; gene++) {
            geneprocess[gene]->GetSitesPostProb(array + i);
            for (int j = 0; j < GeneNsite[gene]; j++) {
                if (array[i + j] < 0) {
                    cerr << "error in slave\n";
                    cerr << i << '\t' << j << '\t' << GeneName[gene] << '\t' << GeneNsite[gene] << '\t'
                         << geneprocess[gene]->GetNsite() << '\n';
                    exit(1);
                }
            }
            i += GetLocalGeneNsite(gene);
        }
        if (i != totnsite) {
            cerr << "error in MultiGeneCodonM2aModel::SlaveTraceSitesPostProb: non "
                    "matching number of sites\n";
            exit(1);
        }

        MPI_Send(array, totnsite, MPI_DOUBLE, 0, TAG1, MPI_COMM_WORLD);
        delete[] array;
    }
};
