
#include "DiscSelACModel.hpp"
#include "IIDGamma.hpp"
#include "IIDNormal.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "Permutation.hpp"

class MultiGeneDiscSelACModel : public MultiGeneProbModel {
  private:
    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;
    std::vector<CodonSequenceAlignment*> alivector;

    string datapath;
    string datafile;
    string treefile;
    string initfile;

    int Ntaxa;
    int Nbranch;

    // branch lengths
    double lambda;
    BranchIIDGamma *branchlength;
    GammaSuffStat hyperlengthsuffstat;

    double blhyperinvshape;
    GammaWhiteNoiseArray *branchlengtharray;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStatBranchArray *lengthhypersuffstatarray;

    // nucleotide rates hyperparameters
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucstat;
    std::vector<double> nucrelrate;
    GTRSubMatrix *nucmatrix;

    SimpleArray<double>* genednds;

    // each gene has its own gamma distribution and its own pi
    int aadistmodel;
    int Gcat;
    double xmin;
    double xmax;
    vector<double> G;

    double Gvarhypermean;
    double Gvarhyperinvshape;
    IIDGamma* Gvararray;
    GammaSuffStat Gvarhypersuffstat;

    double logpsihypermean;
    double logpsihypervar;
    IIDNormal* logpsiarray;
    NormalSuffStat logpsihypersuffstat;

    // sharing the aadist parameters
    double wcom;
    double wpol;
    double wvol;
    // distance matrix between amino acids 
    vector<double> aadist;
    // prior probs for optimal amino acids
    vector<double> aaweights;
    OccupancySuffStat* aaoccupancy;

    SelACProfileBidimArray* selacprofiles;
    AADiffSelCodonMatrixBidimArray* codonmatrices;
    PathSuffStatBidimArray *componentpathsuffstatbidimarray;

    std::vector<DiscSelACModel *> geneprocess;

    double lnL;
    double GeneLogPrior;
    double MeanStatEnt;
    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    // free (0), shrunken (1) or shared (2)
    int blmode;

    // shared (2) or fixed (3)
    int nucmode;
    int aadistmode;
    int gvarmode;
    int psimode;

    Chrono totchrono;
    Chrono paramchrono;
    Chrono blchrono;
    Chrono aachrono;

    vector<double> aadistacc, aadisttot;
    int burnin;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneDiscSelACModel(string indatafile, string intreefile, string ininitfile,
                                     int inGcat, double inxmin, double inxmax, int inaadistmodel,
                                     int inblmode, int innucmode, int inaadistmode, int ingvarmode, int inpsimode,
                                     int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs) {

        datafile = indatafile;
        treefile = intreefile;
        initfile = ininitfile;
        AllocateAlignments(datafile);

        Gcat = inGcat;
        xmin = inxmin;
        xmax = inxmax;

        aadistmodel = inaadistmodel;

        burnin = 0;

        blmode = inblmode;
        nucmode = innucmode;

        aadistmode = inaadistmode;
        gvarmode = ingvarmode;
        psimode = inpsimode;

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

        // nucleotide mutation matrix
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, vector<double>(Nrr, 1.0 / Nrr), ((double)Nrr));
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, vector<double>(Nnuc, 1.0 / Nnuc), ((double)Nnuc));
        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        wcom = grantham_wcom;
        wpol = grantham_wpol;
        wvol = grantham_wvol;

        aadist.assign(Naarr, 1.0);
        if (!aadistmodel)    {
            UpdateGrantham();
        }

        aaweights.assign(Naa, 1.0/Naa);
        aaoccupancy = new OccupancySuffStat(Naa);

        if (! GetMyid())    {
            if (initfile != "None") {
                ifstream is(initfile.c_str());
                for (int a=0; a<Naa; a++)   {
                    is >> aaweights[a];
                }
                for (int i=0; i<Naarr; i++) {
                    is >> aadist[i];
                }
            }
        }

        G.assign(Gcat, 1.0);
        for (int i=0; i<Gcat; i++)  {
            G[i] = exp(xmin + ((double) i)/Gcat*(xmax-xmin));
        }

        Gvarhypermean = 10.0;
        Gvarhyperinvshape = 1.0;
        double alpha = 1.0 / Gvarhyperinvshape;
        double beta = alpha / Gvarhypermean;
        Gvararray = new IIDGamma(GetLocalNgene(), alpha, beta);

        logpsihypermean = 0;
        logpsihypervar = 10.0;
        logpsiarray = new IIDNormal(GetLocalNgene(), logpsihypermean, logpsihypervar);

        selacprofiles = new SelACProfileBidimArray(aadist, G, ((double) Naarr));
        codonmatrices = new AADiffSelCodonMatrixBidimArray(*selacprofiles, *GetCodonStateSpace(), *nucmatrix, 1.0);
        componentpathsuffstatbidimarray = new PathSuffStatBidimArray(Naa, Gcat, GetCodonStateSpace()->GetNstate());

        genednds = new SimpleArray<double>(GetLocalNgene());

        lnL = 0;
        GeneLogPrior = 0;
        MeanStatEnt = 0;

        // aadist params

        if (!GetMyid()) {
            geneprocess.assign(0, (DiscSelACModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (DiscSelACModel *)0);

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
                    geneprocess[gene] = new DiscSelACModel(
                        alivector[gene], tree, aadistmodel, Gcat, xmin, xmax);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    geneprocess[gene] = new DiscSelACModel(GetLocalGeneName(gene), treefile, aadistmodel, Gcat, xmin, xmax);
                }
            }

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->SetAAMode(aadistmode, gvarmode, psimode);
                geneprocess[gene]->Allocate();
            }
        }

        aadistacc.assign(3,0);
        aadisttot.assign(3,0);
    }

    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Naa - i - 1) * i / 2 + j - i - 1
                       : (2 * Naa - j - 1) * j / 2 + i - j - 1;
    }

    void UpdateGrantham()   {
        double tot = 0;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                double tcom = grantham_com[b] - grantham_com[a];
                double tpol = grantham_pol[b] - grantham_pol[a];
                double tvol = grantham_vol[b] - grantham_vol[a];
                double d = sqrt(wcom*tcom*tcom + wpol*tpol*tpol + wvol*tvol*tvol);
                aadist[rrindex(a,b)] = d;
                tot += d;
            }
        }
        tot /= Naarr;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                aadist[rrindex(a,b)] /= tot;
            }
        }
    }

    void UpdateGvar()   {
        double alpha = 1.0 / Gvarhyperinvshape;
        double beta = alpha / Gvarhypermean;
        Gvararray->SetShape(alpha);
        Gvararray->SetScale(beta);
    }

    void UpdateLogPsi() {
        logpsiarray->SetMean(logpsihypermean);
        logpsiarray->SetVar(logpsihypervar);
    }

    void UpdateAADist()  {
        if (! aadistmodel)   {
            UpdateGrantham();
        }
        selacprofiles->Update();
        UpdateCodonMatrices();
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrices that their parameters have changed and that
    //! they should be updated
    void UpdateCodonMatrices() {
        codonmatrices->Corrupt();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that it should be updated
    void UpdateMatrices() {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    void FastUpdate() {
        branchlength->SetScale(lambda);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }

        UpdateGvar();
        UpdateLogPsi();
        UpdateAADist();
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

            MasterSendNucRates();

            MasterSendGvarHyperParameters();
            MasterSendGvar();

            MasterSendLogPsiHyperParameters();
            MasterSendLogPsi();

            MasterSendAADist();

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

        SlaveReceiveNucRates();

        SlaveReceiveGvarHyperParameters();
        SlaveReceiveGvar();

        SlaveReceiveLogPsiHyperParameters();
        SlaveReceiveLogPsi();

        SlaveReceiveAADist();

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

            MasterSendNucRates();

            MasterSendGvarHyperParameters();
            MasterSendGvar();

            MasterSendLogPsiHyperParameters();
            MasterSendLogPsi();

            MasterSendAADist();
        }
    }

    void SlavePostPred(string name) override {
        if (blmode >= 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveBranchLengthsHyperParameters();
            SlaveReceiveGeneBranchLengths();
        }

        SlaveReceiveNucRates();

        SlaveReceiveGvarHyperParameters();
        SlaveReceiveGvar();

        SlaveReceiveLogPsiHyperParameters();
        SlaveReceiveLogPsi();

        SlaveReceiveAADist();

        GenePostPred(name);
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

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\t";
        os << "length\t";
        if (blmode < 2) {
            os << "stdev\t";
        }

        os << "meangvar\t";
        os << "var\t";
        os << "meanlogpsi\t";
        os << "var\t";

        if (! aadistmodel)   {
            os << "w_comp\t";
            os << "w_pol\t";
        }
        os << "distmean\t";
        os << "distvar\t";
        os << "weightent\t";

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
            os << 3 * GetMeanTotalLength() << '\t';
        } else {
            os << 3 * GetMeanLength() << '\t';
            os << 3 * sqrt(GetVarLength()) << '\t';
        }

        os << Gvararray->GetMean() << '\t';
        os << Gvararray->GetVar() << '\t';
        os << logpsiarray->GetMean() << '\t';
        os << logpsiarray->GetVar() << '\t';

        if (! aadistmodel)   {
            os << wcom << '\t';
            os << wpol << '\t';
        }
        os << GetMeanAADist() << '\t';
        os << GetVarAADist() << '\t';
        os << Random::GetEntropy(aaweights) << '\t';

        os << Random::GetEntropy(nucrelratehypercenter) << '\t' << nucrelratehyperinvconc << '\t';
        os << Random::GetEntropy(nucstathypercenter) << '\t' << nucstathyperinvconc << '\n';
        os.flush();
    }

    double GetMeanAADist() const    {
        double m1 = 0;
        for (int i=0; i<Naarr; i++) {
            m1 += aadist[i];
        }
        return m1 / Naarr;
    }

    double GetVarAADist() const {
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<Naarr; i++) {
            m1 += aadist[i];
            m2 += aadist[i]*aadist[i];
        }
        m1 /= Naarr;
        m2 /= Naarr;
        m2 -= m1*m1;
        return m2;
    }


    void TracePsi(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << logpsiarray->GetVal(gene) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TraceAADistHeader(ostream& os) const {
        os << "#";
        for (int a=0; a<Naa; a++)   {
            os << AminoAcids[a] << '\t';
        }
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                os << AminoAcids[a] << AminoAcids[b] << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

    void TraceAADist(ostream& os) const {
        for (int a=0; a<Naa; a++)   {
            os << aaweights[a] << '\t';
        }
        double m = GetMeanAADist();
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                os << aadist[rrindex(a,b)] / m << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

    void Monitor(ostream &os) const override {
        os << totchrono.GetTime() << '\t' << paramchrono.GetTime() 
           << '\t' << blchrono.GetTime() << '\t' << aachrono.GetTime() << '\n';
        os << "prop time in param moves: " << paramchrono.GetTime() / totchrono.GetTime() << '\n';
        os << "sub prop time in bl moves   : " << blchrono.GetTime() / paramchrono.GetTime()
           << '\n';
        os << "sub prop time in aa moves   : " << aachrono.GetTime() / paramchrono.GetTime()
           << '\n';

        os << "aa : ";
       for (size_t i=0; i<aadistacc.size(); i++)    {
          os << aadistacc[i] / aadisttot[i] << '\t';
       }
       os << '\n';
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
        is >> nucrelrate;
        is >> nucstat;

        is >> Gvarhypermean >> Gvarhyperinvshape;
        is >> *Gvararray;
        is >> logpsihypermean >> logpsihypervar;
        is >> *logpsiarray;

        if (! aadistmodel)   {
            is >> wcom >> wpol;
        }
        else    {
            is >> aadist;
        }
        is >> aaweights;

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
        os << nucrelrate<< '\t';
        os << nucstat<< '\t';

        os << Gvarhypermean << '\t' << Gvarhyperinvshape << '\t';
        os << *Gvararray << '\t';

        os << logpsihypermean << '\t' << logpsihypervar << '\t';
        os << *logpsiarray << '\t';

        if (! aadistmodel)   {
            os << wcom << '\t' << wpol << '\t';
        }
        else    {
            os << aadist << '\t';
        }
        os << aaweights << '\t';

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

        total += NucRatesLogPrior();

        total += GvarHyperLogPrior();
        total += GvarLogPrior();

        total += LogPsiHyperLogPrior();
        total += LogPsiLogPrior();

        total += AADistLogPrior();

        return total;
    }

    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(nucrelrate, nucrelratehypercenter,
                                             1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    double GvarHyperLogPrior() const {
        double total = 0;
        total -= 0.1 * Gvarhypermean;
        total -= Gvarhyperinvshape;
        return total;
    }

    double GvarLogPrior() const  {
        return Gvararray->GetLogProb();
    }

    double LogPsiHyperLogPrior() const {
        double total = 0;
        total -= 0.1 * logpsihypermean;
        total -= logpsihypervar;
        return total;
    }

    double LogPsiLogPrior() const  {
        return logpsiarray->GetLogProb();
    }

    double AADistLogPrior() const   {
        double total = 0;
        if (! aadistmodel)   {
            total -= log(wcom) + log(wpol);
        }
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

    double GvarHyperSuffStatLogProb() const  {
        double alpha = 1.0 / Gvarhyperinvshape;
        double beta = alpha / Gvarhypermean;
        return Gvarhypersuffstat.GetLogProb(alpha, beta);
    }

    double LogPsiHyperSuffStatLogProb() const  {
        return logpsihypersuffstat.GetLogProb(logpsihypermean, logpsihypervar);
    }

    //! return log prob of the current substitution mapping, as a function of the
    //! current codon substitution process
    double PathSuffStatLogProb() const {
        return componentpathsuffstatbidimarray->GetLogProb(*codonmatrices);
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

    // log prob for moving nuc rates
    double NucRatesLogProb() const { return NucRatesLogPrior() + PathSuffStatLogProb(); }

    // log prob for moving Gvar hyperparameters
    double GvarHyperLogProb() const { return GvarHyperLogPrior() + GvarHyperSuffStatLogProb(); }

    // log prob for moving logpsi hyperparameters
    double LogPsiHyperLogProb() const { return LogPsiHyperLogPrior() + LogPsiHyperSuffStatLogProb(); }

    // log prob for moving aadist
    double AADistLogProb() const { return AADistLogPrior() + PathSuffStatLogProb(); }

    //-------------------
    // Moves
    //-------------------

    void MasterMove() override {
        totchrono.Start();
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            paramchrono.Start();

            aachrono.Start();

            MasterReceiveGvar();
            MoveGvarHyperParameters();
            MasterSendGvarHyperParameters();

            MasterReceiveLogPsi();
            MoveLogPsiHyperParameters();
            MasterSendLogPsiHyperParameters();

            MasterReceivePathSuffStat();
            MoveAADist();
            MasterSendAADist();

            MasterReceiveAAWeightsSuffStat();
            ResampleAAWeights();
            MasterSendAAWeights();

            aachrono.Stop();

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

            MoveNucRates();
            MasterSendNucRates();

            paramchrono.Stop();
        }

        if (blmode != 2) {
            MasterReceiveGeneBranchLengths();
        }
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
            MoveGeneSelAC();
            // gene aadist
            movechrono.Stop();

            SlaveSendGvar();
            SlaveReceiveGvarHyperParameters();

            SlaveSendLogPsi();
            SlaveReceiveLogPsiHyperParameters();

            SlaveSendPathSuffStat();
            SlaveReceiveAADist();

            SlaveSendAAWeightsSuffStat();
            SlaveReceiveAAWeights();

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

            SlaveReceiveNucRates();
        }

        if (blmode != 2) {
            SlaveSendGeneBranchLengths();
        }
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

    void MoveGeneSelAC() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveSelAC();
            (*logpsiarray)[gene] = geneprocess[gene]->GetLogPsi();
            (*Gvararray)[gene] = geneprocess[gene]->GetGvar();
        }
    }

    void ResampleAAWeights() {
        double total = 0;
        for (int i=0; i<Naa; i++)   {
            aaweights[i] = Random::sGamma(1 + aaoccupancy->GetVal(i));
            total += aaweights[i];
        }
        for (int i=0; i<Naa; i++)   {
            aaweights[i] /= total;
        }
    }

    void SlaveSendAAWeightsSuffStat()    {
        aaoccupancy->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->UpdateAAOccupancies();
            aaoccupancy->Add(*geneprocess[gene]->GetAAOccupancies());
        }
        SlaveSendAdditive(*aaoccupancy);
    }

    void MasterReceiveAAWeightsSuffStat()    {
        aaoccupancy->Clear();
        MasterReceiveAdditive(*aaoccupancy);
    }

    void MoveGvarHyperParameters()   {

        Gvarhypersuffstat.Clear();
        Gvarhypersuffstat.AddSuffStat(*Gvararray);

        ScalingMove(Gvarhypermean, 1.0, 10,
                &MultiGeneDiscSelACModel::GvarHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);
        ScalingMove(Gvarhypermean, 0.3, 10,
                &MultiGeneDiscSelACModel::GvarHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);
        ScalingMove(Gvarhyperinvshape, 1.0, 10,
                &MultiGeneDiscSelACModel::GvarHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);
        ScalingMove(Gvarhyperinvshape, 0.3, 10,
                &MultiGeneDiscSelACModel::GvarHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);

        double alpha = 1.0 / Gvarhyperinvshape;
        double beta = alpha / Gvarhypermean;
        Gvararray->SetShape(alpha);
        Gvararray->SetScale(beta);
    }

    void MoveLogPsiHyperParameters()   {

        logpsihypersuffstat.Clear();
        logpsihypersuffstat.AddSuffStat(*logpsiarray);

        SlidingMove(logpsihypermean, 1.0, 10, 0, 0,
                &MultiGeneDiscSelACModel::LogPsiHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);
        SlidingMove(logpsihypermean, 0.3, 10, 0, 0,
                &MultiGeneDiscSelACModel::LogPsiHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);
        ScalingMove(logpsihypervar, 1.0, 10,
                &MultiGeneDiscSelACModel::LogPsiHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);
        ScalingMove(logpsihypervar, 0.3, 10,
                &MultiGeneDiscSelACModel::LogPsiHyperLogProb, &MultiGeneDiscSelACModel::NoUpdate, this);

        logpsiarray->SetMean(logpsihypermean);
        logpsiarray->SetVar(logpsihypervar);
    }

    double MoveAADist() {
        aadistacc[0] += ProfileMove(aadist, 1.00, 1, 10, &MultiGeneDiscSelACModel::AADistLogProb, &MultiGeneDiscSelACModel::UpdateAADist, this);
        aadisttot[0] ++;
        aadistacc[1] += ProfileMove(aadist, 0.30, 3, 10, &MultiGeneDiscSelACModel::AADistLogProb, &MultiGeneDiscSelACModel::UpdateAADist, this);
        aadisttot[1] ++;
        aadistacc[2] += ProfileMove(aadist, 0.10, 3, 10, &MultiGeneDiscSelACModel::AADistLogProb, &MultiGeneDiscSelACModel::UpdateAADist, this);
        aadisttot[2] ++;
        return 1.0;
    }

    double MoveGranthamWeights()    {
        aadistacc[0] += ScalingMove(wcom, 1.0, 10, &MultiGeneDiscSelACModel::AADistLogProb, &MultiGeneDiscSelACModel::UpdateAADist, this);
        aadisttot[0] ++;
        aadistacc[1] += ScalingMove(wcom, 0.1, 10, &MultiGeneDiscSelACModel::AADistLogProb, &MultiGeneDiscSelACModel::UpdateAADist, this);
        aadisttot[1] ++;
        aadistacc[0] += ScalingMove(wpol, 1.0, 10, &MultiGeneDiscSelACModel::AADistLogProb, &MultiGeneDiscSelACModel::UpdateAADist, this);
        aadisttot[0] ++;
        aadistacc[1] += ScalingMove(wpol, 0.1, 10, &MultiGeneDiscSelACModel::AADistLogProb, &MultiGeneDiscSelACModel::UpdateAADist, this);
        aadisttot[1] ++;
        return 1.0;
    }

    void MoveGeneBranchLengths() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveBranchLengths();
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
    }

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneDiscSelACModel::LambdaHyperLogProb,
                    &MultiGeneDiscSelACModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneDiscSelACModel::LambdaHyperLogProb,
                    &MultiGeneDiscSelACModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, 10);
            BranchLengthsHyperScalingMove(0.3, 10);
        }

        ScalingMove(blhyperinvshape, 1.0, 10,
                    &MultiGeneDiscSelACModel::BranchLengthsHyperLogProb,
                    &MultiGeneDiscSelACModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10,
                    &MultiGeneDiscSelACModel::BranchLengthsHyperLogProb,
                    &MultiGeneDiscSelACModel::NoUpdate, this);

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

    //! MH move on nucleotide rate parameters
    void MoveNucRates() {
        ProfileMove(nucrelrate, 0.1, 1, 3, &MultiGeneDiscSelACModel::NucRatesLogProb,
                    &MultiGeneDiscSelACModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.03, 3, 3, &MultiGeneDiscSelACModel::NucRatesLogProb,
                    &MultiGeneDiscSelACModel::UpdateMatrices, this);
        ProfileMove(nucrelrate, 0.01, 3, 3, &MultiGeneDiscSelACModel::NucRatesLogProb,
                    &MultiGeneDiscSelACModel::UpdateMatrices, this);

        ProfileMove(nucstat, 0.1, 1, 3, &MultiGeneDiscSelACModel::NucRatesLogProb,
                    &MultiGeneDiscSelACModel::UpdateMatrices, this);
        ProfileMove(nucstat, 0.01, 1, 3, &MultiGeneDiscSelACModel::NucRatesLogProb,
                    &MultiGeneDiscSelACModel::UpdateMatrices, this);
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

    void MasterSendGvar()    {
        MasterSendGeneArray(*Gvararray);
    }

    void SlaveReceiveGvar()  {
        SlaveReceiveGeneArray(*Gvararray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetGvar((*Gvararray)[gene]);
        }
    }

    void MasterReceiveGvar() {
        MasterReceiveGeneArray(*Gvararray);
    }

    void SlaveSendGvar() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            (*Gvararray)[gene] = geneprocess[gene]->GetGvar();
        }
        SlaveSendGeneArray(*Gvararray);
    }

    void MasterSendGvarHyperParameters() {
        MasterSendGlobal(Gvarhypermean, Gvarhyperinvshape);
    }

    void SlaveReceiveGvarHyperParameters()   {
        SlaveReceiveGlobal(Gvarhypermean, Gvarhyperinvshape);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetGvarHyperParameters(Gvarhypermean, Gvarhyperinvshape);
        }
    }

    void MasterSendLogPsi()    {
        MasterSendGeneArray(*logpsiarray);
    }

    void SlaveReceiveLogPsi()  {
        SlaveReceiveGeneArray(*logpsiarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetLogPsi((*logpsiarray)[gene]);
        }
    }

    void MasterReceiveLogPsi() {
        MasterReceiveGeneArray(*logpsiarray);
    }

    void SlaveSendLogPsi() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            (*logpsiarray)[gene] = geneprocess[gene]->GetLogPsi();
        }
        SlaveSendGeneArray(*logpsiarray);
    }

    void MasterSendLogPsiHyperParameters() {
        MasterSendGlobal(logpsihypermean, logpsihypervar);
    }

    void SlaveReceiveLogPsiHyperParameters()   {
        SlaveReceiveGlobal(logpsihypermean, logpsihypervar);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetLogPsiHyperParameters(logpsihypermean, logpsihypervar);
        }
    }

    void MasterSendAAWeights()  {
        MasterSendGlobal(aaweights);
    }

    void SlaveReceiveAAWeights()    {
        SlaveReceiveGlobal(aaweights);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetAAWeights(aaweights);
        }
    }

    void MasterSendAADist() {
        if (! aadistmodel)  {
            MasterSendGlobal(wcom, wpol);
        }
        else    {
            MasterSendGlobal(aadist);
        }
    }

    void SlaveReceiveAADist()   {
        if (! aadistmodel)  {
            SlaveReceiveGlobal(wcom, wpol);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetW(wcom, wpol);
            }
        }
        else    {
            SlaveReceiveGlobal(aadist);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetAADist(aadist);
            }
        }
    }

    void SlaveSendPathSuffStat()    {
        // not useful: MoveGeneSelAC already computes those
        // GeneCollectComponentPathSuffStat();
        componentpathsuffstatbidimarray->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            componentpathsuffstatbidimarray->Add(geneprocess[gene]->GetComponentPathSuffStat());
        }
        SlaveSendAdditive(*componentpathsuffstatbidimarray);
    }

    void MasterReceivePathSuffStat()    {
        componentpathsuffstatbidimarray->Clear();
        MasterReceiveAdditive(*componentpathsuffstatbidimarray);
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

    void MasterSendNucRates() {
        MasterSendGlobal(nucrelrate, nucstat);
    }

    void SlaveReceiveNucRates() {
        SlaveReceiveGlobal(nucrelrate, nucstat);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates(nucrelrate, nucstat);
        }
    }

    // log probs

    void SlaveSendLogProbs() {
        GeneLogPrior = 0;
        lnL = 0;
        MeanStatEnt = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
            MeanStatEnt += geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanAAEntropy();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
        SlaveSendAdditive(MeanStatEnt);

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
        MeanStatEnt = 0;
        MasterReceiveAdditive(MeanStatEnt);
        MeanStatEnt /= GetTotNsite();

        moveTime = 0;
        mapTime = 0;
        MasterReceiveAdditive(moveTime);
        MasterReceiveAdditive(mapTime);
        moveTime /= (GetNprocs() - 1);
        mapTime /= (GetNprocs() - 1);
    }
};
