
#include "SelACOmegaModel.hpp"
#include "IIDBernoulliGamma.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "Permutation.hpp"

class MultiGeneSelACOmegaModel : public MultiGeneProbModel {
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

    // each gene has its own gamma distribution and its own pi
    int aadistmodel;
    int Gcat;
    double Ginvshape;

    // sharing the aadist parameters
    double wcom;
    double wpol;
    double wvol;
    // distance matrix between amino acids 
    vector<double> aadist;
    // prior probs for optimal amino acids
    vector<double> aaweight;
    OccupancySuffStat* aaoccupancy;

    double psihypermean;
    double psihyperinvshape;
    IIDGamma* psiarray;
    GammaSuffStat psihypersuffstat;

    std::vector<SelACOmegaModel *> geneprocess;

    double lnL;
    double GeneLogPrior;
    double MeanStatEnt;
    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    int blmode, nucmode, omegamode, omegaprior, modalprior;
    int aadistmode;

    Chrono totchrono;
    Chrono paramchrono;
    Chrono blchrono;
    Chrono aachrono;

    vector<double> aadistacc, aadisttot, gacc, gtot;
    int burnin;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneSelACOmegaModel(string indatafile, string intreefile, int inGcat, int inaadistmodel,
                                     int inblmode, int innucmode, int inaadistmode, int inomegamode,
                                     int inomegaprior, int inmodalprior, double indposompihypermean,
                                     double indposompihyperinvconc, int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs), nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {

        datafile = indatafile;
        treefile = intreefile;
        AllocateAlignments(datafile);

        Gcat = inGcat;
        aadistmodel = inaadistmodel;

        burnin = 0;

        blmode = inblmode;
        nucmode = innucmode;
        aadistmode = inaadistmode;
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

        wcom = grantham_wcom;
        wpol = grantham_wpol;
        wvol = grantham_wvol;

        aadist.assign(Naarr, 1.0);
        if (!aadistmodel)    {
            UpdateGrantham();
        }

        Ginvshape = 1.0;

        aaweight.assign(Naa, 1.0/Naa);
        aaoccupancy = new OccupancySuffStat(Naa);

        psihypermean = 10.0;
        psihyperinvshape = 0.5;
        double alpha = 1.0 / psihyperinvshape;
        double beta = alpha / psihypermean;
        psiarray = new IIDGamma(GetLocalNgene(), alpha, beta);

        genednds = new SimpleArray<double>(GetLocalNgene());

        lnL = 0;
        GeneLogPrior = 0;
        MeanStatEnt = 0;

        // aadist params

        if (!GetMyid()) {
            geneprocess.assign(0, (SelACOmegaModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (SelACOmegaModel *)0);

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
                    if (name == GeneName[index])    {
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
                    geneprocess[gene] = new SelACOmegaModel(
                        alivector[gene], tree, aadistmodel, aadistmode, omegamode, omegaprior, Gcat);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    geneprocess[gene] = new SelACOmegaModel(
                        GetLocalGeneName(gene), treefile, aadistmodel, aadistmode, omegamode, omegaprior, Gcat);
                }
            }

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->Allocate();
            }
        }

        aadistacc.assign(3,0);
        aadisttot.assign(3,0);
        gacc.assign(3,0);
        gtot.assign(3,0);
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

    void UpdateSelAC()  {
        if (! aadistmodel)   {
            UpdateGrantham();
        }
        double alpha = 1.0 / psihyperinvshape;
        double beta = alpha / psihypermean;
        psiarray->SetShape(alpha);
        psiarray->SetScale(beta);
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

        UpdateSelAC();
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

            if (omegamode == 1) {
                MasterSendOmegaHyperParameters();
                MasterSendOmega();
            }

            MasterSendPsiHyperParameters();
            MasterSendPsi();
            MasterSendSelAC();

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

        if (omegamode == 1) {
            SlaveReceiveOmegaHyperParameters();
            SlaveReceiveOmega();
        }

        SlaveReceivePsiHyperParameters();
        SlaveReceivePsi();
        SlaveReceiveSelAC();

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

            if (omegamode == 1) {
                MasterSendOmegaHyperParameters();
                MasterSendOmega();
            }

            MasterSendPsiHyperParameters();
            MasterSendPsi();
            MasterSendSelAC();
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

        if (omegamode == 1) {
            SlaveReceiveOmegaHyperParameters();
            SlaveReceiveOmega();
        }

        SlaveReceivePsiHyperParameters();
        SlaveReceivePsi();
        SlaveReceiveSelAC();

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
        if (omegamode != 3) {
            if (omegaprior == 0) {
                os << "meanomega\t";
                os << "varomega\t";
            } else {
                os << "npos\tposmean\tdposom_pi\tmeandposom\tinvshape\t";
            }
        }

        if (Gcat > 1)   {
            os << "ginvshape\t";
        }
        os << "meanpsi\tinvshape\t";
        if (! aadistmodel)   {
            os << "w_comp\t";
            os << "w_pol\t";
        }
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

        if (Gcat > 1)   {
            os << Ginvshape << '\t';
        }
        os << psihypermean << '\t' << psihyperinvshape << '\t';
        if (! aadistmodel)   {
            os << wcom << '\t';
            os << wpol << '\t';
        }
        os << GetVarAADist() << '\t';
        os << Random::GetEntropy(aaweight) << '\t';

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

    void TracePsi(ostream &os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << psiarray->GetVal(gene) << '\t';
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
            os << aaweight[a] << '\t';
        }
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                os << aadist[rrindex(a,b)] << '\t';
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

        os << "g  : ";
       for (size_t i=0; i<gacc.size(); i++)    {
          os << gacc[i] / gtot[i] << '\t';
       }
       os << '\n';

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

        is >> psihypermean >> psihyperinvshape;
        is >> *psiarray;
        is >> Ginvshape;
        if (! aadistmodel)   {
            is >> wcom >> wpol;
        }
        else    {
            is >> aadist;
        }
        is >> aaweight;

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

        os << psihypermean << '\t' << psihyperinvshape << '\t';
        os << *psiarray << '\t';

        os << Ginvshape << '\t';
        if (! aadistmodel)   {
            os << wcom << '\t' << wpol << '\t';
        }
        else    {
            os << aadist << '\t';
        }
        os << aaweight << '\t';

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

        if (omegamode != 3) {
            total += OmegaHyperLogPrior();
            total += OmegaLogPrior();
        }

        total += AALogPrior();

        return total;
    }

    double PsiHyperLogPrior() const {
        double total = 0;
        total -= 0.1 * psihypermean;
        total -= psihyperinvshape;
        return total;
    }

    double PsiLogPrior() const  {
        return psiarray->GetLogProb();
    }

    double GLogPrior() const    {
        return -Ginvshape;
    }

    double AADistLogPrior() const   {
        double total = 0;
        if (! aadistmodel)   {
            total -= log(wcom) + log(wpol);
        }
        else    {
            for (int i=0; i<Naarr; i++) {
                total -= aadist[i];
            }
        }
        return total;
    }

    double AALogPrior() const   {
        double total = 0;
        total += GLogPrior();
        total += PsiHyperLogPrior();
        total += PsiLogPrior();
        total += AADistLogPrior();
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
            if (modalprior && (omegahyperinvshape > 1.0))  {
                total += Random::INFPROB;
            }
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

    double PsiHyperSuffStatLogProb() const  {
        double ret = 0;
        double alpha = 1.0 / psihyperinvshape;
        double beta = alpha / psihypermean;
        ret = psihypersuffstat.GetLogProb(alpha, beta);
        return ret;
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

    // log prob for moving psi hyperparameters
    double PsiHyperLogProb() const { return PsiHyperLogPrior() + PsiHyperSuffStatLogProb(); }

    //-------------------
    // Moves
    //-------------------

    void MasterMove() override {
        totchrono.Start();
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            paramchrono.Start();

            aachrono.Start();

            MasterMoveSelAC();

            MasterReceivePsi();
            MovePsiHyperParameters();
            MasterSendPsiHyperParameters();

            aachrono.Stop();

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
            // gene aadist
            movechrono.Stop();

            SlaveMoveSelAC();
            MoveGeneSelAC();
            SlaveSendPsi();
            SlaveReceivePsiHyperParameters();

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

    void MoveGeneSelAC() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveSelAC();
            (*psiarray)[gene] = geneprocess[gene]->GetPsi();
        }
    }

    void ResampleAAWeight() {
        double total = 0;
        for (int i=0; i<Naa; i++)   {
            aaweight[i] = Random::sGamma(1 + aaoccupancy->GetVal(i));
            total += aaweight[i];
        }
        for (int i=0; i<Naa; i++)   {
            aaweight[i] /= total;
        }
    }

    void SlaveSendAAWeightSuffStat()    {
        aaoccupancy->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->UpdateOccupancies();
            aaoccupancy->Add(*geneprocess[gene]->GetOccupancies());
        }
        SlaveSendAdditive(*aaoccupancy);
    }

    void MasterReceiveAAWeightSuffStat()    {
        aaoccupancy->Clear();
        MasterReceiveAdditive(*aaoccupancy);
    }

    void MovePsiHyperParameters()   {

        psihypersuffstat.Clear();
        psihypersuffstat.AddSuffStat(*psiarray);

        ScalingMove(psihypermean, 1.0, 10,
                &MultiGeneSelACOmegaModel::PsiHyperLogProb, &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(psihypermean, 0.3, 10,
                &MultiGeneSelACOmegaModel::PsiHyperLogProb, &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(psihyperinvshape, 1.0, 10,
                &MultiGeneSelACOmegaModel::PsiHyperLogProb, &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(psihyperinvshape, 0.3, 10,
                &MultiGeneSelACOmegaModel::PsiHyperLogProb, &MultiGeneSelACOmegaModel::NoUpdate, this);

        double alpha = 1.0 / psihyperinvshape;
        double beta = alpha / psihypermean;
        psiarray->SetShape(alpha);
        psiarray->SetScale(beta);
    }

    void MasterMoveSelAC() {
        if (Gcat > 1)   {
            gacc[0] += MasterMoveG(3,1);
            gtot[0] ++;
            gacc[1] += MasterMoveG(3,0.1);
            gtot[1] ++;
            gacc[2] += MasterMoveG(3,0.01);
            gtot[2] ++;
        }
        if (! aadistmodel)  {
            aadistacc[0] += MasterMoveGranthamWeight(wcom,3,1);
            aadisttot[0] ++;
            aadistacc[1] += MasterMoveGranthamWeight(wcom,3,0.1);
            aadisttot[1] ++;
            aadistacc[0] += MasterMoveGranthamWeight(wpol,3,1);
            aadisttot[0] ++;
            aadistacc[1] += MasterMoveGranthamWeight(wpol,3,0.1);
            aadisttot[1] ++;
        }
        else    {
            aadistacc[0] += MasterMoveAADist(3,1);
            aadisttot[0] ++;
            aadistacc[1] += MasterMoveAADist(3,0.1);
            aadisttot[1] ++;
            aadistacc[2] += MasterMoveAADist(3,0.01);
            aadisttot[2] ++;
        }

        AAPsiCompMove(1.0, 10);
        MasterSendAADist();
        MasterSendPsi();

        MasterReceiveAAWeightSuffStat();
        ResampleAAWeight();
        MasterSendAAWeight();
    }

    void SlaveMoveSelAC()  {
        if (Gcat > 1)   {
            SlaveMoveG(3);
            SlaveMoveG(3);
            SlaveMoveG(3);
        }
        if (! aadistmodel)  {
            SlaveMoveGranthamWeight(3);
            SlaveMoveGranthamWeight(3);
            SlaveMoveGranthamWeight(3);
            SlaveMoveGranthamWeight(3);
        }
        else    {
            SlaveMoveAADist(3);
            SlaveMoveAADist(3);
            SlaveMoveAADist(3);
        }

        SlaveReceiveAADist();
        SlaveReceivePsi();

        SlaveSendAAWeightSuffStat();
        SlaveReceiveAAWeight();
    }

    double AAPsiCompMove(double tuning, int nrep)    {

        double nacc = 0;
        double ntot = 0;

        for (int rep=0; rep<nrep; rep++)    {
            double delta = - AALogPrior();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            for (int i=0; i<GetNgene(); i++)    {
                (*psiarray)[i] /= e;
            }
            for (int i=0; i<Naarr; i++) {
                aadist[i] *= e;
            }
            delta += AALogPrior();
            delta += (Naarr-GetNgene())*m;

            int accepted = (log(Random::Uniform()) < delta);
            if (accepted) {
                nacc++;
            } else {
                for (int i=0; i<GetNgene(); i++)    {
                    (*psiarray)[i] *= e;
                }
                for (int i=0; i<Naarr; i++) {
                    aadist[i] /= e;
                }
            }
            ntot++;
        }

        UpdateSelAC();

        return nacc / ntot;
    }

    double MasterMoveG(int nrep, double tuning)   {

        double acc = 0;
        double tot = 0;

        double logprob1 = MasterReceiveAALogProb();

        for (int rep=0; rep<nrep; rep++)    {

            // propose move
            double delta = -GLogPrior();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            Ginvshape *= e;
            delta += GLogPrior();

            // send new value
            MasterSendG();

            // receive log probs
            double logprob2 = MasterReceiveAALogProb();

            // decide
            delta += logprob2 - logprob1;
            delta += m;
            int accept = (log(Random::Uniform()) < delta);

            // accept or restore
            if (accept) {
                acc++;
                logprob1 = logprob2;
            }
            else    {
                Ginvshape /= e;
            }
            tot++;
        }
        MasterSendG();
        return acc / tot;
    }

    void SlaveMoveG(int nrep)   {
        SlaveSendAALogProb();
        for (int rep=0; rep<nrep; rep++)    {
            SlaveReceiveG();
            SlaveSendAALogProb();
        }
        SlaveReceiveG();
    }

    double MasterMoveAADist(int nrep, double tuning) {

        double acc = 0;
        double tot = 0;

        int Npair = Naa/2;

        vector<double> logprob1(Naa,0);
        vector<double> logprob2(Naa,0);
        vector<double> deltalogprob(Npair,0);
        vector<double> bk(Npair,0);
        int aa[Naa];

        // MasterSendAADist();
        MasterReceiveAALogProbs(logprob1);

        for (int rep=0; rep<nrep; rep++)    {

            // propose move
            Random::DrawFromUrn(aa,Naa,Naa);

            for (int i=0; i<Npair; i++) {
                double& d = aadist[rrindex(aa[2*i], aa[2*i+1])];
                bk[i] = d;
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                deltalogprob[i] = m;
                deltalogprob[i] += d;
                deltalogprob[i] -= logprob1[aa[2*i]] + logprob1[aa[2*i+1]];
                d *= e;
            }

            // send proposed moves
            MasterSendAADist();
            // receive new logprobs
            MasterReceiveAALogProbs(logprob2);

            // make decisions
            for (int i=0; i<Npair; i++) {
                double& d = aadist[rrindex(aa[2*i], aa[2*i+1])];
                deltalogprob[i] -= d;
                deltalogprob[i] += logprob2[aa[2*i]] + logprob2[aa[2*i+1]];

                int accept = (log(Random::Uniform()) < deltalogprob[i]);
                if (accept) {
                    acc++;
                    logprob1[aa[2*i]] = logprob2[aa[2*i]];
                    logprob1[aa[2*i+1]] = logprob2[aa[2*i+1]];
                }
                else    {
                    aadist[rrindex(aa[2*i], aa[2*i+1])] = bk[i];
                }
                tot++;
            }
        }

        MasterSendAADist();
        return acc / tot;
    }

    void SlaveMoveAADist(int nrep)  {
        // SlaveReceiveAADist();
        SlaveSendAALogProbs();
        for (int rep=0; rep<nrep; rep++)    {
            SlaveReceiveAADist();
            SlaveSendAALogProbs();

            /*
            SlaveReceiveAADist();
            SlaveSendAALogProbs();
            */
        }
        SlaveReceiveAADist();
    }

    double MasterMoveGranthamWeight(double& w, int nrep, double tuning) {
        double acc = 0;
        double tot = 0;
        for (int rep=0; rep<nrep; rep++)    {
        }
        return acc/tot;
    }

    void SlaveMoveGranthamWeight(int nrep)    {
        SlaveSendAALogProb();
        for (int rep=0; rep<nrep; rep++)    {
            SlaveReceiveAADist();
            SlaveSendAALogProb();
        }
        SlaveReceiveAADist();
    }


    double MasterReceiveAALogProb() {
        double ret = 0;
        MasterReceiveAdditive(ret);
        return ret;
    }

    void SlaveSendAALogProb()   {
        double ret = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            ret += geneprocess[gene]->PathSuffStatLogProb();
        }
        SlaveSendAdditive(ret);
    }

    void MasterReceiveAALogProbs(vector<double>& logprobs)  {
        for (int a=0; a<Naa; a++)   {
            logprobs[a] = 0;
        }
        MasterReceiveAdditive(logprobs);
    }

    void SlaveSendAALogProbs()  {
        vector<double> logprobs(Naa,0);
        vector<double> tmp(Naa,0);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->AAPathSuffStatLogProbs(tmp);
            for (int a=0; a<Naa; a++)   {
                logprobs[a] += tmp[a];
            }
        }
        SlaveSendAdditive(logprobs);
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

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneSelACOmegaModel::LambdaHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneSelACOmegaModel::LambdaHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, 10);
            BranchLengthsHyperScalingMove(0.3, 10);
        }

        ScalingMove(blhyperinvshape, 1.0, 10,
                    &MultiGeneSelACOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10,
                    &MultiGeneSelACOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);

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
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10,
                    &MultiGeneSelACOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneSelACOmegaModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void MoveOmegaHyperParameters() {
        if (omegaprior == 0) {
            omegahypersuffstat.Clear();
            omegahypersuffstat.AddSuffStat(*omegaarray);

            ScalingMove(omegahypermean, 1.0, 10,
                        &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                        &MultiGeneSelACOmegaModel::NoUpdate, this);
            ScalingMove(omegahypermean, 0.3, 10,
                        &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                        &MultiGeneSelACOmegaModel::NoUpdate, this);
            if (modalprior) {
                SlidingMove(omegahyperinvshape,1.0,10,0,1.0,&MultiGeneSelACOmegaModel::OmegaHyperLogProb,&MultiGeneSelACOmegaModel::NoUpdate,this);
                SlidingMove(omegahyperinvshape,0.3,10,0,1.0,&MultiGeneSelACOmegaModel::OmegaHyperLogProb,&MultiGeneSelACOmegaModel::NoUpdate,this);
                SlidingMove(omegahyperinvshape,0.1,10,0,1.0,&MultiGeneSelACOmegaModel::OmegaHyperLogProb,&MultiGeneSelACOmegaModel::NoUpdate,this);
            } else  {
                ScalingMove(omegahyperinvshape, 1.0, 10,
                            &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                            &MultiGeneSelACOmegaModel::NoUpdate, this);
                ScalingMove(omegahyperinvshape, 0.3, 10,
                            &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                            &MultiGeneSelACOmegaModel::NoUpdate, this);
            }

            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            omegaarray->SetShape(alpha);
            omegaarray->SetScale(beta);
        } else {
            dposomhypersuffstat.Clear();
            dposomhypersuffstat.AddSuffStat(*dposomarray);

            ScalingMove(dposomhypermean, 1.0, 10,
                        &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                        &MultiGeneSelACOmegaModel::NoUpdate, this);
            ScalingMove(dposomhypermean, 0.3, 10,
                        &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                        &MultiGeneSelACOmegaModel::NoUpdate, this);

            if (modalprior) {
                SlidingMove(dposomhyperinvshape,1.0,10,0,1.0,&MultiGeneSelACOmegaModel::OmegaHyperLogProb,&MultiGeneSelACOmegaModel::NoUpdate,this);
                SlidingMove(dposomhyperinvshape,0.3,10,0,1.0,&MultiGeneSelACOmegaModel::OmegaHyperLogProb,&MultiGeneSelACOmegaModel::NoUpdate,this);
                SlidingMove(dposomhyperinvshape,0.1,10,0,1.0,&MultiGeneSelACOmegaModel::OmegaHyperLogProb,&MultiGeneSelACOmegaModel::NoUpdate,this);
            } else  {
                ScalingMove(dposomhyperinvshape, 1.0, 10,
                            &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                            &MultiGeneSelACOmegaModel::NoUpdate, this);
                ScalingMove(dposomhyperinvshape, 0.3, 10,
                            &MultiGeneSelACOmegaModel::OmegaHyperLogProb,
                            &MultiGeneSelACOmegaModel::NoUpdate, this);
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
        double postalpha = Random::sGamma(pialpha + n1);
        double postbeta = Random::sGamma(pibeta + n0);
        dposompi = postalpha / (postalpha + postbeta);
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

    void MasterSendPsi()    {
        MasterSendGeneArray(*psiarray);
    }

    void SlaveReceivePsi()  {
        SlaveReceiveGeneArray(*psiarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetPsi((*psiarray)[gene]);
        }
    }

    void MasterReceivePsi() {
        MasterReceiveGeneArray(*psiarray);
    }

    void SlaveSendPsi() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            (*psiarray)[gene] = geneprocess[gene]->GetPsi();
        }
        SlaveSendGeneArray(*psiarray);
    }

    void MasterSendPsiHyperParameters() {
        MasterSendGlobal(psihypermean, psihyperinvshape);
    }

    void SlaveReceivePsiHyperParameters()   {
        SlaveReceiveGlobal(psihypermean, psihyperinvshape);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetPsiHyperParameters(psihypermean, psihyperinvshape);
        }
    }

    void MasterSendG()  {
        MasterSendGlobal(Ginvshape);
    }

    void SlaveReceiveG()    {
        SlaveReceiveGlobal(Ginvshape);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetG(Ginvshape);
        }
    }

    void MasterSendW()  {
        MasterSendGlobal(wcom,wpol);
    }

    void SlaveReceiveW()    {
        SlaveReceiveGlobal(wcom,wpol);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetW(wcom, wpol);
        }
    }

    void MasterSendAAWeight()  {
        MasterSendGlobal(aaweight);
    }

    void SlaveReceiveAAWeight()    {
        SlaveReceiveGlobal(aaweight);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetAAWeight(aaweight);
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

    void MasterSendSelAC() {
        MasterSendGlobal(Ginvshape, aaweight);
        if (! aadistmodel)  {
            MasterSendGlobal(wcom, wpol);
        }
        else    {
            MasterSendGlobal(aadist);
        }
    }

    void SlaveReceiveSelAC()  {
        SlaveReceiveGlobal(Ginvshape, aaweight);
        if (! aadistmodel)  {
            SlaveReceiveGlobal(wcom, wpol);
        }
        else    {
            SlaveReceiveGlobal(aadist);
        }
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetG(Ginvshape);
            geneprocess[gene]->SetAAWeight(aaweight);
            if (! aadistmodel)  {
                geneprocess[gene]->SetW(wcom, wpol);
            }
            else    {
                geneprocess[gene]->SetAADist(aadist);
            }
        }
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
