
#include "AAMutSelSparseOmegaModel.hpp"
#include "IIDBernoulliGamma.hpp"
#include "IIDBernoulliCauchy.hpp"
#include "IIDBeta.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "Permutation.hpp"

class MultiGeneAAMutSelSparseOmegaModel : public MultiGeneProbModel {
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
    double maxdposom;
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
    IIDBernoulliGamma *gammadposomarray;
    IIDBernoulliCauchy *cauchydposomarray;
    BernoulliGammaSuffStat dposomhypersuffstat;

    double epsilonhypermean;
    double epsilonhyperinvconc;
    IIDBeta* epsilonarray;
    BetaSuffStat epsilonhypersuffstat;

    double pihypermean;
    double pihyperinvconc;
    IIDBeta* piarray;
    BetaSuffStat pihypersuffstat;

    std::vector<AAMutSelSparseOmegaModel *> geneprocess;

    double lnL;
    double GeneLogPrior;
    double MeanStatEnt;
    double MeanWidth;
    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    int blmode, nucmode, omegamode, omegaprior, modalprior;

    Chrono totchrono;
    Chrono paramchrono;
    Chrono blchrono;

    int burnin;
    int chainsize;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneAAMutSelSparseOmegaModel(string indatafile, string intreefile, 
                                     int inblmode, int innucmode,
                                     int inomegamode, int inomegaprior, int inmodalprior,
                                     double indposompihypermean, double indposompihyperinvconc,
                                     double inmaxdposom, 
                                     double inepsilonhypermean, double inepsilonhyperinvconc,
                                     double inpihypermean, double inpihyperinvconc,
                                     int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs), nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {

        datafile = indatafile;
        treefile = intreefile;
        AllocateAlignments(datafile);

        burnin = 10;
        chainsize = 0;

        blmode = inblmode;
        nucmode = innucmode;
        omegamode = inomegamode;
        omegaprior = inomegaprior;
        modalprior = inmodalprior;
        dposompihypermean = indposompihypermean;
        dposompihyperinvconc = indposompihyperinvconc;
        maxdposom = inmaxdposom;

        epsilonhypermean = inepsilonhypermean;
        epsilonhyperinvconc = inepsilonhyperinvconc;
        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;

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

    void SetChainSize(int insize) {
        chainsize = insize;
        if (myid) {
            for (int gene=0; gene<GetNgene(); gene++)   {
                geneprocess[gene]->SetChainSize(insize);
            }
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
        } else if ((omegaprior == 1) || (omegaprior == 2)) {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            gammadposomarray = new IIDBernoulliGamma(GetLocalNgene(), dposompi, alpha, beta);
            cauchydposomarray = 0;
            for (int i = 0; i < GetLocalNgene(); i++) {
                (*gammadposomarray)[i] = 0;
            }
            if (maxdposom)   {
                for (int i = 0; i < GetLocalNgene(); i++) {
                    if ((*gammadposomarray)[i] > maxdposom)    {
                        (*gammadposomarray)[i] = maxdposom;
                    }
                }
            }
        } else if (omegaprior == 3) {
            double gamma = 1.0 / dposomhyperinvshape;
            cauchydposomarray = new IIDBernoulliCauchy(GetLocalNgene(), dposompi, gamma);
            gammadposomarray = 0;
            for (int i = 0; i < GetLocalNgene(); i++) {
                (*cauchydposomarray)[i] = 0;
            }
        } else {
            cerr << "error: unrecognized omega prior\n";
            exit(1);
        }

        genednds = new SimpleArray<double>(GetLocalNgene());

        if (epsilonhyperinvconc)   {
            double alpha = epsilonhypermean / epsilonhyperinvconc;
            double beta = (1-epsilonhypermean) / epsilonhyperinvconc;
            epsilonarray = new IIDBeta(GetNgene(), alpha, beta);
        }
        if (pihyperinvconc)   {
            double alpha = pihypermean / pihyperinvconc;
            double beta = (1-pihypermean) / pihyperinvconc;
            piarray = new IIDBeta(GetNgene(), alpha, beta);
        }

        lnL = 0;
        GeneLogPrior = 0;
        MeanStatEnt = 0;
        MeanWidth = 0;


        if (!GetMyid()) {
            geneprocess.assign(0, (AAMutSelSparseOmegaModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (AAMutSelSparseOmegaModel *)0);

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
                    double epsilon = epsilonhyperinvconc ? -1.0 : epsilonhypermean;
                    double pi = pihyperinvconc ? -1.0 : pihypermean;
                    geneprocess[gene] = new AAMutSelSparseOmegaModel(
                        alivector[gene], tree, omegamode, omegaprior, 3, 0, epsilon, pi);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    double epsilon = epsilonhyperinvconc ? -1.0 : epsilonhypermean;
                    double pi = pihyperinvconc ? -1.0 : pihypermean;
                    geneprocess[gene] = new AAMutSelSparseOmegaModel(GetLocalGeneName(gene),
                            treefile, omegamode, omegaprior, 3, 0, epsilon, pi);
                }
            }

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->SetMaxDPosOm(maxdposom);
                geneprocess[gene]->Allocate();
            }
        }
    }

    void FastUpdate() {
        UpdateBranchLengths();
        UpdateNucRates();
        if (epsilonhyperinvconc)   {
            UpdateEpsilon();
        }
        if (pihyperinvconc)    {
            UpdatePi();
        }
        UpdateOmegaMixture();
    }

    void UpdateBranchLengths()  {
        branchlength->SetScale(lambda);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }
    }

    void UpdateNucRates()   {
        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    void UpdateEpsilon() {
        double alpha = epsilonhypermean / epsilonhyperinvconc;
        double beta = (1-epsilonhypermean) / epsilonhyperinvconc;
        epsilonarray->SetAlpha(alpha);
        epsilonarray->SetBeta(beta);
    }

    void UpdatePi() {
        double alpha = pihypermean / pihyperinvconc;
        double beta = (1-pihypermean) / pihyperinvconc;
        piarray->SetAlpha(alpha);
        piarray->SetBeta(beta);
    }

    void UpdateOmegaMixture()   {
        if (omegaprior == 0) {
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            omegaarray->SetShape(alpha);
            omegaarray->SetScale(beta);
        } else if ((omegaprior == 1) || (omegaprior == 2)) {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            gammadposomarray->SetPi(dposompi);
            gammadposomarray->SetShape(alpha);
            gammadposomarray->SetScale(beta);
        } else if (omegaprior == 3) {
            cauchydposomarray->SetPi(dposompi);
            cauchydposomarray->SetGamma(1.0 / dposomhyperinvshape);
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

            if (epsilonhyperinvconc)   {
                MasterSendEpsilonHyperParameters();
                MasterSendEpsilon();
            }

            if (pihyperinvconc)    {
                MasterSendPiHyperParameters();
                MasterSendPi();
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

        if (epsilonhyperinvconc)   {
            SlaveReceiveEpsilonHyperParameters();
            SlaveReceiveEpsilon();
        }

        if (pihyperinvconc)    {
            SlaveReceivePiHyperParameters();
            SlaveReceivePi();
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

            if (epsilonhyperinvconc)   {
                MasterSendEpsilonHyperParameters();
                MasterSendEpsilon();
            }

            if (pihyperinvconc)    {
                MasterSendPiHyperParameters();
                MasterSendPi();
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

        if (epsilonhyperinvconc)   {
            SlaveReceiveEpsilonHyperParameters();
            SlaveReceiveEpsilon();
        }

        if (pihyperinvconc)    {
            SlaveReceivePiHyperParameters();
            SlaveReceivePi();
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

    void TraceEpsilon(ostream& os) const   {
        for (int gene = 0; gene < Ngene; gene++) {
            os << epsilonarray->GetVal(gene) << '\t';
        }
        os << '\n';
    }

    void TracePi(ostream& os) const {
        for (int gene = 0; gene < Ngene; gene++) {
            os << piarray->GetVal(gene) << '\t';
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
            } else if ((omegaprior == 1) || (omegaprior == 2))  {
                os << "npos\tposmean\tdposom_pi\tmeandposom\tinvshape\t";
            } else  {
                os << "npos\tposmean\tdposom_pi\tcauchyinvshape\t";
            }
        }
        if (epsilonhyperinvconc)    {
            os << "eps\tinvconc\t";
        }
        if (pihyperinvconc) {
            os << "pi\tinvconc\t";
        }
        os << "width\t";
        os << "aastatent\t";
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
            } else if ((omegaprior == 1) || (omegaprior == 2))  {
                os << gammadposomarray->GetNpos() << '\t';
                os << gammadposomarray->GetPosMean() << '\t';
                os << dposompi << '\t' << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
            } else if (omegaprior == 3) {
                os << cauchydposomarray->GetNpos() << '\t';
                os << cauchydposomarray->GetPosMean() << '\t';
                os << dposompi << '\t' << dposomhyperinvshape << '\t';
            }
        }
        if (epsilonhyperinvconc)    {
            os << epsilonhypermean << '\t' << epsilonhyperinvconc << '\t';
        }
        if (pihyperinvconc) {
            os << pihypermean << '\t' << pihyperinvconc << '\t';
        }
        os << MeanWidth << '\t';
        os << MeanStatEnt << '\t';
        os << Random::GetEntropy(nucrelratehypercenter) << '\t' << nucrelratehyperinvconc << '\t';
        os << Random::GetEntropy(nucstathypercenter) << '\t' << nucstathyperinvconc << '\n';
        os.flush();
    }

    void TraceOmega(ostream &os) const {
        if (omegaprior == 0) {
            for (int gene = 0; gene < Ngene; gene++) {
                os << omegaarray->GetVal(gene) << '\t';
            }
        } else if (omegaprior == 1) {
            for (int gene = 0; gene < Ngene; gene++) {
                os << 1.0 + gammadposomarray->GetVal(gene) << '\t';
            }
        } else if (omegaprior == 2) {
            for (int gene = 0; gene < Ngene; gene++) {
                os << exp(gammadposomarray->GetVal(gene)) << '\t';
            }
        } else if (omegaprior == 3) {
            for (int gene = 0; gene < Ngene; gene++) {
                os << 1.0 + cauchydposomarray->GetVal(gene) << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

    void Monitor(ostream &os) const override {
        os << totchrono.GetTime() << '\t' << paramchrono.GetTime() << '\t' << blchrono.GetTime() << '\n';
        os << "prop time in param moves: " << paramchrono.GetTime() / totchrono.GetTime() << '\n';
        os << "sub prop time in bl moves   : " << blchrono.GetTime() / paramchrono.GetTime() << '\n';
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
            } else if ((omegaprior == 1) || (omegaprior == 2))  {
                is >> dposompi;
                is >> dposomhypermean;
                is >> dposomhyperinvshape;
                is >> *gammadposomarray;
            } else if (omegaprior == 3) {
                is >> dposompi;
                is >> dposomhyperinvshape;
                is >> *cauchydposomarray;
            }
        }
        if (epsilonhyperinvconc)    {
            is >> epsilonhypermean >> epsilonhyperinvconc;
            is >> *epsilonarray;
        }
        if (pihyperinvconc) {
            is >> pihypermean >> pihyperinvconc;
            is >> *piarray;
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
            } else if ((omegaprior == 1) || (omegaprior == 2))  {
                os << dposompi << '\t';
                os << dposomhypermean << '\t';
                os << dposomhyperinvshape << '\t';
                os << *gammadposomarray << '\t';
            } else if (omegaprior == 3) {
                os << dposompi << '\t';
                os << dposomhyperinvshape << '\t';
                os << *cauchydposomarray << '\t';
            }
        }
        if (epsilonhyperinvconc)    {
            os << epsilonhypermean << '\t' << epsilonhyperinvconc << '\t';
            os << *epsilonarray << '\t';
        }
        if (pihyperinvconc) {
            os << pihypermean << '\t' << pihyperinvconc << '\t';
            os << *piarray;
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

        if (omegamode != 3) {
            total += OmegaHyperLogPrior();
        }

        if (epsilonhyperinvconc)    {
            total += EpsilonHyperLogPrior();
        }
        if (pihyperinvconc) {
            total += PiHyperLogPrior();
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
        } else if ((omegaprior == 1) || (omegaprior == 2)) {
            double pialpha = dposompihypermean / dposompihyperinvconc;
            double pibeta = (1 - dposompihypermean) / dposompihyperinvconc;
            total += (pialpha - 1) * log(1.0 - dposompi) + (pibeta - 1) * log(dposompi);
            total -= dposomhypermean;
            total -= dposomhyperinvshape;
            if (modalprior && (dposomhyperinvshape > 1.0))  {
                total += Random::INFPROB;
            }
        } else if (omegaprior == 3) {
            double pialpha = dposompihypermean / dposompihyperinvconc;
            double pibeta = (1 - dposompihypermean) / dposompihyperinvconc;
            total += (pialpha - 1) * log(1.0 - dposompi) + (pibeta - 1) * log(dposompi);
            total -= dposomhyperinvshape;
        }
        return total;
    }

    double EpsilonHyperLogPrior() const {
        return -epsilonhyperinvconc;
    }

    double PiHyperLogPrior() const {
        return -pihyperinvconc;
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
        } else if ((omegaprior == 1) || (omegaprior == 2))  {
            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            ret = dposomhypersuffstat.GetLogProb(dposompi, alpha, beta);
        } else  {
            cerr << "error: prior not valid in omega hyper suff stat log prob\n";
            exit(1);
        }
        return ret;
    }

    double EpsilonHyperSuffStatLogProb() const  {
        double alpha = epsilonhypermean / epsilonhyperinvconc;
        double beta = (1-epsilonhypermean) / epsilonhyperinvconc;
        return epsilonhypersuffstat.GetLogProb(alpha,beta);
    }

    double PiHyperSuffStatLogProb() const   {
        double alpha = pihypermean / pihyperinvconc;
        double beta = (1-pihypermean) / pihyperinvconc;
        return pihypersuffstat.GetLogProb(alpha,beta);
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
    double GammaOmegaHyperLogProb() const { return OmegaHyperLogPrior() + OmegaHyperSuffStatLogProb(); }

    double CauchyOmegaHyperLogProb() const { return OmegaHyperLogPrior() + cauchydposomarray->GetPosLogProb();}

    void CauchyOmegaUpdate() {
        cauchydposomarray->SetPi(dposompi);
        cauchydposomarray->SetGamma(1.0 / dposomhyperinvshape);
    }

    double EpsilonHyperLogProb() const { return EpsilonHyperLogPrior() + EpsilonHyperSuffStatLogProb();}

    double PiHyperLogProb() const { return PiHyperLogPrior() + PiHyperSuffStatLogProb();}

    //-------------------
    // Moves
    //-------------------

    void MasterMove() override {
        totchrono.Start();
        int nrep = 20;

        for (int rep = 0; rep < nrep; rep++) {
            paramchrono.Start();

            if (epsilonhyperinvconc)    {
                MasterReceiveEpsilon();
                MoveEpsilonHyperParameters();
                MasterSendEpsilonHyperParameters();
            }

            if (pihyperinvconc) {
                MasterReceivePi();
                MovePiHyperParameters();
                MasterSendPiHyperParameters();
            }

            if ((chainsize >= burnin) && (omegamode != 3)) {
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

        chainsize++;
    }

    // slave move
    void SlaveMove() override {
        movechrono.Start();
        mapchrono.Start();
        GeneResampleSub(1.0);
        mapchrono.Stop();
        movechrono.Stop();

        int nrep = 20;

        for (int rep = 0; rep < nrep; rep++) {
            movechrono.Start();
            GeneCollectPathSuffStat();
            MoveGeneAA();
            movechrono.Stop();

            if (epsilonhyperinvconc)    {
                SlaveSendEpsilon();
                SlaveReceiveEpsilonHyperParameters();
            }

            if (pihyperinvconc) {
                SlaveSendPi();
                SlaveReceivePiHyperParameters();
            }

            if ((chainsize >= burnin) && (omegamode != 3)) {
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
        }
    }

    void MoveGeneOmegas() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveOmega();
            if (omegaprior == 0) {
                (*omegaarray)[gene] = geneprocess[gene]->GetOmega();
            } else if (omegaprior == 1) {
                (*gammadposomarray)[gene] = geneprocess[gene]->GetOmega() - 1;
            } else if (omegaprior == 2) {
                (*gammadposomarray)[gene] = log(geneprocess[gene]->GetOmega());
            } else if (omegaprior == 3) {
                (*cauchydposomarray)[gene] = geneprocess[gene]->GetOmega() - 1;
            }
        }
    }

    void MoveGeneAA() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveAA(10);
            if (epsilonhyperinvconc)    {
                (*epsilonarray)[gene] = geneprocess[gene]->GetEpsilon();
            }
            if (pihyperinvconc) {
                (*piarray)[gene] = geneprocess[gene]->GetPi();
            }
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

    void ResampleBranchLengths() { branchlength->GibbsResample(*lengthpathsuffstatarray); }

    void MoveLambda() {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneAAMutSelSparseOmegaModel::LambdaHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneAAMutSelSparseOmegaModel::LambdaHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, 10);
            BranchLengthsHyperScalingMove(0.3, 10);
        }

        ScalingMove(blhyperinvshape, 1.0, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);

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
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelSparseOmegaModel::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);

        UpdateNucRates();
    }

    void MoveEpsilonHyperParameters()   {
        epsilonhypersuffstat.Clear();
        epsilonarray->AddSuffStat(epsilonhypersuffstat);
        SlidingMove(epsilonhypermean, 1.0, 10, 0, 1.0, &MultiGeneAAMutSelSparseOmegaModel::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        SlidingMove(epsilonhypermean, 0.1, 10, 0, 1.0, &MultiGeneAAMutSelSparseOmegaModel::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(epsilonhyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseOmegaModel::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(epsilonhyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseOmegaModel::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);

        double alpha = epsilonhypermean / epsilonhyperinvconc;
        double beta = (1.0 - epsilonhypermean) / epsilonhyperinvconc;
        epsilonarray->SetAlpha(alpha);
        epsilonarray->SetBeta(beta);
    }

    void MovePiHyperParameters()   {
        pihypersuffstat.Clear();
        piarray->AddSuffStat(pihypersuffstat);
        SlidingMove(pihypermean, 1.0, 10, 0, 1.0, &MultiGeneAAMutSelSparseOmegaModel::PiHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        SlidingMove(pihypermean, 0.1, 10, 0, 1.0, &MultiGeneAAMutSelSparseOmegaModel::PiHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(pihyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseOmegaModel::PiHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
        ScalingMove(pihyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseOmegaModel::PiHyperLogProb,
                &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);

        double alpha = pihypermean / pihyperinvconc;
        double beta = (1.0 - pihypermean) / pihyperinvconc;
        piarray->SetAlpha(alpha);
        piarray->SetBeta(beta);
    }

    void MoveOmegaHyperParameters() {
        if (omegaprior == 0) {
            omegahypersuffstat.Clear();
            omegahypersuffstat.AddSuffStat(*omegaarray);

            ScalingMove(omegahypermean, 1.0, 10,
                        &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                        &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
            ScalingMove(omegahypermean, 0.3, 10,
                        &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                        &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
            if (modalprior) {
                SlidingMove(omegahyperinvshape,1.0,10,0,1.0,&MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,&MultiGeneAAMutSelSparseOmegaModel::NoUpdate,this);
                SlidingMove(omegahyperinvshape,0.3,10,0,1.0,&MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,&MultiGeneAAMutSelSparseOmegaModel::NoUpdate,this);
                SlidingMove(omegahyperinvshape,0.1,10,0,1.0,&MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,&MultiGeneAAMutSelSparseOmegaModel::NoUpdate,this);
            } else  {
                ScalingMove(omegahyperinvshape, 1.0, 10,
                            &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                            &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
                ScalingMove(omegahyperinvshape, 0.3, 10,
                            &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                            &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
            }

            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            omegaarray->SetShape(alpha);
            omegaarray->SetScale(beta);

        } else if ((omegaprior == 1) || (omegaprior == 2))  {
            dposomhypersuffstat.Clear();
            dposomhypersuffstat.AddSuffStat(*gammadposomarray);

            ScalingMove(dposomhypermean, 1.0, 10,
                        &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                        &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
            ScalingMove(dposomhypermean, 0.3, 10,
                        &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                        &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);

            if (modalprior) {
                SlidingMove(dposomhyperinvshape,1.0,10,0,1.0,&MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,&MultiGeneAAMutSelSparseOmegaModel::NoUpdate,this);
                SlidingMove(dposomhyperinvshape,0.3,10,0,1.0,&MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,&MultiGeneAAMutSelSparseOmegaModel::NoUpdate,this);
                SlidingMove(dposomhyperinvshape,0.1,10,0,1.0,&MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,&MultiGeneAAMutSelSparseOmegaModel::NoUpdate,this);
            } else  {
                ScalingMove(dposomhyperinvshape, 1.0, 10,
                            &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                            &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
                ScalingMove(dposomhyperinvshape, 0.3, 10,
                            &MultiGeneAAMutSelSparseOmegaModel::GammaOmegaHyperLogProb,
                            &MultiGeneAAMutSelSparseOmegaModel::NoUpdate, this);
            }

            // if (chainsize >= burnin) {
            if (dposompihyperinvconc) {
                ResampleDPosOmPi();
            }
            // }

            double alpha = 1.0 / dposomhyperinvshape;
            double beta = alpha / dposomhypermean;
            gammadposomarray->SetPi(dposompi);
            gammadposomarray->SetShape(alpha);
            gammadposomarray->SetScale(beta);

        } else if (omegaprior == 3) {

            ScalingMove(dposomhyperinvshape, 1.0, 10,
                        &MultiGeneAAMutSelSparseOmegaModel::CauchyOmegaHyperLogProb,
                        &MultiGeneAAMutSelSparseOmegaModel::CauchyOmegaUpdate, this);
            ScalingMove(dposomhyperinvshape, 0.3, 10,
                        &MultiGeneAAMutSelSparseOmegaModel::CauchyOmegaHyperLogProb,
                        &MultiGeneAAMutSelSparseOmegaModel::CauchyOmegaUpdate, this);

            // if (chainsize >= burnin) {
            if (dposompihyperinvconc) {
                ResampleDPosOmPi();
            }
            // }

            CauchyOmegaUpdate();
        }
    }

    int GetNpos() const {
        if (! omegaprior)   {
            cerr << "error in get npos: not under mixture model\n";
            exit(1);
        }
        int ret = 0;
        if ((omegaprior == 1) || (omegaprior == 2)) {
            ret = gammadposomarray->GetNpos();
        }
        else    {
            ret = cauchydposomarray->GetNpos();
        }
        return ret;
    }

    void ResampleDPosOmPi() {
        int n1 = GetNpos();
        int n0 = Ngene- n1;
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
        UpdateBranchLengths();
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
        UpdateNucRates();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter,
                                                          nucrelratehyperinvconc,
                                                          nucstathypercenter, nucstathyperinvconc);
        }
    }

    // omega (and hyperparameters)

    void SlaveSendOmega() {
        if (omegaprior == 0) {
            SlaveSendGeneArray(*omegaarray);
        } else if (omegaprior == 1) {
            SlaveSendGeneArray(*gammadposomarray);
        } else if (omegaprior == 2) {
            SlaveSendGeneArray(*gammadposomarray);
        } else if (omegaprior == 3) {
            SlaveSendGeneArray(*cauchydposomarray);
        }
    }

    void MasterReceiveOmega() {
        if (omegaprior == 0) {
            MasterReceiveGeneArray(*omegaarray);
        } else if ((omegaprior == 1) || (omegaprior == 2))  {
            MasterReceiveGeneArray(*gammadposomarray);
        } else if (omegaprior == 3) {
            MasterReceiveGeneArray(*cauchydposomarray);
        }
    }

    void MasterSendOmega() {
        if (omegaprior == 0) {
            MasterSendGeneArray(*omegaarray);
        } else if ((omegaprior == 1) || (omegaprior == 2))  {
            MasterSendGeneArray(*gammadposomarray);
        } else if (omegaprior == 3) {
            MasterSendGeneArray(*cauchydposomarray);
        }
    }

    void SlaveReceiveOmega() {
        if (omegaprior == 0) {
            SlaveReceiveGeneArray(*omegaarray);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmega((*omegaarray)[gene]);
            }
        } else if (omegaprior == 1) {
            SlaveReceiveGeneArray(*gammadposomarray);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmega((*gammadposomarray)[gene] + 1);
            }
        } else if (omegaprior == 2) {
            SlaveReceiveGeneArray(*gammadposomarray);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmega(exp((*gammadposomarray)[gene]));
            }
        } else if (omegaprior == 3) {
            SlaveReceiveGeneArray(*cauchydposomarray);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmega((*cauchydposomarray)[gene] + 1);
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
            double alpha = 1.0 / omegahyperinvshape;
            double beta = alpha / omegahypermean;
            omegaarray->SetShape(alpha);
            omegaarray->SetScale(beta);
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape);
            }
        } else {
            SlaveReceiveGlobal(dposompi);
            SlaveReceiveGlobal(dposomhypermean, dposomhyperinvshape);
            if ((omegaprior == 1) || (omegaprior == 2)) {
                double alpha = 1.0 / dposomhyperinvshape;
                double beta = alpha / dposomhypermean;
                gammadposomarray->SetPi(dposompi);
                gammadposomarray->SetShape(alpha);
                gammadposomarray->SetScale(beta);
            } else if (omegaprior == 3) {
                cauchydposomarray->SetPi(dposompi);
                cauchydposomarray->SetGamma(1.0 / dposomhyperinvshape);
            }
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetDPosOmHyperParameters(dposompi, dposomhypermean,
                                                            dposomhyperinvshape);
            }
        }
    }

    void MasterSendEpsilon()    {
        MasterSendGeneArray(*epsilonarray);
    }

    void SlaveReceiveEpsilon()  {
        SlaveReceiveGeneArray(*epsilonarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetEpsilon((*epsilonarray)[gene]);
        }
    }

    void SlaveSendEpsilon() {
        SlaveSendGeneArray(*epsilonarray);
    }

    void MasterReceiveEpsilon() {
        MasterReceiveGeneArray(*epsilonarray);
    }

    void MasterSendEpsilonHyperParameters() {
        MasterSendGlobal(epsilonhypermean, epsilonhyperinvconc);
    }

    void SlaveReceiveEpsilonHyperParameters()   {
        SlaveReceiveGlobal(epsilonhypermean, epsilonhyperinvconc);
    }

    void MasterSendPi() {
        MasterSendGeneArray(*piarray);
    }

    void SlaveReceivePi()   {
        SlaveReceiveGeneArray(*piarray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetPi((*piarray)[gene]);
        }
    }

    void SlaveSendPi()  {
        SlaveSendGeneArray(*piarray);
    }

    void MasterReceivePi()  {
        MasterReceiveGeneArray(*piarray);
    }

    void MasterSendPiHyperParameters()  {
        MasterSendGlobal(pihypermean, pihyperinvconc);
    }

    void SlaveReceivePiHyperParameters()    {
        SlaveReceiveGlobal(pihypermean, pihyperinvconc);
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
        MeanWidth = 0;
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
            MeanStatEnt += geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanAAEntropy();
            MeanWidth += geneprocess[gene]->GetNsite() * geneprocess[gene]->GetMeanWidth();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
        SlaveSendAdditive(MeanStatEnt);
        SlaveSendAdditive(MeanWidth);

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
        MeanWidth = 0;
        MasterReceiveAdditive(MeanStatEnt);
        MeanStatEnt /= GetTotNsite();
        MasterReceiveAdditive(MeanWidth);
        MeanWidth /= GetTotNsite();

        moveTime = 0;
        mapTime = 0;
        MasterReceiveAdditive(moveTime);
        MasterReceiveAdditive(mapTime);
        moveTime /= (GetNprocs() - 1);
        mapTime /= (GetNprocs() - 1);
    }
};
