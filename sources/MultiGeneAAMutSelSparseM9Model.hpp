#include "AAMutSelSparseM9Model.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"
#include "Permutation.hpp"

class MultiGeneAAMutSelSparseM9Model : public MultiGeneProbModel {
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

    vector<double> mixhyperparam;

    double& posmeanhypermean;
    double& posmeanhyperinvshape;
    IIDGamma* posmeanarray;
    GammaSuffStat posmeansuffstat;

    double& posinvshapehypermean;
    double& posinvshapehyperinvshape;
    IIDGamma* posinvshapearray;
    GammaSuffStat posinvshapesuffstat;

    double& poswhypermean;
    double& poswhyperinvconc;
    IIDBernoulliBeta* poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pospihypermean;
    double pospihyperinvconc;
    double& pospi;

    SimpleArray<double>* genednds;
    SimpleArray<double>* geneomegaarray;

    double epsilonhypermean;
    double epsilonhyperinvconc;
    IIDBeta* epsilonarray;
    BetaSuffStat epsilonhypersuffstat;

    double pihypermean;
    double pihyperinvconc;
    IIDBeta* piarray;
    BetaSuffStat pihypersuffstat;

    std::vector<AAMutSelSparseM9Model *> geneprocess;

    double lnL;
    double GeneLogPrior;
    double MeanStatEnt;
    double MeanWidth;
    double moveTime;
    double mapTime;
    Chrono movechrono;
    Chrono mapchrono;

    int blmode, nucmode, modalprior;

    Chrono totchrono;
    Chrono paramchrono;
    Chrono blchrono;

    int burnin;
    int chainsize;

  public:
    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneAAMutSelSparseM9Model(string indatafile, string intreefile, 
                                     int inblmode, int innucmode, int inmodalprior,
                                     double inepsilonhypermean, double inepsilonhyperinvconc,
                                     double inpihypermean, double inpihyperinvconc,
                                     double inpospihypermean, double inpospihyperinvconc,
                                     int inmyid, int innprocs)
        : MultiGeneProbModel(inmyid, innprocs),
        
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc),
          mixhyperparam(7, 0),
          posmeanhypermean(mixhyperparam[0]),
          posmeanhyperinvshape(mixhyperparam[1]),
          posinvshapehypermean(mixhyperparam[2]),
          posinvshapehyperinvshape(mixhyperparam[3]),
          poswhypermean(mixhyperparam[4]),
          poswhyperinvconc(mixhyperparam[5]),
          pospi(mixhyperparam[6])   {

        datafile = indatafile;
        treefile = intreefile;
        AllocateAlignments(datafile);

        burnin = 10;
        chainsize = 0;

        blmode = inblmode;
        nucmode = innucmode;
        modalprior = inmodalprior;

        epsilonhypermean = inepsilonhypermean;
        epsilonhyperinvconc = inepsilonhyperinvconc;
        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;

        pospihypermean = inpospihypermean;
        pospihyperinvconc = inpospihyperinvconc;
        pospi = pospihypermean;

        posmeanhypermean = 1.0;
        posmeanhyperinvshape = 1.0;
        posinvshapehypermean = 1.0;
        posinvshapehyperinvshape = 1.0;

        if (! pospi)   {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }
        else    {
            poswhypermean = 0.1;
            poswhyperinvconc = 1;
        }

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
            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene]->SetChainSize(insize);
            }
        }
    }

    void SetModalMixturePrior(int in)  {
        modalprior = in;
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

        double posmeanalpha = 1.0 / posmeanhyperinvshape;
        double posmeanbeta = posmeanalpha / posmeanhypermean;
        posmeanarray = new IIDGamma(GetLocalNgene(), posmeanalpha, posmeanbeta);

        double posinvshapealpha = 1.0 / posinvshapehyperinvshape;
        double posinvshapebeta = posinvshapealpha / posinvshapehypermean;
        posinvshapearray = new IIDGamma(GetLocalNgene(), posinvshapealpha, posinvshapebeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        poswarray = new IIDBernoulliBeta(GetLocalNgene(), pospi, poswalpha, poswbeta);
        poswarray->OffsetZeros(0.01);

        genednds = new SimpleArray<double>(GetLocalNgene());
        geneomegaarray = new SimpleArray<double>(GetLocalNgene());

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
            geneprocess.assign(0, (AAMutSelSparseM9Model *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (AAMutSelSparseM9Model *)0);

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
                    geneprocess[gene] = new AAMutSelSparseM9Model(
                        alivector[gene], tree, 3, 0, epsilon, pi, pospi);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    double epsilon = epsilonhyperinvconc ? -1.0 : epsilonhypermean;
                    double pi = pihyperinvconc ? -1.0 : pihypermean;
                    geneprocess[gene] = new AAMutSelSparseM9Model(GetLocalGeneName(gene),
                            treefile, 3, 0, epsilon, pi, pospi);
                }
            }

            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetBLMode(blmode);
                geneprocess[gene]->SetNucMode(nucmode);
                geneprocess[gene]->Allocate();
            }
        }
    }

    void SetMixtureHyperParameters(
            double inpospi, 
            double inposmeanhypermean, double inposmeanhyperinvshape,
            double inposinvshapehypermean, double inposinvshapehyperinvshape,
            double inposwhypermean, double inposwhyperinvconc)  {

        pospi = inpospi;

        posmeanhypermean = inposmeanhypermean;
        posmeanhyperinvshape = inposmeanhyperinvshape;
        posinvshapehypermean = inposinvshapehypermean;
        posinvshapehyperinvshape = inposinvshapehyperinvshape;

        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;

        if (!pospi) {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }
    }

    void SetMixtureArrays() {

        double posmeanalpha = 1.0 / posmeanhyperinvshape;
        double posmeanbeta = posmeanalpha / posmeanhypermean;
        posmeanarray->SetShape(posmeanalpha);
        posmeanarray->SetScale(posmeanbeta);

        double posinvshapealpha = 1.0 / posinvshapehyperinvshape;
        double posinvshapebeta = posinvshapealpha / posinvshapehypermean;
        posinvshapearray->SetShape(posinvshapealpha);
        posinvshapearray->SetScale(posinvshapebeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        poswarray->SetPi(pospi);
        poswarray->SetAlpha(poswalpha);
        poswarray->SetBeta(poswbeta);
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
        SetMixtureArrays();
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

            MasterSendMixtureHyperParameters();
            MasterSendMixture();
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

        SlaveReceiveMixtureHyperParameters();
        SlaveReceiveMixture();
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

            MasterSendMixtureHyperParameters();
            MasterSendMixture();
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

        SlaveReceiveMixtureHyperParameters();
        SlaveReceiveMixture();
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

    void MasterTraceSitePredictedDNDS(ostream &os) {
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
                        os << array[i++] << '\t';
                    }
                }
            }
            if (i != totnsite) {
                cerr << "error in MultiGeneCodonM2aModel::MasterTraceSiteOmega: non "
                        "matching number of sites\n";
                exit(1);
            }
            delete[] array;
        }
        os << '\n';
        os.flush();
    }

    void SlaveTraceSitePredictedDNDS() {
        int ngene = GetLocalNgene();
        int totnsite = GetLocalTotNsite();
        double *array = new double[totnsite];
        int i = 0;
        for (int gene = 0; gene < ngene; gene++) {
            geneprocess[gene]->GetSitePredictedDNDS(array + i);
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
            cerr << "error in SlaveTraceSitePredictedDNDS: non "
                    "matching number of sites\n";
            exit(1);
        }

        MPI_Send(array, totnsite, MPI_DOUBLE, 0, TAG1, MPI_COMM_WORLD);
        delete[] array;
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

    double GetMeanOmega() const {
        double m1 = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)  {
            m1 += geneomegaarray->GetVal(gene);
        }
        m1 /= GetLocalNgene();
        return m1;
    }

    double GetVarOmega() const {
        double m1 = 0;
        double m2 = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)  {
            double tmp = geneomegaarray->GetVal(gene);
            m1 += tmp;
            m2 += tmp*tmp;
        }
        m1 /= GetLocalNgene();
        m2 /= GetLocalNgene();
        m2 -= m1*m1;
        return m2;
    }

    // summary statistics for tracing MCMC
    int GetNpos() const {
        return GetNgene() - poswarray->GetNullSet(); 
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprior\tlnL\t";
        os << "length\t";
        if (blmode < 2) {
            os << "stdev\t";
        }
        os << "\tpospi";
        os << "\tnposfrac";
        os << "\tmeanom\tvarom";
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
        os << '\t' << pospi;
        os << '\t' << GetNpos();
        os << '\t' << GetMeanOmega() << '\t' << GetVarOmega();
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

    void TracePosWeight(ostream &os) const  {
        for (int gene = 0; gene < Ngene; gene++) {
            os << poswarray->GetVal(gene) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TraceGeneOmega(ostream &os) const  {
        for (int gene = 0; gene < Ngene; gene++) {
            os << geneomegaarray[gene] << '\t';
        }
        os << '\n';
        os.flush();
    }

    void MasterTraceSiteOmega(ostream &os) {
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
                        os << array[i++] << '\t';
                    }
                }
            }
            if (i != totnsite) {
                cerr << "error in MultiGeneCodonM2aModel::MasterTraceSiteOmega: non "
                        "matching number of sites\n";
                exit(1);
            }
            delete[] array;
        }
        os << '\n';
        os.flush();
    }

    void SlaveTraceSiteOmega() {
        int ngene = GetLocalNgene();
        int totnsite = GetLocalTotNsite();
        double *array = new double[totnsite];
        int i = 0;
        for (int gene = 0; gene < ngene; gene++) {
            geneprocess[gene]->GetSiteOmega(array + i);
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
            cerr << "error in MultiGeneCodonM2aModel::SlaveTraceSiteOmega: non "
                    "matching number of sites\n";
            exit(1);
        }

        MPI_Send(array, totnsite, MPI_DOUBLE, 0, TAG1, MPI_COMM_WORLD);
        delete[] array;
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

        is >> pospi;
        is >> posmeanhypermean >> posmeanhyperinvshape;
        is >> posinvshapehypermean >> posinvshapehyperinvshape;
        is >> poswhypermean >> poswhyperinvconc;

        is >> *posmeanarray;
        is >> *posinvshapearray;
        is >> *poswarray;

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

        os << pospi << '\t';
        os << posmeanhypermean << '\t' << posmeanhyperinvshape << '\t';
        os << posinvshapehypermean << '\t' << posinvshapehyperinvshape << '\t';
        os << poswhypermean << '\t' << poswhyperinvconc << '\t';

        os << *posmeanarray << '\t';
        os << *posinvshapearray << '\t';
        os << *poswarray << '\t';

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

        total += MixtureHyperLogPrior();

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

    // mixture
    double MixtureHyperLogPrior() const {
        double total = 0;
        total += PosPiHyperLogPrior();
        total += PosMeanHyperLogPrior();
        total += PosInvShapeHyperLogPrior();
        total += PosWHyperLogPrior();
        return total;
    }

    double PosPiHyperLogPrior() const  {
        double total = 0;
        if (pospi) {
            // beta distribution for pi, if not 0
            double pospialpha = pospihypermean / pospihyperinvconc;
            double pospibeta = (1 - pospihypermean) / pospihyperinvconc;
            total += (pospialpha - 1) * log(1.0 - pospi) + (pospibeta - 1) * log(pospi);
        }
        return total;
    }

    double PosMeanHyperLogPrior() const {
        double total = 0;
        total -= posmeanhypermean;
        total -= posmeanhyperinvshape;
        return total;
    }

    double PosInvShapeHyperLogPrior() const {
        double total = 0;
        total -= posinvshapehypermean;
        total -= posinvshapehyperinvshape;
        return total;
    }

    double PosWHyperLogPrior() const    {
        double total = - poswhyperinvconc;
        // distribution across genes should be modal
        double alpha = poswhypermean / poswhyperinvconc;
        double beta = (1 - poswhypermean) / poswhyperinvconc;
        if (modalprior && ((alpha < 1) || (beta < 1))) {
            total += log(0);
            // total += Random::INFPROB;
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

    // suff stat for gene-specific mixture parameters, as a function of mixture
    // hyperparameters
    double MixtureHyperSuffStatLogProb() const {
        double total = 0;

        double posmeanalpha = 1.0 / posmeanhyperinvshape;
        double posmeanbeta = posmeanalpha / posmeanhypermean;
        total += posmeansuffstat.GetLogProb(posmeanalpha, posmeanbeta);

        double posinvshapealpha = 1.0 / posinvshapehyperinvshape;
        double posinvshapebeta = posinvshapealpha / posinvshapehypermean;
        total += posinvshapesuffstat.GetLogProb(posinvshapealpha, posinvshapebeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        total += poswsuffstat.GetLogProb(pospi, poswalpha, poswbeta);

        return total;
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

    // log prob for moving mixture hyper params
    double MixtureHyperLogProb() const {
        return MixtureHyperLogPrior() + MixtureHyperSuffStatLogProb();
    }

    double EpsilonHyperLogProb() const { return EpsilonHyperLogPrior() + EpsilonHyperSuffStatLogProb();}

    double PiHyperLogProb() const { return PiHyperLogPrior() + PiHyperSuffStatLogProb();}

    //-------------------
    // Moves
    //-------------------

    void MasterMove() override {
        totchrono.Start();
        int nrep = 15;

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

            // mixture hyperparameters
            MasterReceiveMixtureHyperSuffStat();
            movechrono.Start();
            MoveMixtureHyperParameters();
            movechrono.Stop();
            MasterSendMixtureHyperParameters();

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

        MasterReceiveMixture();
        MasterReceiveGeneOmega();

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

        int nrep = 15;

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

            MoveGeneOmega();
            // mixture hyperparameters
            SlaveSendMixtureHyperSuffStat();
            SlaveReceiveMixtureHyperParameters();

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

        SlaveSendMixture();
        SlaveSendGeneOmega();

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

    void MoveGeneOmega() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveOmega();
            geneprocess[gene]->MoveOmegaHyperParameters();

            // collect new parameter values across genes (at the level of the slave)
            geneprocess[gene]->GetMixtureParameters(
                    (*poswarray)[gene],
                    (*posmeanarray)[gene], (*posinvshapearray)[gene]);

            (*geneomegaarray)[gene] = geneprocess[gene]->GetMeanOmega();
        }
    }

    void MoveGeneAA() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->MoveAA(2);
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
        ScalingMove(lambda, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::LambdaHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneAAMutSelSparseM9Model::LambdaHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    void MoveBranchLengthsHyperParameters() {
        for (int j = 0; j < Nbranch; j++) {
            BranchLengthsHyperScalingMove(1.0, 10);
            BranchLengthsHyperScalingMove(0.3, 10);
        }

        ScalingMove(blhyperinvshape, 1.0, 10,
                    &MultiGeneAAMutSelSparseM9Model::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10,
                    &MultiGeneAAMutSelSparseM9Model::BranchLengthsHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

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
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10,
                    &MultiGeneAAMutSelSparseM9Model::NucRatesHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

        UpdateNucRates();
    }

    void MoveEpsilonHyperParameters()   {
        epsilonhypersuffstat.Clear();
        epsilonarray->AddSuffStat(epsilonhypersuffstat);
        SlidingMove(epsilonhypermean, 1.0, 10, 0, 1.0, &MultiGeneAAMutSelSparseM9Model::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        SlidingMove(epsilonhypermean, 0.1, 10, 0, 1.0, &MultiGeneAAMutSelSparseM9Model::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(epsilonhyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(epsilonhyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::EpsilonHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

        double alpha = epsilonhypermean / epsilonhyperinvconc;
        double beta = (1.0 - epsilonhypermean) / epsilonhyperinvconc;
        epsilonarray->SetAlpha(alpha);
        epsilonarray->SetBeta(beta);
    }

    void MovePiHyperParameters()   {
        pihypersuffstat.Clear();
        piarray->AddSuffStat(pihypersuffstat);
        SlidingMove(pihypermean, 1.0, 10, 0, 1.0, &MultiGeneAAMutSelSparseM9Model::PiHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        SlidingMove(pihypermean, 0.1, 10, 0, 1.0, &MultiGeneAAMutSelSparseM9Model::PiHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(pihyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::PiHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(pihyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::PiHyperLogProb,
                &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

        double alpha = pihypermean / pihyperinvconc;
        double beta = (1.0 - pihypermean) / pihyperinvconc;
        piarray->SetAlpha(alpha);
        piarray->SetBeta(beta);
    }

    // moving mixture hyper params
    void MoveMixtureHyperParameters()   {

        ScalingMove(posmeanhypermean, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(posmeanhypermean, 0.3, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(posmeanhyperinvshape, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(posmeanhyperinvshape, 0.3, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

        ScalingMove(posinvshapehypermean, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(posinvshapehypermean, 0.3, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(posinvshapehyperinvshape, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
        ScalingMove(posinvshapehyperinvshape, 0.3, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                    &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

        MovePoswHyper();

        if (chainsize >= burnin) {
            if (pospihyperinvconc) {
                ResamplePosPi();
            }
        }
        SetMixtureArrays();
    }

    void MovePoswHyper()    {
            SlidingMove(poswhypermean, 1.0, 10, 0, 1, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                        &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
            SlidingMove(poswhypermean, 0.3, 10, 0, 1, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                        &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
            SlidingMove(poswhypermean, 0.1, 10, 0, 1, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                        &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 1.0, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                        &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 0.3, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                        &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 0.1, 10, &MultiGeneAAMutSelSparseM9Model::MixtureHyperLogProb,
                        &MultiGeneAAMutSelSparseM9Model::NoUpdate, this);

        for (int rep=0; rep<10; rep++)	{
            PoswCompMove(1.0);
            PoswCompMove(0.3);
            PoswCompMove(0.1);
        }
    }

    int PoswCompMove(double tuning) {
        double bkposwhypermean = poswhypermean;
        double bkposwhyperinvconc = poswhyperinvconc;
        double deltalogprob = - MixtureHyperLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        if (poswhypermean * e > 1.0)	{
            return 0;
        }
        poswhypermean *= e;
        poswhyperinvconc *= e;
        deltalogprob += MixtureHyperLogProb();
        deltalogprob += 2*m;

        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted) {
            poswhypermean = bkposwhypermean;
            poswhyperinvconc = bkposwhyperinvconc;
        }
        return accepted;
    }

    // special function for moving pi
    void ResamplePosPi()   {
        int n0 = poswsuffstat.GetN0();
        int n1 = poswsuffstat.GetN1();
        if ((n0 + n1) != Ngene) {
            cerr << "error in resample pi\n";
            exit(1);
        }
        double pospialpha = pospihypermean / pospihyperinvconc;
        double pospibeta = (1 - pospihypermean) / pospihyperinvconc;
        double postalpha = Random::sGamma(pospialpha + n1);
        double postbeta = Random::sGamma(pospibeta + n0);
        pospi = postalpha / (postalpha + postbeta);
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

    void SlaveSendMixture() {
        SlaveSendGeneArray(*posmeanarray, *posinvshapearray);
        SlaveSendGeneArray(*poswarray);
    }

    void MasterReceiveMixture() {
        MasterReceiveGeneArray(*posmeanarray, *posinvshapearray);
        MasterReceiveGeneArray(*poswarray);
    }

    void SlaveReceiveMixture()  {
        SlaveReceiveGeneArray(*posmeanarray, *posinvshapearray);
        SlaveReceiveGeneArray(*poswarray);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetMixtureParameters(
                    (*poswarray)[gene],
                    (*posmeanarray)[gene], (*posinvshapearray)[gene]);
        }
    }

    void MasterSendMixture()    {
        MasterSendGeneArray(*posmeanarray, *posinvshapearray);
        MasterSendGeneArray(*poswarray);
    }

    void MasterSendMixtureHyperParameters() {
        MasterSendGlobal(mixhyperparam); 
    }

    void SlaveReceiveMixtureHyperParameters()   {
        SlaveReceiveGlobal(mixhyperparam);
        SetMixtureArrays();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetMixtureHyperParameters(
                    pospi, 
                    posmeanhypermean, posmeanhyperinvshape,
                    posinvshapehypermean, posinvshapehyperinvshape,
                    poswhypermean, poswhyperinvconc);
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

    void SlaveSendMixtureHyperSuffStat()    {

        posmeansuffstat.Clear();
        posmeansuffstat.AddSuffStat(*posmeanarray);
        SlaveSendAdditive(posmeansuffstat);

        posinvshapesuffstat.Clear();
        posinvshapesuffstat.AddSuffStat(*posinvshapearray);
        SlaveSendAdditive(posinvshapesuffstat);

        poswsuffstat.Clear();
        poswarray->AddSuffStat(poswsuffstat);
        SlaveSendAdditive(poswsuffstat);
    }

    void MasterReceiveMixtureHyperSuffStat()    {

        posmeansuffstat.Clear();
        MasterReceiveAdditive(posmeansuffstat);

        posinvshapesuffstat.Clear();
        MasterReceiveAdditive(posinvshapesuffstat);

        poswsuffstat.Clear();
        MasterReceiveAdditive(poswsuffstat);
    }

    void SlaveSendGeneOmega() { 
        SlaveSendGeneArray(*geneomegaarray); 
    }

    void MasterReceiveGeneOmega() {
        MasterReceiveGeneArray(*geneomegaarray); 
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
