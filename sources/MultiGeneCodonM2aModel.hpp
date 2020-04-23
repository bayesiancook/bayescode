
/**
 * \brief A multi-gene version of CodonM2aModel
 *
 * Branch lengths and nucrates can be either:
 * - (2): shared across all genes
 * - (1): gene specific, with hyperparameters estimated across genes (shrinkage)
 * - (0): gene-specific, with fixed hyperparameters (no shrinkage); in that
 * case, the hyperparameters are set up so as to implement vague priors
 *
 * These three alternative modes can be tuned separately for branch lengths and
 * nuc rates by setting the variables blmode and nucmode to 0, 1 or 2. By
 * default, bl and nucrates are shared across genes.
 *
 * For each gene, the 3-component mixture of omega's across sites has four
 * parameters (see CodonM2aModel.hpp):
 * - 0 < purom < 1
 * - 0 < dposom < +infty
 * - 0 < purw < 1
 * - 0 <= posw < 1
 *
 * These 4 parameters are always gene-specific (they cannot be shared by all
 * genes); the priors over these parameters are as follows:
 * - purom ~ Beta(puromhypermean,puromhyperinvconc)
 * - dposom ~ Gamma(dposomhypermean,dposomhyperinvshape)
 * - purw ~ Beta purwhypermean,purwhyperinvconc)
 * - posw ~  mixture of point mass at 0 with prob 1-pi, and
 * Beta(poswhypermean,poswhyperinvconc) with prob pi
 *
 * In turn, the 9 hyperparameters of these priors can be either:
 * - (1) : estimated across genes
 * - (0): fixed (such as given by the user)
 *
 * again, separately for each set of hyperparameters, by setting the following
 * variables to 0 or 1: purommode, dposommode, purwmode, poswmode. By default,
 * the 9 mixture hyperparameters are estimated (shrunken) across genes (mode 1).
 * The hyperpriors over these parameters are then as follows:
 * - Uniform(0,1) for puromhypermean, purwhypermean and poswhypermean
 * - Exponential(1) for puromhyperinvconc, dposomhypermean, dposomhyperinvshape,
 * purwhyperinvconc, poswhyperinvconc
 * - Beta(pihypermean=0.1,pihyperinvconc=0.2) for pi
 *
 * Note that, when poswmode == 0, only poswhypermean and poswhyperinvconc are
 * fixed, on the other hand, pi is free. pi can be fixed by setting
 * pihyperinvconc to 0 (in which case pi is set equal to pihypermean).
 */

#include "Chrono.hpp"
#include "CodonM2aModel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "IIDDirichlet.hpp"
#include "MultiGeneMPIModule.hpp"
#include "MultiGeneProbModel.hpp"
#include "Parallel.hpp"

class MultiGeneCodonM2aModel : public MultiGeneProbModel {

    //-------------------
    // Data structures
    // ------------------

    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;
    std::vector<CodonSequenceAlignment*> alivector;

    string datapath;
    string datafile;
    string treefile;

    int Ntaxa;
    int Nbranch;

    double lambda;
    BranchIIDGamma *branchlength;
    GammaSuffStat hyperlengthsuffstat;

    double blhyperinvshape;
    GammaWhiteNoiseArray *branchlengtharray;
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;
    GammaSuffStatBranchArray *lengthhypersuffstatarray;

    PoissonSuffStatTreeArray* lengthpathsuffstattreearray;

    vector<double> mixhyperparam;

    double &puromhypermean;
    double &puromhyperinvconc;
    IIDBeta *puromarray;
    BetaSuffStat puromsuffstat;

    double &dposomhypermean;
    double &dposomhyperinvshape;
    IIDGamma *dposomarray;
    GammaSuffStat dposomsuffstat;

    double &purwhypermean;
    double &purwhyperinvconc;
    IIDBeta *purwarray;
    BetaSuffStat purwsuffstat;

    double &poswhypermean;
    double &poswhyperinvconc;
    IIDBernoulliBeta *poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double &pi;

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

    std::vector<CodonM2aModel *> geneprocess;

    double lnL;
    double GeneLogPrior;
    double GeneBLLogPrior;
    double GeneNucRatesLogPrior;
    double GeneOmegaLogPrior;
    double moveTime;
    double mapTime;

    int burnin;

    // 0: free (fixed hyper parameters)
    // 1: free and shrinkage (free hyper parameters)
    // 2: shared across genes
    int blmode;
    int blsamplemode;
    int nucmode;
    int purommode;
    int dposommode;
    int purwmode;
    int poswmode;

    int modalprior;

    Chrono movechrono;
    Chrono mapchrono;
    Chrono totchrono;


  public:
    //-------------------
    // Constructors
    // ------------------

    MultiGeneCodonM2aModel(string indatapath, string indatafile, string intreefile, double inpihypermean,
                           double inpihyperinvconc, int inmyid, int innprocs) :
          MultiGeneProbModel(inmyid, innprocs),
          mixhyperparam(9, 0),
          puromhypermean(mixhyperparam[0]),
          puromhyperinvconc(mixhyperparam[1]),
          dposomhypermean(mixhyperparam[2]),
          dposomhyperinvshape(mixhyperparam[3]),
          purwhypermean(mixhyperparam[4]),
          purwhyperinvconc(mixhyperparam[5]),
          poswhypermean(mixhyperparam[6]),
          poswhyperinvconc(mixhyperparam[7]),
          pi(mixhyperparam[8]),
          nucrelratesuffstat(Nrr),
          nucstatsuffstat(Nnuc) {

        burnin = 0;

        // 0 : gathering branch lengths across genes and then using gamma suff stats to resample hyperparams
        // 1 : gathering suff stats across genes, then resampling hyperparams based on integrated bls
        blsamplemode = 0;

        modalprior = 1;

        pihypermean = inpihypermean;
        pihyperinvconc = inpihyperinvconc;
        pi = pihypermean;

        datapath = indatapath;
        datafile = indatafile;
        treefile = intreefile;
        AllocateAlignments(datafile, datapath);

        // all datafiles have all taxa (with missing data if needed) in same order
        // makes it easier to register tree with data, etc.

        refcodondata = new CodonSequenceAlignment(refdata, true);
        taxonset = refdata->GetTaxonSet();
        Ntaxa = refdata->GetNtaxa();

        // get tree from file (newick format)
        tree = new Tree(datapath + treefile);

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
        } else {
            branchlength->SetAllBranches(1.0 / lambda);
            branchlengtharray =
                new GammaWhiteNoiseArray(GetLocalNgene(), *tree, *branchlength, 1.0 / blhyperinvshape);

            if (blsamplemode == 1)   {
                lengthpathsuffstatarray = 0;
                lengthpathsuffstattreearray = new PoissonSuffStatTreeArray(*tree,GetLocalNgene());
                lengthhypersuffstatarray = 0;
            }
            else    {
                lengthpathsuffstatarray = 0;
                // lengthpathsuffstattreearray = 0;
                if (myid)   {
                    lengthpathsuffstattreearray = new PoissonSuffStatTreeArray(*tree,GetLocalNgene());
                }
                lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
            }
        }

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

        double puromalpha = puromhypermean / puromhyperinvconc;
        double purombeta = (1 - puromhypermean) / puromhyperinvconc;
        puromarray = new IIDBeta(GetLocalNgene(), puromalpha, purombeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray = new IIDGamma(GetLocalNgene(), dposomalpha, dposombeta);

        double purwalpha = purwhypermean / purwhyperinvconc;
        double purwbeta = (1 - purwhypermean) / purwhyperinvconc;
        purwarray = new IIDBeta(GetLocalNgene(), purwalpha, purwbeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        poswarray = new IIDBernoulliBeta(GetLocalNgene(), pi, poswalpha, poswbeta);

        if (!GetMyid()) {
            geneprocess.assign(0, (CodonM2aModel *)0);
        } else {
            geneprocess.assign(GetLocalNgene(), (CodonM2aModel *)0);

            ifstream is((datapath + datafile).c_str());
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
                    geneprocess[gene] = new CodonM2aModel(alivector[gene], tree, pi);
                }
            }
            else    {
                for (int gene = 0; gene < GetLocalNgene(); gene++) {
                    geneprocess[gene] = new CodonM2aModel(datapath, GetLocalGeneName(gene), treefile, pi);
                }
            }
            for (int gene = 0; gene < GetLocalNgene(); gene++) {
                geneprocess[gene]->SetAcrossGenesModes(blmode, nucmode);
                geneprocess[gene]->Allocate();
            }
        }
    }

    //-------------------
    // Accessors
    // ------------------

    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)refcodondata->GetStateSpace();
    }

    //-------------------
    // Setting and updating
    // ------------------

    // called upon constructing the model
    // mode == 2: global (only for branch lengths and nuc rates)
    // mode == 1: gene specific, with hyperparameters estimated across genes
    // mode == 0: gene-specific, with fixed hyperparameters
    void SetAcrossGenesModes(int inblmode, int innucmode, int inpurommode, int indposommode,
                             int inpurwmode, int inposwmode)    {
        blmode = inblmode;
        nucmode = innucmode;
        purommode = inpurommode;
        dposommode = indposommode;
        purwmode = inpurwmode;
        poswmode = inposwmode;
    }

    void SetBLSamplingMode(int inmode)  {
        blsamplemode = inmode;
    }

    void SetModalMixturePrior(int in)  {
        modalprior = in;
    }

    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc,
                                   double indposomhypermean, double indposomhyperinvshape,
                                   double inpurwhypermean, double inpurwhyperinvconc,
                                   double inposwhypermean, double inposwhyperinvconc)   {
        puromhypermean = inpuromhypermean;
        puromhyperinvconc = inpuromhyperinvconc;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
        purwhypermean = inpurwhypermean;
        purwhyperinvconc = inpurwhyperinvconc;
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;
    }

    void SetBurnin(double in)   {
        burnin = in;
    }

    void UpdateNucMatrix()  {
        nucmatrix->CopyStationary((*nucstatarray)[0]);
        nucmatrix->CorruptMatrix();
    }

    void SetMixtureArrays() {
        double puromalpha = puromhypermean / puromhyperinvconc;
        double purombeta = (1 - puromhypermean) / puromhyperinvconc;
        puromarray->SetAlpha(puromalpha);
        puromarray->SetBeta(purombeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        dposomarray->SetShape(dposomalpha);
        dposomarray->SetScale(dposombeta);

        double purwalpha = purwhypermean / purwhyperinvconc;
        double purwbeta = (1 - purwhypermean) / purwhyperinvconc;
        purwarray->SetAlpha(purwalpha);
        purwarray->SetBeta(purwbeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        poswarray->SetPi(pi);
        poswarray->SetAlpha(poswalpha);
        poswarray->SetBeta(poswbeta);
    }

    void FastUpdate()   {
        branchlength->SetScale(lambda);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }
        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
        SetMixtureArrays();
    }

    void MasterUpdate() override    {
        FastUpdate();
        if (nprocs > 1) {
            MasterSendBranchLengthsHyperParameters();
            MasterSendNucRatesHyperParameters();
            MasterSendMixtureHyperParameters();

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

            MasterSendMixture();
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {
        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveMixtureHyperParameters();
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

        SlaveReceiveMixture();
        GeneUpdate();
        SlaveSendLogProbs();
    }

    void GeneUpdate()   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->Update();
        }
    }

    void NoUpdate() {}

    void MasterPostPred(string name) override   {
        FastUpdate();
        if (nprocs > 1) {
            MasterSendBranchLengthsHyperParameters();
            MasterSendNucRatesHyperParameters();
            MasterSendMixtureHyperParameters();

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

            MasterSendMixture();
            // MasterReceiveLogProbs();
        }

        ofstream os((name + ".ppredparams").c_str());
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            os << GetLocalGeneName(gene) << '\t' << 0 << '\t' << poswarray->GetVal(gene) << '\t' << 1.0 + dposomarray->GetVal(gene) << '\n';
        }
    }

    void SlavePostPred(string name) override    {
        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveMixtureHyperParameters();
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

        SlaveReceiveMixture();
        GenePostPred(name);
        // SlaveSendLogProbs();
    }

    void GenePostPred(string name)  {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->PostPred(name + GetLocalGeneName(gene));
        }
    }


    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream &os) const override    {
        os << "#logprior\tlnL";
        if (blmode == 2) {
            os << "\tlength";
        } else {
            os << "\tmeanlength\tstdev";
        }
        os << "\tpi";
        os << "\tnposfrac";
        os << "\tpurommean\tinvconc";
        os << "\tdposommean\tinvshape";
        os << "\tpurwmean\tinvconc";
        os << "\tposwmean\tinvconc";
        os << "\tstatent";
        os << "\trrent";
        if (nucmode != 2) {
            os << "\tstdevrr\tcenter\thyperinvconc";
            os << "\tstdevstat\tcenter\thyperinvconc";
        }
        os << "\tgenelogprior\tbl\tnuc\tomega";
        os << "\tblinvshape\tbllogprob";
        os << '\n';
    }

    void Trace(ostream &os) const override  {
        os << GetLogPrior();
        os << '\t' << GetLogLikelihood();
        if (blmode == 2) {
            os << '\t' << GetMeanTotalLength();
        } else {
            os << '\t' << GetMeanLength();
            os << '\t' << sqrt(GetVarLength());
        }
        os << '\t' << pi;
        os << '\t' << GetNpos();
        os << '\t' << puromhypermean << '\t' << puromhyperinvconc;
        os << '\t' << dposomhypermean << '\t' << dposomhyperinvshape;
        os << '\t' << purwhypermean << '\t' << purwhyperinvconc;
        os << '\t' << poswhypermean << '\t' << poswhyperinvconc;
        os << '\t' << nucstatarray->GetMeanEntropy();
        os << '\t' << nucrelratearray->GetMeanEntropy();
        if (nucmode != 2) {
            os << '\t' << sqrt(GetVarNucRelRate()) << '\t' << Random::GetEntropy(nucrelratehypercenter)
               << '\t' << nucrelratehyperinvconc;
            os << '\t' << sqrt(GetVarNucStat()) << '\t' << Random::GetEntropy(nucstathypercenter)
               << '\t' << nucstathyperinvconc;
        }
        os << '\t' << GeneLogPrior << '\t' << GeneBLLogPrior << '\t' << GeneNucRatesLogPrior << '\t'
           << GeneOmegaLogPrior;
        os << '\t' << blhyperinvshape << '\t' << branchlength->GetLogProb();
        os << '\n';
        os.flush();
    }

    void TracePosWeight(ostream &os) const  {
        for (int gene = 0; gene < Ngene; gene++) {
            os << poswarray->GetVal(gene) << '\t';
        }
        os << '\n';
        os.flush();
    }

    void TracePosOm(ostream &os) const  {
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

    void MasterTraceSitesPostProb(ostream &os)  {
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

    void SlaveTraceSitesPostProb()  {
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

    void Monitor(ostream &os) const override {}

    void MasterFromStream(istream &is) override {

        is >> puromhypermean >> puromhyperinvconc;
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> purwhypermean >> purwhyperinvconc;
        is >> poswhypermean >> poswhyperinvconc;
        is >> pi;
        is >> *puromarray;
        is >> *dposomarray;
        is >> *purwarray;
        is >> *poswarray;

        is >> nucrelratehypercenter;
        is >> nucrelratehyperinvconc;
        is >> nucstathypercenter;
        is >> nucstathyperinvconc;
        is >> *nucrelratearray;
        is >> *nucstatarray;

        if (blmode == 2) {
            is >> lambda;
            is >> *branchlength;
        } else {
            is >> lambda;
            is >> *branchlength;
            is >> blhyperinvshape;
            is >> *branchlengtharray;
        }
    }

    void MasterToStream(ostream &os) const override {

        os << puromhypermean << '\t' << puromhyperinvconc << '\t';
        os << dposomhypermean << '\t' << dposomhyperinvshape << '\t';
        os << purwhypermean << '\t' << purwhyperinvconc << '\t';
        os << poswhypermean << '\t' << poswhyperinvconc << '\t';
        os << pi << '\t';
        os << *puromarray << '\t';
        os << *dposomarray << '\t';
        os << *purwarray << '\t';
        os << *poswarray << '\t';

        os << nucrelratehypercenter << '\t';
        os << nucrelratehyperinvconc << '\t';
        os << nucstathypercenter << '\t';
        os << nucstathyperinvconc << '\t';
        os << *nucrelratearray << '\t';
        os << *nucstatarray << '\t';

        if (blmode == 2) {
            os << lambda << '\t';
            os << *branchlength << '\n';
        } else {
            os << lambda << '\t';
            os << *branchlength << '\t';
            os << blhyperinvshape << '\t';
            os << *branchlengtharray << '\n';
        }
    }


    // summary statistics for tracing MCMC
    int GetNpos() const {
        return GetNgene() - poswarray->GetNullSet(); 
    }

    double GetMeanTotalLength() const   {
        double tot = 0;
        for (int j = 0; j < Nbranch; j++) {
            tot += branchlength->GetVal(j);
        }
        return tot;
    }


    double GetMeanLength() const    {
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

    double GetVarNucStat() const    {
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

    //-------------------
    // Log Probs
    // ------------------

    //-------------------
    // Log Priors
    // ------------------

    // total log prior
    double GetLogPrior() const  {
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
            total += GlobalNucRatesLogPrior();
        } else if (nucmode == 1) {
            total += GeneNucRatesHyperLogPrior();
        } else {
            // nothing: everything accounted for by gene component
        }

        // mixture
        total += MixtureHyperLogPrior();

        return total;
    }

    // branch lengths
    // exponential of mean 10 for lambda
    double LambdaHyperLogPrior() const { return -lambda / 10; }

    double GlobalBranchLengthsLogPrior() const {
        return LambdaHyperLogPrior() + branchlength->GetLogProb();
    }

    // exponential of mean 1 for blhyperinvshape
    double BranchLengthsHyperInvShapeLogPrior() const { return -blhyperinvshape; }

    double GeneBranchLengthsHyperLogPrior() const {
        return BranchLengthsHyperInvShapeLogPrior() + branchlength->GetLogProb();
    }

    // nuc rates
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

    // mixture
    double MixtureHyperLogPrior() const {
        double total = 0;
        if (pi) {
            // beta distribution for pi, if not 0
            double pialpha = pihypermean / pihyperinvconc;
            double pibeta = (1 - pihypermean) / pihyperinvconc;
            total += (pialpha - 1) * log(1.0 - pi) + (pibeta - 1) * log(pi);
        }
        // exponential of mean 1 for purom and purw hyperinvconc
        total -= puromhyperinvconc;
        total -= purwhyperinvconc;

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

    //-------------------
    // Log Likelihood
    // ------------------

    double GetLogLikelihood() const { return lnL; }

    //-------------------
    // Suff Stat Log Probs
    // ------------------

    // suff stat for global branch lengths, as a function of lambda
    double LambdaHyperSuffStatLogProb() const {
        return hyperlengthsuffstat.GetLogProb(1.0, lambda);
    }

    // suff stat for gene-specific branch lengths, as a function of bl
    // hyperparameters
    double BranchLengthsHyperSuffStatLogProb() const {
        return lengthhypersuffstatarray->GetLogProb(*branchlength, blhyperinvshape);
    }

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

    // suff stat for gene-specific mixture parameters, as a function of mixture
    // hyperparameters
    double MixtureHyperSuffStatLogProb() const {
        double total = 0;

        double puromalpha = puromhypermean / puromhyperinvconc;
        double purombeta = (1 - puromhypermean) / puromhyperinvconc;
        total += puromsuffstat.GetLogProb(puromalpha, purombeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        total += dposomsuffstat.GetLogProb(dposomalpha, dposombeta);

        double purwalpha = purwhypermean / purwhyperinvconc;
        double purwbeta = (1 - purwhypermean) / purwhyperinvconc;
        total += purwsuffstat.GetLogProb(purwalpha, purwbeta);

        double poswalpha = poswhypermean / poswhyperinvconc;
        double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
        total += poswsuffstat.GetLogProb(pi, poswalpha, poswbeta);
        return total;
    }

    //-------------------
    // Log Probs for specific MH Moves
    // ------------------

    // logprob for moving lambda
    double LambdaHyperLogProb() const {
        return LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
    }

    // logprob for moving hyperparameters of gene-specific branchlengths
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperInvShapeLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // logprob for moving branch-length hyperparameters, integrated over gene-specific branchlengths
    double BranchLengthsIntegratedHyperLogProb() const {
        return BranchLengthsHyperInvShapeLogPrior() + GeneBranchLengthsHyperSuffStatLogProb();
    }

    // log prob of substitution mappings across genes and branches analytically integrated over branch lengths
    double GeneBranchLengthsHyperSuffStatLogProb() const    {
        double total = 0;
        for (int j=0; j<Nbranch; j++)   {
            total += GeneBranchLengthsHyperSuffStatLogProb(j);
        }
        return total;
    }

    // log prob of substitution mappings across genes for a given branch analytically integrated over branch length
    double GeneBranchLengthsHyperSuffStatLogProb(int branch) const  {
        double total = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)    {
            total += GeneBranchLengthsHyperSuffStatLogProb(gene,branch);
        }
        return total;
    }

    // log prob of substitution mapping for a given gene and over a given branch analytically integrated over branch length
    double GeneBranchLengthsHyperSuffStatLogProb(int gene, int branch) const    {
        const PoissonSuffStat &suffstat = lengthpathsuffstattreearray->GetVal(gene).GetVal(branch);
        int count = suffstat.GetCount();
        double b = suffstat.GetBeta();

        double alpha = 1.0 / blhyperinvshape;
        double beta = alpha / branchlength->GetVal(branch);

        return alpha * log(beta) - Random::logGamma(alpha) + Random::logGamma(alpha + count) - (alpha + count) * log(beta + b);
    }

    // log prob for moving mixture hyper params
    double MixtureHyperLogProb() const {
        return MixtureHyperLogPrior() + MixtureHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates hyper params
    double NucRatesHyperLogProb() const {
        return GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates
    double NucRatesLogProb() const { return GlobalNucRatesLogPrior() + NucRatesSuffStatLogProb(); }

    //-------------------
    // Moves
    // ------------------

    // general move schedule
    void MasterMove() override  {
        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            // mixture hyperparameters
            MasterReceiveMixtureHyperSuffStat();
            movechrono.Start();
            MoveMixtureHyperParameters();
            movechrono.Stop();
            MasterSendMixtureHyperParameters();

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                MasterReceiveBranchLengthsSuffStat();
                movechrono.Start();
                ResampleBranchLengths();
                MoveLambda();
                movechrono.Stop();
                MasterSendGlobalBranchLengths();
            } else if (blmode == 1) {
                if (blsamplemode == 1)   {
                    MasterReceiveGeneBranchLengthsSuffStat();
                    movechrono.Start();
                    MoveBranchLengthsHyperParametersIntegrated();
                    movechrono.Stop();
                    MasterSendBranchLengthsHyperParameters();
                    // MasterReceiveGeneBranchLengths();
                }
                else    {
                    MasterReceiveBranchLengthsHyperSuffStat();
                    movechrono.Start();
                    MoveBranchLengthsHyperParameters();
                    movechrono.Stop();
                    MasterSendBranchLengthsHyperParameters();
                }
            }

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                MasterReceiveNucPathSuffStat();
                movechrono.Start();
                MoveNucRates();
                movechrono.Stop();
                MasterSendGlobalNucRates();
            } else if (nucmode == 1) {
                MasterReceiveNucRatesHyperSuffStat();
                movechrono.Start();
                MoveNucRatesHyperParameters();
                movechrono.Stop();
                MasterSendNucRatesHyperParameters();
            }
        }
        burnin++;
        if (blmode != 2) {
            MasterReceiveGeneBranchLengths();
        }
        if (nucmode != 2) {
            MasterReceiveGeneNucRates();
        }
        MasterReceiveMixture();
        MasterReceiveLogProbs();
    }

    void SlaveMove() override   {
        movechrono.Start();
        mapchrono.Start();
        GeneResampleSub(1.0);
        mapchrono.Stop();
        movechrono.Stop();

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            // gene specific mixture parameters
            // possibly branch lengths and nuc rates (if mode == 1 or 2)
            movechrono.Start();
            MoveGeneParameters(1.0);
            movechrono.Stop();

            // mixture hyperparameters
            SlaveSendMixtureHyperSuffStat();
            SlaveReceiveMixtureHyperParameters();
            ResampleEmptyOmegas();

            // global branch lengths, or gene branch lengths hyperparameters
            if (blmode == 2) {
                SlaveSendBranchLengthsSuffStat();
                SlaveReceiveGlobalBranchLengths();
            } else if (blmode == 1) {
                if (blsamplemode == 1)   {
                    // integrated move
                    SlaveSendGeneBranchLengthsSuffStat();
                    SlaveReceiveBranchLengthsHyperParameters();
                    ResampleGeneBranchLengths();
                }
                else    {
                    // conditional move
                    SlaveSendBranchLengthsHyperSuffStat();
                    SlaveReceiveBranchLengthsHyperParameters();
                    GeneResampleEmptyBranches();
                }
            }

            // global nucrates, or gene nucrates hyperparameters
            if (nucmode == 2) {
                SlaveSendNucPathSuffStat();
                SlaveReceiveGlobalNucRates();
            } else if (nucmode == 1) {
                SlaveSendNucRatesHyperSuffStat();
                SlaveReceiveNucRatesHyperParameters();
            }
        }
        burnin++;

        // collect current state
        if (blmode != 2) {
            SlaveSendGeneBranchLengths();
        }
        if (nucmode != 2) {
            SlaveSendGeneNucRates();
        }
        SlaveSendMixture();
        SlaveSendLogProbs();
    }

    // moving gene-specific parameters
    void GeneResampleSub(double frac)   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void MoveGeneParameters(int nrep)   {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            // move gene-specific parameters
            geneprocess[gene]->MoveParameters(nrep);

            // collect new parameter values across genes (at the level of the slave)
            geneprocess[gene]->GetMixtureParameters((*puromarray)[gene], (*dposomarray)[gene],
                                                    (*purwarray)[gene], (*poswarray)[gene]);
            if (blmode != 2) {
                geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
            }
            if (nucmode != 2) {
                geneprocess[gene]->GetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
            }
        }
    }

    // moving global branch lengths and lambda hyperparam
    void ResampleBranchLengths()    {
        branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

    void MoveLambda()   {
        hyperlengthsuffstat.Clear();
        hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda, 1.0, 10, &MultiGeneCodonM2aModel::LambdaHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(lambda, 0.3, 10, &MultiGeneCodonM2aModel::LambdaHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        branchlength->SetScale(lambda);
    }

    // resampling gene-specific branch lengths
    void ResampleGeneBranchLengths()    {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleBranchLengths();
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
    }

    // moving hyperparams for gene-specific branch lengths
    void MoveBranchLengthsHyperParameters() {
        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);

        ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneCodonM2aModel::BranchLengthsHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneCodonM2aModel::BranchLengthsHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);

        branchlengtharray->SetShape(1.0 / blhyperinvshape);
        MoveLambda();
    }

    double BranchLengthsHyperScalingMove(double tuning, int nrep)   {
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

    // moving hyperparams for gene-specific branch lengths, integrated version
    void MoveBranchLengthsHyperParametersIntegrated()   {
        BranchLengthsHyperScalingMoveIntegrated(1.0, 10);
        BranchLengthsHyperScalingMoveIntegrated(0.3, 10);

        ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneCodonM2aModel::BranchLengthsIntegratedHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneCodonM2aModel::BranchLengthsIntegratedHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);

        branchlengtharray->SetShape(1.0 / blhyperinvshape);
        MoveLambda();
    }

    double BranchLengthsHyperScalingMoveIntegrated(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int j = 0; j < Nbranch; j++) {
                double deltalogprob =
                    -branchlength->GetLogProb(j) - 
                    GeneBranchLengthsHyperSuffStatLogProb(j);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*branchlength)[j] *= e;
                deltalogprob +=
                    branchlength->GetLogProb(j) +
                    GeneBranchLengthsHyperSuffStatLogProb(j);
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

    // moving global nuc rates
    void MoveNucRates() {
        vector<double> &nucrelrate = (*nucrelratearray)[0];
        ProfileMove(nucrelrate, 0.1, 1, 10, &MultiGeneCodonM2aModel::NucRatesLogProb,
                    &MultiGeneCodonM2aModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.03, 3, 10, &MultiGeneCodonM2aModel::NucRatesLogProb,
                    &MultiGeneCodonM2aModel::UpdateNucMatrix, this);
        ProfileMove(nucrelrate, 0.01, 3, 10, &MultiGeneCodonM2aModel::NucRatesLogProb,
                    &MultiGeneCodonM2aModel::UpdateNucMatrix, this);

        vector<double> &nucstat = (*nucstatarray)[0];
        ProfileMove(nucstat, 0.1, 1, 10, &MultiGeneCodonM2aModel::NucRatesLogProb,
                    &MultiGeneCodonM2aModel::UpdateNucMatrix, this);
        ProfileMove(nucstat, 0.01, 1, 10, &MultiGeneCodonM2aModel::NucRatesLogProb,
                    &MultiGeneCodonM2aModel::UpdateNucMatrix, this);
    }

    // moving gene-specific nuc rates hyper params
    void MoveNucRatesHyperParameters()  {
        ProfileMove(nucrelratehypercenter, 1.0, 1, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.3, 1, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ProfileMove(nucrelratehypercenter, 0.1, 3, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 1.0, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.3, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(nucrelratehyperinvconc, 0.03, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);

        ProfileMove(nucstathypercenter, 1.0, 1, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.3, 1, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ProfileMove(nucstathypercenter, 0.1, 2, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 1.0, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.3, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);
        ScalingMove(nucstathyperinvconc, 0.03, 10, &MultiGeneCodonM2aModel::NucRatesHyperLogProb,
                    &MultiGeneCodonM2aModel::NoUpdate, this);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    }

    // moving mixture hyper params
    void MoveMixtureHyperParameters()   {
        if (purommode == 1) {
            SlidingMove(puromhypermean, 1.0, 10, 0, 1, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            SlidingMove(puromhypermean, 0.3, 10, 0, 1, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(puromhyperinvconc, 1.0, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(puromhyperinvconc, 0.3, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
        }

        if (dposommode == 1) {
            ScalingMove(dposomhypermean, 1.0, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(dposomhypermean, 0.3, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(dposomhyperinvshape, 1.0, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(dposomhyperinvshape, 0.3, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
        }

        if (poswmode == 1) {
            MovePoswHyper();
        }

        if (purwmode == 1) {
            SlidingMove(purwhypermean, 1.0, 10, 0, 1, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            SlidingMove(purwhypermean, 0.3, 10, 0, 1, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(purwhyperinvconc, 1.0, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(purwhyperinvconc, 0.3, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
        }

        if (burnin > 10) {
            if (pihyperinvconc) {
                ResamplePi();
            }
        }
        SetMixtureArrays();
    }

    void MovePoswHyper()    {
            SlidingMove(poswhypermean, 1.0, 10, 0, 1, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            SlidingMove(poswhypermean, 0.3, 10, 0, 1, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            SlidingMove(poswhypermean, 0.1, 10, 0, 1, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 1.0, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 0.3, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);
            ScalingMove(poswhyperinvconc, 0.1, 10, &MultiGeneCodonM2aModel::MixtureHyperLogProb,
                        &MultiGeneCodonM2aModel::NoUpdate, this);

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
    void ResamplePi()   {
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
    // MPI send/receive
    // ------------------

    // global parameters

    void MasterSendGlobalBranchLengths()    {
        MasterSendGlobal(*branchlength); 
    }

    void SlaveReceiveGlobalBranchLengths()  {
        SlaveReceiveGlobal(*branchlength);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengths(*branchlength);
        }
    }

    void MasterSendGlobalNucRates() {
        MasterSendGlobal(nucrelratearray->GetVal(0), nucstatarray->GetVal(0));
    }

    void SlaveReceiveGlobalNucRates()   {
        SlaveReceiveGlobal((*nucrelratearray)[0], (*nucstatarray)[0]);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[0], (*nucstatarray)[0]);
        }
    }

    // gene-specific parameters

    void MasterSendGeneBranchLengths()  {
        MasterSendGeneArray(*branchlengtharray);
    }

    void SlaveReceiveGeneBranchLengths()    {
        SlaveReceiveGeneArray(*branchlengtharray);
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengths(branchlengtharray->GetVal(gene));
        }
    }

    void SlaveSendGeneBranchLengths()   {
        SlaveSendGeneArray(*branchlengtharray);
    }

    void MasterReceiveGeneBranchLengths()   {
        MasterReceiveGeneArray(*branchlengtharray);
    }

    void MasterSendGeneNucRates()   {
        MasterSendGeneArray(*nucrelratearray, *nucstatarray);
    }

    void SlaveReceiveGeneNucRates() {
        SlaveReceiveGeneArray(*nucrelratearray, *nucstatarray);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
        }
    }

    void SlaveSendGeneNucRates()    {
        SlaveSendGeneArray(*nucrelratearray, *nucstatarray);
    }

    void MasterReceiveGeneNucRates()    {
        MasterReceiveGeneArray(*nucrelratearray, *nucstatarray);
    }

    void SlaveSendMixture() {
        SlaveSendGeneArray(*puromarray, *dposomarray);
        SlaveSendGeneArray(*purwarray, *poswarray);
    }

    void MasterReceiveMixture() {
        MasterReceiveGeneArray(*puromarray, *dposomarray);
        MasterReceiveGeneArray(*purwarray, *poswarray);
    }

    void SlaveReceiveMixture()  {
        SlaveReceiveGeneArray(*puromarray, *dposomarray);
        SlaveReceiveGeneArray(*purwarray, *poswarray);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetMixtureParameters((*puromarray)[gene], (*dposomarray)[gene],
                                                    (*purwarray)[gene], (*poswarray)[gene]);
        }
    }

    void MasterSendMixture()    {
        MasterSendGeneArray(*puromarray, *dposomarray);
        MasterSendGeneArray(*purwarray, *poswarray);
    }

    // hyper parameters

    void MasterSendBranchLengthsHyperParameters()   {
        MasterSendGlobal(*branchlength, blhyperinvshape);
    }

    void SlaveReceiveBranchLengthsHyperParameters() {
        SlaveReceiveGlobal(*branchlength, blhyperinvshape);
        if (blmode == 1) {
            branchlengtharray->SetShape(1.0 / blhyperinvshape);
        }
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
        }
    }

    void GeneResampleEmptyBranches()    {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->ResampleEmptyBranches();
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
    }

    void MasterSendNucRatesHyperParameters()    {
        MasterSendGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        MasterSendGlobal(nucstathypercenter, nucstathyperinvconc);
    }

    void SlaveReceiveNucRatesHyperParameters()  {
        SlaveReceiveGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
        SlaveReceiveGlobal(nucstathypercenter, nucstathyperinvconc);

        nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
        nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter, nucrelratehyperinvconc,
                                                          nucstathypercenter, nucstathyperinvconc);
        }
    }

    void MasterSendMixtureHyperParameters() {
        MasterSendGlobal(mixhyperparam); 
    }

    void SlaveReceiveMixtureHyperParameters()   {
        SlaveReceiveGlobal(mixhyperparam);
        SetMixtureArrays();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetMixtureHyperParameters(
                puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape, pi,
                purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc);
        }
    }

    void ResampleEmptyOmegas()  {
        // resample dposom from the prior, for those genes for which posw == 0
        dposomarray->PriorResample(*poswarray);
        // notify the gene processes of this change
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetMixtureParameters((*puromarray)[gene], (*dposomarray)[gene],
                                                    (*purwarray)[gene], (*poswarray)[gene]);
        }
    }

    // suff stats

    // branch lengths
    void SlaveSendBranchLengthsSuffStat()   {
        lengthpathsuffstatarray->Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectLengthSuffStat();
            lengthpathsuffstatarray->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
        }
        SlaveSendAdditive(*lengthpathsuffstatarray);
    }

    void MasterReceiveBranchLengthsSuffStat()   {
        lengthpathsuffstatarray->Clear();
        MasterReceiveAdditive(*lengthpathsuffstatarray);
    }

    void CollectGeneBranchLengthsSuffStat() {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectLengthSuffStat();
            // (*lengthpathsuffstattreearray)[gene].Clear();
            // (*lengthpathsuffstattreearray)[gene]->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
            (*lengthpathsuffstattreearray)[gene].BranchArray<PoissonSuffStat>::Copy(*geneprocess[gene]->GetLengthPathSuffStatArray());
        }
    }

    void SlaveSendGeneBranchLengthsSuffStat()   {
        CollectGeneBranchLengthsSuffStat();
        SlaveSendGeneArray(*lengthpathsuffstattreearray);
    }

    void MasterReceiveGeneBranchLengthsSuffStat()   {
        MasterReceiveGeneArray(*lengthpathsuffstattreearray);
    }

    void SlaveSendNucPathSuffStat() {
        nucpathsuffstat.Clear();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->CollectComponentPathSuffStat();
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
        }

        SlaveSendAdditive(nucpathsuffstat);
    }

    void MasterReceiveNucPathSuffStat() {
        nucpathsuffstat.Clear();
        MasterReceiveAdditive(nucpathsuffstat);
    }

    // hyper suffstats

    void SlaveSendBranchLengthsHyperSuffStat()  {
        CollectGeneBranchLengthsSuffStat();
        lengthhypersuffstatarray->Clear();
        lengthhypersuffstatarray->AddSuffStat(*branchlengtharray,*lengthpathsuffstattreearray);
        // lengthhypersuffstatarray->AddSuffStat(*branchlengtharray);
        SlaveSendAdditive(*lengthhypersuffstatarray);
    }

    void MasterReceiveBranchLengthsHyperSuffStat()  {
        lengthhypersuffstatarray->Clear();
        MasterReceiveAdditive(*lengthhypersuffstatarray);
    }

    void SlaveSendNucRatesHyperSuffStat()   {
        nucrelratesuffstat.Clear();
        nucrelratearray->AddSuffStat(nucrelratesuffstat);
        SlaveSendAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        nucstatarray->AddSuffStat(nucstatsuffstat);
        SlaveSendAdditive(nucstatsuffstat);
    }

    void MasterReceiveNucRatesHyperSuffStat()   {
        nucrelratesuffstat.Clear();
        MasterReceiveAdditive(nucrelratesuffstat);

        nucstatsuffstat.Clear();
        MasterReceiveAdditive(nucstatsuffstat);
    }

    void SlaveSendMixtureHyperSuffStat()    {
        puromsuffstat.Clear();
        puromarray->AddSuffStat(puromsuffstat);
        SlaveSendAdditive(puromsuffstat);

        dposomsuffstat.Clear();
        dposomsuffstat.AddSuffStat(*dposomarray, *poswarray);
        SlaveSendAdditive(dposomsuffstat);

        purwsuffstat.Clear();
        purwarray->AddSuffStat(purwsuffstat);
        SlaveSendAdditive(purwsuffstat);

        poswsuffstat.Clear();
        poswarray->AddSuffStat(poswsuffstat);
        SlaveSendAdditive(poswsuffstat);
    }

    void MasterReceiveMixtureHyperSuffStat()    {
        puromsuffstat.Clear();
        MasterReceiveAdditive(puromsuffstat);
        dposomsuffstat.Clear();
        MasterReceiveAdditive(dposomsuffstat);
        purwsuffstat.Clear();
        MasterReceiveAdditive(purwsuffstat);
        poswsuffstat.Clear();
        MasterReceiveAdditive(poswsuffstat);
    }

    // log likelihoods

    void SlaveSendLogProbs()    {
        GeneLogPrior = 0;
        GeneBLLogPrior = 0;
        GeneNucRatesLogPrior = 0;
        GeneOmegaLogPrior = 0;
        lnL = 0;
        moveTime = movechrono.GetTime();
        mapTime = mapchrono.GetTime();
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            GeneBLLogPrior += geneprocess[gene]->BranchLengthsLogPrior();
            GeneNucRatesLogPrior += geneprocess[gene]->NucRatesLogPrior();
            GeneOmegaLogPrior += geneprocess[gene]->OmegaLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(GeneBLLogPrior);
        SlaveSendAdditive(GeneNucRatesLogPrior);
        SlaveSendAdditive(GeneOmegaLogPrior);
        SlaveSendAdditive(lnL);
        SlaveSendAdditive(moveTime);
        SlaveSendAdditive(mapTime);
    }

    void MasterReceiveLogProbs()    {
        GeneLogPrior = 0;
        MasterReceiveAdditive(GeneLogPrior);
        GeneBLLogPrior = 0;
        MasterReceiveAdditive(GeneBLLogPrior);
        GeneNucRatesLogPrior = 0;
        MasterReceiveAdditive(GeneNucRatesLogPrior);
        GeneOmegaLogPrior = 0;
        MasterReceiveAdditive(GeneOmegaLogPrior);
        lnL = 0;
        MasterReceiveAdditive(lnL);
        moveTime = 0;
        mapTime = 0;
        MasterReceiveAdditive(moveTime);
        MasterReceiveAdditive(mapTime);
        moveTime /= (GetNprocs() - 1);
        mapTime /= (GetNprocs() - 1);
    }

    double GetSlaveMoveTime() const { return moveTime; }

    double GetSlaveMapTime() const { return mapTime; }

    double GetMasterMoveTime() const { return movechrono.GetTime(); }

    void StartChrono() { totchrono.Start(); }

};
