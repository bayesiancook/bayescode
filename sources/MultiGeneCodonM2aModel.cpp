
#include "MultiGeneCodonM2aModel.hpp"

//-------------------
// Constructing and setting the model
// ------------------

MultiGeneCodonM2aModel::MultiGeneCodonM2aModel(string indatapath, string indatafile, string intreefile,
                                               double inpihypermean, double inpihyperinvconc,
                                               int inmyid, int innprocs)
    :

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

void MultiGeneCodonM2aModel::Allocate() {
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
                if (GeneAlloc[gene] == myid)    {
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

void MultiGeneCodonM2aModel::FastUpdate() {
    branchlength->SetScale(lambda);
    if (blmode == 1) {
        branchlengtharray->SetShape(1.0 / blhyperinvshape);
    }
    nucrelratearray->SetConcentration(1.0 / nucrelratehyperinvconc);
    nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
    SetMixtureArrays();
}

void MultiGeneCodonM2aModel::MasterUpdate() {
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

void MultiGeneCodonM2aModel::SlaveUpdate() {
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

void MultiGeneCodonM2aModel::GeneUpdate() {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->Update();
    }
}

void MultiGeneCodonM2aModel::MasterPostPred(string name) {
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

void MultiGeneCodonM2aModel::SlavePostPred(string name) {
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

void MultiGeneCodonM2aModel::GenePostPred(string name) {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->PostPred(name + GetLocalGeneName(gene));
    }
}

void MultiGeneCodonM2aModel::SetAcrossGenesModes(int inblmode, int innucmode, int inpurommode,
                                                 int indposommode, int inpurwmode, int inposwmode) {
    blmode = inblmode;
    nucmode = innucmode;
    purommode = inpurommode;
    dposommode = indposommode;
    purwmode = inpurwmode;
    poswmode = inposwmode;
}

void MultiGeneCodonM2aModel::SetMixtureHyperParameters(
    double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean,
    double indposomhyperinvshape, double inpurwhypermean, double inpurwhyperinvconc,
    double inposwhypermean, double inposwhyperinvconc) {
    puromhypermean = inpuromhypermean;
    puromhyperinvconc = inpuromhyperinvconc;
    dposomhypermean = indposomhypermean;
    dposomhyperinvshape = indposomhyperinvshape;
    purwhypermean = inpurwhypermean;
    purwhyperinvconc = inpurwhyperinvconc;
    poswhypermean = inposwhypermean;
    poswhyperinvconc = inposwhyperinvconc;
}

void MultiGeneCodonM2aModel::UpdateNucMatrix() {
    nucmatrix->CopyStationary((*nucstatarray)[0]);
    nucmatrix->CorruptMatrix();
}

void MultiGeneCodonM2aModel::SetMixtureArrays() {
    double puromalpha = puromhypermean / puromhyperinvconc;
    double purombeta = (1 - puromhypermean) / puromhyperinvconc;
    puromarray->SetAlpha(puromalpha);
    puromarray->SetBeta(purombeta);

    double dposomalpha = 1.0 / dposomhyperinvshape;
    double dposombeta = dposomalpha / dposomhypermean;
    dposomarray->SetShape(dposomalpha);
    dposomarray->SetScale(dposombeta);
    if (myid) {
        dposomarray->PriorResample(*poswarray);
        // necessary after changing some dposom values
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->SetMixtureParameters((*puromarray)[gene], (*dposomarray)[gene],
                                                    (*purwarray)[gene], (*poswarray)[gene]);
        }
    }

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

//-------------------
// Streams
// ------------------

void MultiGeneCodonM2aModel::MasterFromStream(istream &is) {

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

void MultiGeneCodonM2aModel::MasterToStream(ostream &os) const {

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

//-------------------
// Traces and Monitors
// ------------------

void MultiGeneCodonM2aModel::TraceHeader(ostream &os) const {
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

void MultiGeneCodonM2aModel::Trace(ostream &os) const {
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

void MultiGeneCodonM2aModel::TracePosWeight(ostream &os) const {
    for (int gene = 0; gene < Ngene; gene++) {
        os << poswarray->GetVal(gene) << '\t';
    }
    os << '\n';
    os.flush();
}

void MultiGeneCodonM2aModel::TracePosOm(ostream &os) const {
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

int MultiGeneCodonM2aModel::GetNpos() const { return GetNgene() - poswarray->GetNullSet(); }

double MultiGeneCodonM2aModel::GetMeanTotalLength() const {
    double tot = 0;
    for (int j = 0; j < Nbranch; j++) {
        tot += branchlength->GetVal(j);
    }
    return tot;
}

double MultiGeneCodonM2aModel::GetMeanLength() const {
    if (blmode == 2) {
        cerr << "error: in getvarlength\n";
        exit(1);
    }

    return branchlengtharray->GetMeanLength();
}

double MultiGeneCodonM2aModel::GetVarLength() const {
    if (blmode == 2) {
        cerr << "error: in getvarlength\n";
        exit(1);
    }

    return branchlengtharray->GetVarLength();
}

double MultiGeneCodonM2aModel::GetVarNucRelRate() const {
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

double MultiGeneCodonM2aModel::GetVarNucStat() const {
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
// Log Priors and likelihood
// ------------------

double MultiGeneCodonM2aModel::GetLogPrior() const {
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
    // already accounted for by gene component
    // total += MixtureLogPrior();

    return total;
}

double MultiGeneCodonM2aModel::MixtureHyperLogPrior() const {
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
// Moves
// ------------------

void MultiGeneCodonM2aModel::MasterMove() {
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

// slave move
void MultiGeneCodonM2aModel::SlaveMove() {
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
        SetMixtureArrays();

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

void MultiGeneCodonM2aModel::GeneResampleSub(double frac) {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->ResampleSub(frac);
    }
}

void MultiGeneCodonM2aModel::MoveGeneParameters(int nrep) {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->MoveParameters(nrep);
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

void MultiGeneCodonM2aModel::ResampleBranchLengths() {
    branchlength->GibbsResample(*lengthpathsuffstatarray);
}

void MultiGeneCodonM2aModel::ResampleGeneBranchLengths()   {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->ResampleBranchLengths();
        geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
    }
}

void MultiGeneCodonM2aModel::MoveLambda() {
    hyperlengthsuffstat.Clear();
    hyperlengthsuffstat.AddSuffStat(*branchlength);
    ScalingMove(lambda, 1.0, 10, &MultiGeneCodonM2aModel::LambdaHyperLogProb,
                &MultiGeneCodonM2aModel::NoUpdate, this);
    ScalingMove(lambda, 0.3, 10, &MultiGeneCodonM2aModel::LambdaHyperLogProb,
                &MultiGeneCodonM2aModel::NoUpdate, this);
    branchlength->SetScale(lambda);
}

void MultiGeneCodonM2aModel::MoveBranchLengthsHyperParametersIntegrated() {

    BranchLengthsHyperScalingMoveIntegrated(1.0, 10);
    BranchLengthsHyperScalingMoveIntegrated(0.3, 10);

    ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneCodonM2aModel::BranchLengthsIntegratedHyperLogProb,
                &MultiGeneCodonM2aModel::NoUpdate, this);
    ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneCodonM2aModel::BranchLengthsIntegratedHyperLogProb,
                &MultiGeneCodonM2aModel::NoUpdate, this);

    branchlengtharray->SetShape(1.0 / blhyperinvshape);
    MoveLambda();
}

double MultiGeneCodonM2aModel::GeneBranchLengthsHyperSuffStatLogProb() const {
    double total = 0;
    for (int j=0; j<Nbranch; j++)   {
        total += GeneBranchLengthsHyperSuffStatLogProb(j);
    }
    return total;
}

double MultiGeneCodonM2aModel::GeneBranchLengthsHyperSuffStatLogProb(int branch) const {
    double total = 0;
    for (int gene=0; gene<GetLocalNgene(); gene++)    {
        total += GeneBranchLengthsHyperSuffStatLogProb(gene,branch);
    }
    return total;
}

double MultiGeneCodonM2aModel::GeneBranchLengthsHyperSuffStatLogProb(int gene, int branch)  const {

    const PoissonSuffStat &suffstat = lengthpathsuffstattreearray->GetVal(gene).GetVal(branch);
    int count = suffstat.GetCount();
    double b = suffstat.GetBeta();

    double alpha = 1.0 / blhyperinvshape;
    double beta = alpha / branchlength->GetVal(branch);
    /*
    double alpha = branchlengtharray->GetVal(gene).GetAlpha(branch);
    double beta = branchlengtharray->GetVal(gene).GetBeta(branch);
    */

    return alpha * log(beta) - Random::logGamma(alpha) + Random::logGamma(alpha + count) - (alpha + count) * log(beta + b);
}


double MultiGeneCodonM2aModel::BranchLengthsHyperScalingMoveIntegrated(double tuning, int nrep) {
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

void MultiGeneCodonM2aModel::MoveBranchLengthsHyperParameters() {

    BranchLengthsHyperScalingMove(1.0, 10);
    BranchLengthsHyperScalingMove(0.3, 10);

    ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneCodonM2aModel::BranchLengthsHyperLogProb,
                &MultiGeneCodonM2aModel::NoUpdate, this);
    ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneCodonM2aModel::BranchLengthsHyperLogProb,
                &MultiGeneCodonM2aModel::NoUpdate, this);

    branchlengtharray->SetShape(1.0 / blhyperinvshape);
    MoveLambda();
}

double MultiGeneCodonM2aModel::BranchLengthsHyperScalingMove(double tuning, int nrep) {
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

void MultiGeneCodonM2aModel::MoveNucRatesHyperParameters() {
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

void MultiGeneCodonM2aModel::MoveNucRates() {
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

void MultiGeneCodonM2aModel::MoveMixtureHyperParameters() {
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
}

void MultiGeneCodonM2aModel::MovePoswHyper()	{

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

int MultiGeneCodonM2aModel::PoswCompMove(double tuning)	{

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


void MultiGeneCodonM2aModel::ResamplePi() {
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
// MPI Send Receive
// ------------------

// branch lengths

void MultiGeneCodonM2aModel::MasterSendGlobalBranchLengths() { MasterSendGlobal(*branchlength); }

void MultiGeneCodonM2aModel::SlaveReceiveGlobalBranchLengths() {
    SlaveReceiveGlobal(*branchlength);
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetBranchLengths(*branchlength);
    }
}

void MultiGeneCodonM2aModel::MasterSendBranchLengthsHyperParameters() {
    MasterSendGlobal(*branchlength, blhyperinvshape);
}

void MultiGeneCodonM2aModel::SlaveReceiveBranchLengthsHyperParameters() {
    SlaveReceiveGlobal(*branchlength, blhyperinvshape);
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
    }
}

void MultiGeneCodonM2aModel::GeneResampleEmptyBranches()   {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->ResampleEmptyBranches();
        geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
    }
}

void MultiGeneCodonM2aModel::MasterSendGeneBranchLengths() {
    MasterSendGeneArray(*branchlengtharray);
}

void MultiGeneCodonM2aModel::SlaveReceiveGeneBranchLengths() {
    SlaveReceiveGeneArray(*branchlengtharray);
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetBranchLengths(branchlengtharray->GetVal(gene));
    }
}

void MultiGeneCodonM2aModel::SlaveSendGeneBranchLengths() {
    SlaveSendGeneArray(*branchlengtharray);
}

void MultiGeneCodonM2aModel::MasterReceiveGeneBranchLengths() {
    MasterReceiveGeneArray(*branchlengtharray);
}

void MultiGeneCodonM2aModel::SlaveSendBranchLengthsSuffStat() {
    lengthpathsuffstatarray->Clear();
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->CollectLengthSuffStat();
        lengthpathsuffstatarray->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
    }
    SlaveSendAdditive(*lengthpathsuffstatarray);
}

void MultiGeneCodonM2aModel::CollectGeneBranchLengthsSuffStat()    {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->CollectLengthSuffStat();
        // (*lengthpathsuffstattreearray)[gene].Clear();
        // (*lengthpathsuffstattreearray)[gene]->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
        (*lengthpathsuffstattreearray)[gene].BranchArray<PoissonSuffStat>::Copy(*geneprocess[gene]->GetLengthPathSuffStatArray());
    }
}

void MultiGeneCodonM2aModel::SlaveSendGeneBranchLengthsSuffStat() {
    CollectGeneBranchLengthsSuffStat();
    SlaveSendGeneArray(*lengthpathsuffstattreearray);
}


void MultiGeneCodonM2aModel::MasterReceiveBranchLengthsSuffStat() {
    lengthpathsuffstatarray->Clear();
    MasterReceiveAdditive(*lengthpathsuffstatarray);
}

void MultiGeneCodonM2aModel::MasterReceiveGeneBranchLengthsSuffStat() {
    MasterReceiveGeneArray(*lengthpathsuffstattreearray);
}

void MultiGeneCodonM2aModel::SlaveSendBranchLengthsHyperSuffStat() {
    CollectGeneBranchLengthsSuffStat();
    lengthhypersuffstatarray->Clear();
    lengthhypersuffstatarray->AddSuffStat(*branchlengtharray,*lengthpathsuffstattreearray);
    // lengthhypersuffstatarray->AddSuffStat(*branchlengtharray);
    SlaveSendAdditive(*lengthhypersuffstatarray);
}

void MultiGeneCodonM2aModel::MasterReceiveBranchLengthsHyperSuffStat() {
    lengthhypersuffstatarray->Clear();
    MasterReceiveAdditive(*lengthhypersuffstatarray);
}

// Nuc Rates

void MultiGeneCodonM2aModel::MasterSendGlobalNucRates() {
    MasterSendGlobal(nucrelratearray->GetVal(0), nucstatarray->GetVal(0));
}

void MultiGeneCodonM2aModel::SlaveReceiveGlobalNucRates() {
    SlaveReceiveGlobal((*nucrelratearray)[0], (*nucstatarray)[0]);

    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetNucRates((*nucrelratearray)[0], (*nucstatarray)[0]);
    }
}

void MultiGeneCodonM2aModel::MasterSendGeneNucRates() {
    MasterSendGeneArray(*nucrelratearray, *nucstatarray);
}

void MultiGeneCodonM2aModel::SlaveReceiveGeneNucRates() {
    SlaveReceiveGeneArray(*nucrelratearray, *nucstatarray);

    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetNucRates((*nucrelratearray)[gene], (*nucstatarray)[gene]);
    }
}

void MultiGeneCodonM2aModel::SlaveSendGeneNucRates() {
    SlaveSendGeneArray(*nucrelratearray, *nucstatarray);
}

void MultiGeneCodonM2aModel::MasterReceiveGeneNucRates() {
    MasterReceiveGeneArray(*nucrelratearray, *nucstatarray);
}

void MultiGeneCodonM2aModel::SlaveReceiveNucRatesHyperParameters() {
    SlaveReceiveGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
    SlaveReceiveGlobal(nucstathypercenter, nucstathyperinvconc);

    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter, nucrelratehyperinvconc,
                                                      nucstathypercenter, nucstathyperinvconc);
    }
}

void MultiGeneCodonM2aModel::SlaveSendNucRatesHyperSuffStat() {
    nucrelratesuffstat.Clear();
    nucrelratearray->AddSuffStat(nucrelratesuffstat);
    SlaveSendAdditive(nucrelratesuffstat);

    nucstatsuffstat.Clear();
    nucstatarray->AddSuffStat(nucstatsuffstat);
    SlaveSendAdditive(nucstatsuffstat);
}

void MultiGeneCodonM2aModel::MasterReceiveNucRatesHyperSuffStat() {
    nucrelratesuffstat.Clear();
    MasterReceiveAdditive(nucrelratesuffstat);

    nucstatsuffstat.Clear();
    MasterReceiveAdditive(nucstatsuffstat);
}

void MultiGeneCodonM2aModel::SlaveSendNucPathSuffStat() {
    nucpathsuffstat.Clear();
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->CollectComponentPathSuffStat();
        geneprocess[gene]->CollectNucPathSuffStat();
        nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
    }

    SlaveSendAdditive(nucpathsuffstat);
}

void MultiGeneCodonM2aModel::MasterReceiveNucPathSuffStat() {
    nucpathsuffstat.Clear();
    MasterReceiveAdditive(nucpathsuffstat);
}

// Mixture params and hyperparams

void MultiGeneCodonM2aModel::MasterSendMixtureHyperParameters() { MasterSendGlobal(mixhyperparam); }

void MultiGeneCodonM2aModel::SlaveReceiveMixtureHyperParameters() {
    SlaveReceiveGlobal(mixhyperparam);

    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetMixtureHyperParameters(
            puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape, pi,
            purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc);
    }
}

void MultiGeneCodonM2aModel::MasterSendNucRatesHyperParameters() {
    MasterSendGlobal(nucrelratehypercenter, nucrelratehyperinvconc);
    MasterSendGlobal(nucstathypercenter, nucstathyperinvconc);
}
void MultiGeneCodonM2aModel::SlaveSendMixtureHyperSuffStat() {
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

void MultiGeneCodonM2aModel::MasterReceiveMixtureHyperSuffStat() {
    puromsuffstat.Clear();
    MasterReceiveAdditive(puromsuffstat);
    dposomsuffstat.Clear();
    MasterReceiveAdditive(dposomsuffstat);
    purwsuffstat.Clear();
    MasterReceiveAdditive(purwsuffstat);
    poswsuffstat.Clear();
    MasterReceiveAdditive(poswsuffstat);
}

void MultiGeneCodonM2aModel::SlaveSendMixture() {
    SlaveSendGeneArray(*puromarray, *dposomarray);
    SlaveSendGeneArray(*purwarray, *poswarray);
}

void MultiGeneCodonM2aModel::MasterReceiveMixture() {
    MasterReceiveGeneArray(*puromarray, *dposomarray);
    MasterReceiveGeneArray(*purwarray, *poswarray);
}

void MultiGeneCodonM2aModel::MasterSendMixture() {
    MasterSendGeneArray(*puromarray, *dposomarray);
    MasterSendGeneArray(*purwarray, *poswarray);
}

void MultiGeneCodonM2aModel::SlaveReceiveMixture() {
    SlaveReceiveGeneArray(*puromarray, *dposomarray);
    SlaveReceiveGeneArray(*purwarray, *poswarray);

    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetMixtureParameters((*puromarray)[gene], (*dposomarray)[gene],
                                                (*purwarray)[gene], (*poswarray)[gene]);
    }
}

void MultiGeneCodonM2aModel::SlaveSendLogProbs() {
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

void MultiGeneCodonM2aModel::MasterReceiveLogProbs() {
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

void MultiGeneCodonM2aModel::MasterTraceSitesPostProb(ostream &os) {
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
            exit(1);
        }
        delete[] array;
    }
    os << '\n';
    os.flush();
}

void MultiGeneCodonM2aModel::SlaveTraceSitesPostProb() {
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
