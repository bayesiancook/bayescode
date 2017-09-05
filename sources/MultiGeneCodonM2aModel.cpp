
#include "MultiGeneCodonM2aModel.hpp"


//-------------------
// Constructing and setting the model
// ------------------

MultiGeneCodonM2aModel::MultiGeneCodonM2aModel(string datafile, string intreefile, double inpihypermean, double inpihyperinvconc, int inmyid, int innprocs) : 

    MultiGeneMPIModule(inmyid,innprocs), 
    mixhyperparam(9,0), 
    puromhypermean(mixhyperparam[0]), 
    puromhyperinvconc(mixhyperparam[1]),
    dposomhypermean(mixhyperparam[2]),
    dposomhyperinvshape(mixhyperparam[3]),
    purwhypermean(mixhyperparam[4]),
    purwhyperinvconc(mixhyperparam[5]),
    poswhypermean(mixhyperparam[6]),
    poswhyperinvconc(mixhyperparam[7]),
    pi(mixhyperparam[8]),
    nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {

    burnin = 0;

    pihypermean = inpihypermean;
    pihyperinvconc = inpihyperinvconc;
    pi = pihypermean;

    AllocateAlignments(datafile);
    treefile = intreefile;

    // all datafiles have all taxa (with missing data if needed) in same order
    // makes it easier to register tree with data, etc.

    refcodondata = new CodonSequenceAlignment(refdata, true);
    taxonset = refdata->GetTaxonSet();
    Ntaxa = refdata->GetNtaxa();

    // get tree from file (newick format)
    tree = new Tree(treefile);

    // check whether tree and data fits together
    tree->RegisterWith(taxonset);

    tree->SetIndices();
    Nbranch = tree->GetNbranch();

    if (! myid) {
        std::cerr << "number of taxa : " << Ntaxa << '\n';
        std::cerr << "number of branches : " << Nbranch << '\n';
        std::cerr << "-- Tree and data fit together\n";
    }
}

void MultiGeneCodonM2aModel::Allocate() {

    nucstatacc1 = nucstatacc2 = nucstattot1 = nucstattot2 = 0;
    nucrracc1 = nucrracc2 = nucrrtot1 = nucrrtot2 = 0;

    lambda = 10;
    branchlength = new BranchIIDGamma(*tree,1.0,lambda);
    blhyperinvshape = 1.0;
    if (blmode == 2)    {
        lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        lengthhypersuffstatarray = 0;
    }
    else    {
        branchlength->SetAllBranches(1.0/lambda);
        branchlengtharray = new GammaWhiteNoiseArray(GetLocalNgene(),*tree,*branchlength,1.0/blhyperinvshape);
        lengthsuffstatarray = 0;
        lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
    }

    nucrelratehypercenter.assign(Nrr,1.0/Nrr);
    nucrelratehyperinvconc = 1.0 / Nrr;

    nucstathypercenter.assign(Nnuc,1.0/Nnuc);
    nucstathyperinvconc = 1.0 / Nnuc;

    if (nucmode == 2)   {
        nucrelratearray = new IIDDirichlet(1,nucrelratehypercenter,1.0/nucrelratehyperinvconc);
        nucstatarray = new IIDDirichlet(1,nucstathypercenter,1.0/nucstathyperinvconc);
        nucmatrix = new GTRSubMatrix(Nnuc,(*nucrelratearray)[0],(*nucstatarray)[0],true);
    }
    else    {
        nucrelratearray = new IIDDirichlet(GetLocalNgene(),nucrelratehypercenter,1.0/nucrelratehyperinvconc);
        nucstatarray = new IIDDirichlet(GetLocalNgene(),nucstathypercenter,1.0/nucstathyperinvconc);
        nucmatrix = 0;
    }

    double puromalpha = puromhypermean / puromhyperinvconc;
    double purombeta = (1-puromhypermean) / puromhyperinvconc;
    puromarray = new IIDBeta(GetLocalNgene(),puromalpha,purombeta);

    double dposomalpha = 1.0 / dposomhyperinvshape;
    double dposombeta = dposomalpha / dposomhypermean;
    dposomarray = new IIDGamma(GetLocalNgene(),dposomalpha,dposombeta);

    double purwalpha = purwhypermean / purwhyperinvconc;
    double purwbeta = (1-purwhypermean) / purwhyperinvconc;
    purwarray = new IIDBeta(GetLocalNgene(),purwalpha,purwbeta);

    double poswalpha = poswhypermean / poswhyperinvconc;
    double poswbeta = (1-poswhypermean) / poswhyperinvconc;
    poswarray = new IIDBernoulliBeta(GetLocalNgene(),pi,poswalpha,poswbeta);

    if (! GetMyid())    {
        geneprocess.assign(0,(CodonM2aModel*) 0);
    }
    else    {
        geneprocess.assign(GetLocalNgene(),(CodonM2aModel*) 0);

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene] = new CodonM2aModel(GetLocalGeneName(gene),treefile,pi);
            geneprocess[gene]->SetAcrossGenesModes(blmode,nucmode);
        }
    }
}

void MultiGeneCodonM2aModel::Unfold()   {

    if (! GetMyid())    {

        MasterSendBranchLengthsHyperParameters();
        MasterSendNucRatesHyperParameters();
        MasterSendMixtureHyperParameters();

        if (blmode == 2)    {
            MasterSendGlobalBranchLengths();
        }
        else    {
            MasterSendGeneBranchLengths();
        }

        if (nucmode == 2)   {
            MasterSendGlobalNucRates();
        }
        else    {
            MasterSendGeneNucRates();
        }

        MasterSendMixture();
        MasterReceiveLogProbs();
    }
    else    {

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->Allocate();
        }

        SlaveReceiveBranchLengthsHyperParameters();
        SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveMixtureHyperParameters();

        if (blmode == 2)    {
            SlaveReceiveGlobalBranchLengths();
        }
        else    {
            SlaveReceiveGeneBranchLengths();
        }

        if (nucmode == 2)   {
            SlaveReceiveGlobalNucRates();
        }
        else    {
            SlaveReceiveGeneNucRates();
        }

        SlaveReceiveMixture();

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->UpdateMatrices();
            geneprocess[gene]->Unfold();
        }

        SlaveSendLogProbs();
    }
}

void MultiGeneCodonM2aModel::SetAcrossGenesModes(int inblmode, int innucmode, int inpurommode, int indposommode, int inpurwmode, int inposwmode)  {
    blmode = inblmode;
    nucmode = innucmode;
    purommode = inpurommode;
    dposommode = indposommode;
    purwmode = inpurwmode;
    poswmode = inposwmode;
}

void MultiGeneCodonM2aModel::SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape, double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc)   {
    puromhypermean = inpuromhypermean;
    puromhyperinvconc = inpuromhyperinvconc;
    dposomhypermean = indposomhypermean;
    dposomhyperinvshape = indposomhyperinvshape;
    purwhypermean = inpurwhypermean;
    purwhyperinvconc = inpurwhyperinvconc;
    poswhypermean = inposwhypermean;
    poswhyperinvconc = inposwhyperinvconc;
}

void MultiGeneCodonM2aModel::UpdateNucMatrix()	{
    nucmatrix->CopyStationary((*nucstatarray)[0]);
    nucmatrix->CorruptMatrix();
}

void MultiGeneCodonM2aModel::SetMixtureArrays()    {

    double puromalpha = puromhypermean / puromhyperinvconc;
    double purombeta = (1-puromhypermean) / puromhyperinvconc;
    puromarray->SetAlpha(puromalpha);
    puromarray->SetBeta(purombeta);

    double dposomalpha = 1.0 / dposomhyperinvshape;
    double dposombeta = dposomalpha / dposomhypermean;
    dposomarray->SetShape(dposomalpha);
    dposomarray->SetScale(dposombeta);
    dposomarray->PriorResample(*poswarray);
    // necessary after changing some dposom values
    for (int gene=0; gene<GetLocalNgene(); gene++)    {
        geneprocess[gene]->SetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
    }

    double purwalpha = purwhypermean / purwhyperinvconc;
    double purwbeta = (1-purwhypermean) / purwhyperinvconc;
    purwarray->SetAlpha(purwalpha);
    purwarray->SetBeta(purwbeta);

    double poswalpha = poswhypermean / poswhyperinvconc;
    double poswbeta = (1-poswhypermean) / poswhyperinvconc;
    poswarray->SetPi(pi);
    poswarray->SetAlpha(poswalpha);
    poswarray->SetBeta(poswbeta);
}


//-------------------
// Traces and Monitors
// ------------------


void MultiGeneCodonM2aModel::TraceHeader(ostream& os)   {

    os << "#logprior\tlnL";
    if (blmode == 2)    {
        os << "\tlength";
    }
    else    {
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
    if (nucmode != 2)   {
        os << "\tstdevrr\tcenter\thyperinvconc";
        os << "\tstdevstat\tcenter\thyperinvconc";
        if (nucmode == 1)   {
            os << "\tnucrracc1\tnucrracc2\tnucrracc3";
            os << "\tnucstatacc1\tnucstatacc2\tnucstatacc3";
            os << "\tnucrrlogprob\tnucstatlogprob";
            os << "\tnucrrsuffstatlogprob\tnucstatsuffstatlogprob";
        }
    }
    os << '\n';
}

void MultiGeneCodonM2aModel::Trace(ostream& os)    {
    os << GetLogPrior();
    os << '\t' << GetLogLikelihood();
    if (blmode == 2)    {
        os << '\t' << GetMeanTotalLength();
    }
    else    {
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
    if (nucmode != 2)   {
        os << '\t' << sqrt(GetVarNucRelRate()) << '\t' << Random::GetEntropy(nucrelratehypercenter) << '\t' << nucrelratehyperinvconc;
        os << '\t' << sqrt(GetVarNucStat()) << '\t' << Random::GetEntropy(nucstathypercenter) << '\t' << nucstathyperinvconc;
        if (nucmode == 1)   {
            os << '\t' << 100*((double) nucrracc1) / nucrrtot1;
            os << '\t' << 100*((double) nucrracc2) / nucrrtot2;
            os << '\t' << 100*((double) nucrracc3) / nucrrtot3;
            os << '\t' << 100*((double) nucstatacc1) / nucstattot1;
            os << '\t' << 100*((double) nucstatacc2) / nucstattot2;
            os << '\t' << 100*((double) nucstatacc3) / nucstattot3;
            nucrracc1 = nucrracc2 = nucrracc3 = nucrrtot1 = nucrrtot2 = nucrrtot3 = 0;
            nucstatacc1 = nucstatacc2 = nucstatacc3 = nucstattot1 = nucstattot2 = nucstattot3 = 0;
            os << '\t' << nucrelratearray->GetLogProb();
            os << '\t' << nucstatarray->GetLogProb();
            os << '\t' << nucrelratesuffstat.GetLogProb(nucrelratehypercenter,1.0/nucrelratehyperinvconc);
            os << '\t' << nucstatsuffstat.GetLogProb(nucstathypercenter,1.0/nucstathyperinvconc);
        }
    }
    os << '\n';
    os.flush();
}

void MultiGeneCodonM2aModel::TracePosWeight(ostream& os) {

    for (int gene=0; gene<Ngene; gene++)    {
        os << (*poswarray)[gene] << '\t';
    }
    os << '\n';
    os.flush();
}

void MultiGeneCodonM2aModel::TracePosOm(ostream& os) {

    for (int gene=0; gene<Ngene; gene++)    {
        os << 1 + (*dposomarray)[gene] << '\t';
    }
    os << '\n';
    os.flush();
}


int MultiGeneCodonM2aModel::GetNpos()    {
    return GetNgene() - poswarray->GetNullSet();
}

double MultiGeneCodonM2aModel::GetMeanTotalLength()	{
    double tot = 0;
    for (int j=1; j<Nbranch; j++)	{
        tot += branchlength->GetVal(j);
    }
    return tot;
}

double MultiGeneCodonM2aModel::GetMeanLength()   {

    if (blmode == 2)    {
        cerr << "error: in getvarlength\n";
        exit(1);
    }

    return branchlengtharray->GetMeanLength();
}

double MultiGeneCodonM2aModel::GetVarLength()   {

    if (blmode == 2)    {
        cerr << "error: in getvarlength\n";
        exit(1);
    }

    return branchlengtharray->GetVarLength();
}

void MultiGeneCodonM2aModel::SetNucRelRateCenterToMean()   {

    for (int j=0; j<Nrr; j++)   {
        double mean = 0;
        for (int g=0; g<Ngene; g++) {
            double tmp = (*nucrelratearray)[g][j];
            mean += tmp;
        }
        mean /= Ngene;
        nucrelratehypercenter[j] = mean;
    }
    nucrelratehyperinvconc = GetVarNucRelRate();
}

void MultiGeneCodonM2aModel::SetNucStatCenterToMean()   {

    for (int j=0; j<Nnuc; j++)   {
        double mean = 0;
        for (int g=0; g<Ngene; g++) {
            double tmp = (*nucstatarray)[g][j];
            mean += tmp;
        }
        mean /= Ngene;
        nucstathypercenter[j] = mean;
    }
    nucstathyperinvconc = GetVarNucStat();
}

double MultiGeneCodonM2aModel::GetVarNucRelRate()   {

    if (nucmode == 2)   {
        cerr << "error in getvarnucrelrate\n";
        exit(1);
    }

    double tot = 0;
    for (int j=0; j<Nrr; j++)   {
        double mean = 0;
        double var = 0;
        for (int g=0; g<Ngene; g++) {
            double tmp = (*nucrelratearray)[g][j];
            mean += tmp;
            var += tmp*tmp;
        }
        mean /= Ngene;
        var /= Ngene;
        var -= mean*mean;
        tot += var;
    }
    tot /= Nrr;
    return tot;
}

double MultiGeneCodonM2aModel::GetVarNucStat()  {

    if (nucmode == 2)   {
        cerr << "error in getvarnucstat\n";
        exit(1);
    }

    double tot = 0;
    for (int j=0; j<Nnuc; j++)   {
        double mean = 0;
        double var = 0;
        for (int g=0; g<Ngene; g++) {
            double tmp = (*nucstatarray)[g][j];
            mean += tmp;
            var += tmp*tmp;
        }
        mean /= Ngene;
        var /= Ngene;
        var -= mean*mean;
        tot += var;
    }
    tot /= Nnuc;
    return tot;
}


//-------------------
// Log Priors and likelihood
// ------------------

double MultiGeneCodonM2aModel::GetLogPrior()    {

    // gene contributions
    double total = GeneLogPrior;

    // branch lengths
    if (blmode == 2)    {
        total += GlobalBranchLengthsLogPrior();
    }
    else if (blmode == 1)   {
        total += GeneBranchLengthsHyperLogPrior();
    }
    else    {
        // nothing: everything accounted for by gene component
    }

    // nuc rates
    if (nucmode == 2)   {
        total += GlobalNucRatesLogPrior();
    }
    else if (nucmode == 1)  {
        total += GeneNucRatesHyperLogPrior();
    }
    else    {
        // nothing: everything accounted for by gene component
    }

    // mixture
    total += MixtureHyperLogPrior();
    // already accounted for by gene component
    // total += MixtureLogPrior();

    return total;
}

double MultiGeneCodonM2aModel::GlobalBranchLengthsLogPrior()    {

    double total = 0;
    total += LambdaHyperLogPrior();
    total += branchlength->GetLogProb();
    return total;
}

double MultiGeneCodonM2aModel::GeneBranchLengthsHyperLogPrior() {

    double total = 0;
    total += BranchLengthsHyperInvShapeLogPrior();
    total += branchlength->GetLogProb();
    return total;
}

double MultiGeneCodonM2aModel::LambdaHyperLogPrior()    {
    return -lambda/10;
}

double MultiGeneCodonM2aModel::BranchLengthsHyperInvShapeLogPrior()    {

    return -blhyperinvshape;
}

double MultiGeneCodonM2aModel::GeneNucRatesHyperLogPrior()  {

    double total = 0;
    if (nucmode == 1)   {
        total -= nucrelratehyperinvconc;
        total -= nucstathyperinvconc;
    }
    if (isnan(total))   {
        cerr << "error in NucRatesHyperLogPrior: nan\n";
        exit(1);
    }
    if (isinf(total))   {
        cerr << "error in NucRatesHyperLogPrior: inf\n";
        exit(1);
    }
    return total;
}

double MultiGeneCodonM2aModel::GlobalNucRatesLogPrior()    {

    double total = 0;
    // Dirichlet prior on relrates and stationaries
    total += nucrelratearray->GetLogProb();
    total += nucstatarray->GetLogProb();
    if (isnan(total))   {
        cerr << "error in NucRatesLogPrior: nan\n";
        exit(1);
    }
    if (isinf(total))   {
        cerr << "error in NucRatesLogPrior: inf\n";
        exit(1);
    }
    return total;
}

double MultiGeneCodonM2aModel::MixtureHyperLogPrior()   {

    double total = 0;
    if (pi) {
        double pialpha = pihypermean / pihyperinvconc;
        double pibeta = (1-pihypermean) / pihyperinvconc;
        total += (pialpha-1) * log(1.0 - pi) + (pibeta-1) * log(pi);
    }
    total -= puromhyperinvconc;
    total -= purwhyperinvconc;
    total -= 10*poswhyperinvconc;
    total -= dposomhypermean;
    total -= 10*dposomhyperinvshape;
    return total;
}

//-------------------
// SuffStatLogProbs
// ------------------

double MultiGeneCodonM2aModel::NucRatesSuffStatLogProb() {

    return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
}

//-------------------
// Hyper SuffStatLogProbs
// ------------------

double MultiGeneCodonM2aModel::LambdaHyperSuffStatLogProb()	{

    return lambdasuffstat.GetLogProb(1.0,lambda);
}

double MultiGeneCodonM2aModel::BranchLengthsHyperSuffStatLogProb()   {

    return lengthhypersuffstatarray->GetLogProb(*branchlength,blhyperinvshape);
}

double MultiGeneCodonM2aModel::NucRatesHyperSuffStatLogProb()   {

    double total = 0;
    total += nucrelratesuffstat.GetLogProb(nucrelratehypercenter,1.0/nucrelratehyperinvconc);
    total += nucstatsuffstat.GetLogProb(nucstathypercenter,1.0/nucstathyperinvconc);
    if (isnan(total))   {
        cerr << "error in NucRatesHyperSuffStatLogProb: nan\n";
        exit(1);
    }
    if (isinf(total))   {
        cerr << "error in NucRatesHyperSuffStatLogProb: inf\n";
        exit(1);
    }
    return total;
}

double MultiGeneCodonM2aModel::MixtureHyperSuffStatLogProb()   {

    double total = 0;

    double puromalpha = puromhypermean / puromhyperinvconc;
    double purombeta = (1-puromhypermean) / puromhyperinvconc;
    total += puromsuffstat.GetLogProb(puromalpha,purombeta);

    double dposomalpha = 1.0 / dposomhyperinvshape;
    double dposombeta = dposomalpha / dposomhypermean;
    total += dposomsuffstat.GetLogProb(dposomalpha,dposombeta);

    double purwalpha = purwhypermean / purwhyperinvconc;
    double purwbeta = (1-purwhypermean) / purwhyperinvconc;
    total += purwsuffstat.GetLogProb(purwalpha,purwbeta);

    double poswalpha = poswhypermean / poswhyperinvconc;
    double poswbeta = (1-poswhypermean) / poswhyperinvconc;
    total += poswsuffstat.GetLogProb(pi,poswalpha,poswbeta);

    if (isnan(total))   {
        cerr << "hyper suff stat log prob is nan\n";
        cerr << puromsuffstat.GetLogProb(puromalpha,purombeta) << '\n';
        cerr << dposomsuffstat.GetLogProb(dposomalpha,dposombeta) << '\n';
        cerr << purwsuffstat.GetLogProb(purwalpha,purwbeta) << '\n';
        cerr << poswsuffstat.GetLogProb(pi,poswalpha,poswbeta) << '\n';
        cerr << pi << '\t' << poswalpha << '\t' << poswbeta << '\n';
        exit(1);
    }

    return total;
}

//-------------------
// Moves
// ------------------


void MultiGeneCodonM2aModel::MasterMove() {

    int nrep = 30;

    for (int rep=0; rep<nrep; rep++)	{

        // mixture hyperparameters
        MasterReceiveMixtureHyperSuffStat();
        MasterMoveMixtureHyperParameters();
        MasterSendMixtureHyperParameters();

        // global branch lengths, or gene branch lengths hyperparameters
        if (blmode == 2)    {
            MasterReceiveBranchLengthsSuffStat();
            MasterResampleBranchLengths();
            MasterMoveLambda();
            MasterSendGlobalBranchLengths();
        }
        else if (blmode == 1)    {
            MasterReceiveBranchLengthsHyperSuffStat();
            MasterMoveBranchLengthsHyperParameters();
            MasterSendBranchLengthsHyperParameters();
        }

        // global nucrates, or gene nucrates hyperparameters
        if (nucmode == 2)   {
            MasterReceiveNucPathSuffStat();
            MasterMoveNucRates();
            MasterSendGlobalNucRates();
        }
        else if (nucmode == 1)  {
            MasterReceiveNucRatesHyperSuffStat();
            MasterMoveNucRatesHyperParameters();
            MasterSendNucRatesHyperParameters();
        }
    }
    burnin++;
    if (blmode != 2)    {
        MasterReceiveGeneBranchLengths();
    }
    if (nucmode != 2)   {
        MasterReceiveGeneNucRates();
    }
    MasterReceiveMixture();
    MasterReceiveLogProbs();
}

// slave move
void MultiGeneCodonM2aModel::SlaveMove() {

    SlaveResampleSub(1.0);

    int nrep = 30;

    for (int rep=0; rep<nrep; rep++)	{

        // gene specific mixture parameters
        // possibly branch lengths and nuc rates (if mode == 1 or 2)
        SlaveMoveGeneParameters(1.0);

        // mixture hyperparameters
        SlaveSendMixtureHyperSuffStat();
        SlaveReceiveMixtureHyperParameters();
        SetMixtureArrays();

        // global branch lengths, or gene branch lengths hyperparameters
        if (blmode == 2)    {
            SlaveSendBranchLengthsSuffStat();
            SlaveReceiveGlobalBranchLengths();
        }
        else if (blmode == 1)   {
            SlaveSendBranchLengthsHyperSuffStat();
            SlaveReceiveBranchLengthsHyperParameters();
        }

        // global nucrates, or gene nucrates hyperparameters
        if (nucmode == 2)   {
            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalNucRates();
        }
        else if (nucmode == 1)  {
            SlaveSendNucRatesHyperSuffStat();
            SlaveReceiveNucRatesHyperParameters();
        }
    }
    burnin++;

    // collect current state 
    if (blmode != 2)    {
        SlaveSendGeneBranchLengths();
    }
    if (nucmode != 2)   {
        SlaveSendGeneNucRates();
    }
    SlaveSendMixture();
    SlaveSendLogProbs();
}

void MultiGeneCodonM2aModel::SlaveResampleSub(double frac)  {
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->ResampleSub(frac);
    }
}

void MultiGeneCodonM2aModel::SlaveMoveGeneParameters(int nrep)  {
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->MoveParameters(nrep);
        geneprocess[gene]->GetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
        if (blmode != 2) {
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
        if (nucmode != 2)    {
            geneprocess[gene]->GetNucRates((*nucrelratearray)[gene],(*nucstatarray)[gene]);
        }
    }
}

void MultiGeneCodonM2aModel::MasterResampleBranchLengths()    {
    branchlength->GibbsResample(*lengthsuffstatarray);
}

void MultiGeneCodonM2aModel::MasterMoveLambda()	{

    lambdasuffstat.Clear();
    branchlength->AddSuffStat(lambdasuffstat);
    MoveLambda(1.0,10);
    MoveLambda(0.3,10);
    branchlength->SetScale(lambda);
}

double MultiGeneCodonM2aModel::MoveLambda(double tuning, int nrep)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - LambdaHyperLogPrior() - LambdaHyperSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        lambda *= e;
        deltalogprob += LambdaHyperLogPrior() + LambdaHyperSuffStatLogProb();
        deltalogprob += m;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            lambda /= e;
        }
        ntot++;
    }
    return nacc/ntot;
}

void MultiGeneCodonM2aModel::MasterMoveBranchLengthsHyperParameters()   {

    for (int j=0; j<Nbranch; j++)   {
        BranchLengthsHyperScalingMove(1.0,10);
        BranchLengthsHyperScalingMove(0.3,10);
    }
    BranchLengthsHyperInvShapeMove(1.0,10);
    BranchLengthsHyperInvShapeMove(0.3,10);
    branchlengtharray->SetShape(1.0 / blhyperinvshape);
    MasterMoveLambda();
}

double MultiGeneCodonM2aModel::BranchLengthsHyperScalingMove(double tuning, int nrep)  {

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        for (int j=0; j<Nbranch; j++)   {
            double deltalogprob = - branchlength->GetLogProb(j) - lengthhypersuffstatarray->GetVal(j).GetLogProb(1.0/blhyperinvshape, 1.0 / blhyperinvshape / branchlength->GetVal(j));
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            (*branchlength)[j] *= e;
            deltalogprob += branchlength->GetLogProb(j) + lengthhypersuffstatarray->GetVal(j).GetLogProb(1.0/blhyperinvshape, 1.0 / blhyperinvshape / branchlength->GetVal(j));
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted)	{
                nacc ++;
            }
            else	{
                (*branchlength)[j] /= e;
            }
            ntot++;
        }
    }
    return nacc/ntot;
}

double MultiGeneCodonM2aModel::BranchLengthsHyperInvShapeMove(double tuning, int nrep)    {

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - BranchLengthsHyperInvShapeLogPrior() - BranchLengthsHyperSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        blhyperinvshape *= e;
        deltalogprob += BranchLengthsHyperInvShapeLogPrior() + BranchLengthsHyperSuffStatLogProb();
        deltalogprob += m;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            blhyperinvshape /= e;
        }
        ntot++;
    }
    return nacc/ntot;
}

double MultiGeneCodonM2aModel::NucRatesHyperProfileMove(vector<double>& x, double tuning, int n, int nrep)	{

    double nacc = 0;
    double ntot = 0;
    vector<double> bk(x.size(),0);
    for (int rep=0; rep<nrep; rep++)	{
        bk = x;
        double deltalogprob = - GeneNucRatesHyperLogPrior() - NucRatesHyperSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(x,x.size(),tuning,n);
        deltalogprob += GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
        deltalogprob += loghastings;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            x = bk;
        }
        ntot++;
    }
    return nacc/ntot;
}

double MultiGeneCodonM2aModel::NucRatesHyperScalingMove(double& x, double tuning, int nrep)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - GeneNucRatesHyperLogPrior() - NucRatesHyperSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        x *= e;
        deltalogprob += GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
        deltalogprob += m;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            x /= e;
        }
        ntot++;
    }
    return nacc/ntot;
}

void MultiGeneCodonM2aModel::MasterMoveNucRatesHyperParameters()    {

    /*
    if (burnin < 10)    {

        SetNucRelRateCenterToMean();
        SetNucStatCenterToMean();
    }
    else    {
    */

    NucRatesHyperProfileMove(nucrelratehypercenter,1.0,1,10);
    NucRatesHyperProfileMove(nucrelratehypercenter,0.3,1,10);
    NucRatesHyperProfileMove(nucrelratehypercenter,0.1,3,10);
    nucrracc1 += NucRatesHyperScalingMove(nucrelratehyperinvconc,1.0,10);
    nucrracc2 += NucRatesHyperScalingMove(nucrelratehyperinvconc,0.3,10);
    nucrracc3 += NucRatesHyperScalingMove(nucrelratehyperinvconc,0.03,10);
    nucrrtot1++;
    nucrrtot2++;
    nucrrtot3++;


    NucRatesHyperProfileMove(nucstathypercenter,1.0,1,10);
    NucRatesHyperProfileMove(nucstathypercenter,0.3,1,10);
    NucRatesHyperProfileMove(nucstathypercenter,0.1,2,10);
    /*
    nucstatacc1 += NucRatesHyperProfileMove(nucstathypercenter,GetVarNucStat(),0,10);
    nucstatacc2 += NucRatesHyperProfileMove(nucstathypercenter,0.3*GetVarNucStat(),0,10);
    nucstatacc3 += NucRatesHyperProfileMove(nucstathypercenter,0.1*GetVarNucStat(),0,10);
    */
    NucRatesHyperScalingMove(nucstathyperinvconc,1.0,10);
    NucRatesHyperScalingMove(nucstathyperinvconc,0.3,10);
    NucRatesHyperScalingMove(nucstathyperinvconc,0.03,10);
    /*
    nucstatacc1 += NucRatesHyperScalingMove(nucstathyperinvconc,1.0,10);
    nucstatacc2 += NucRatesHyperScalingMove(nucstathyperinvconc,0.3,10);
    nucstatacc3 += NucRatesHyperScalingMove(nucstathyperinvconc,0.03,10);
    */
    nucstattot1 ++;
    nucstattot2 ++;
    nucstattot3 ++;

    // }

    nucrelratearray->SetConcentration(1.0/nucrelratehyperinvconc);
    nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
}

void MultiGeneCodonM2aModel::MasterMoveMixtureHyperParameters()  {

    if (purommode == 1) {
        MixtureHyperSlidingMove(puromhypermean,1.0,10,0,1);
        MixtureHyperSlidingMove(puromhypermean,0.3,10,0,1);
        MixtureHyperScalingMove(puromhyperinvconc,1.0,10);
        MixtureHyperScalingMove(puromhyperinvconc,0.3,10);
    }

    if (dposommode == 1)    {
        MixtureHyperSlidingMove(dposomhypermean,3.0,10,0.5,10);
        MixtureHyperSlidingMove(dposomhypermean,1.0,10,0.5,10);
        MixtureHyperSlidingMove(dposomhypermean,0.3,10,0.5,10);
        MixtureHyperScalingMove(dposomhyperinvshape,1.0,10);
        MixtureHyperScalingMove(dposomhyperinvshape,0.3,10);
    }

    if (poswmode == 1)  {
        MixtureHyperSlidingMove(poswhypermean,1.0,10,0,1);
        MixtureHyperSlidingMove(poswhypermean,0.3,10,0,1);
        MixtureHyperScalingMove(poswhyperinvconc,1.0,10);
        MixtureHyperScalingMove(poswhyperinvconc,0.3,10);
    }

    if (purwmode == 1)  {
        MixtureHyperSlidingMove(purwhypermean,1.0,10,0,1);
        MixtureHyperSlidingMove(purwhypermean,0.3,10,0,1);
        MixtureHyperScalingMove(purwhyperinvconc,1.0,10);
        MixtureHyperScalingMove(purwhyperinvconc,0.3,10);
    }

    if (burnin > 10)    {
        if (pihyperinvconc)    {
            ResamplePi();
        }
    }
}

void MultiGeneCodonM2aModel::ResamplePi()   {

    int n0 = poswsuffstat.GetN0();
    int n1 = poswsuffstat.GetN1();
    if ((n0+n1) != Ngene)   {
        cerr << "error in resample pi\n";
        exit(1);
    }
    double pialpha = pihypermean / pihyperinvconc;
    double pibeta = (1-pihypermean) / pihyperinvconc;
    double a0 = Random::sGamma(pialpha + n0);
    double a1 = Random::sGamma(pibeta + n1);
    pi = a1 / (a0 + a1);
}

double MultiGeneCodonM2aModel::MixtureHyperSlidingMove(double& x, double tuning, int nrep, double min, double max)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - MixtureHyperLogPrior() - MixtureHyperSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        x += m;
        if (max > min)  {
            while ((x < min) || (x > max))  {
                if (x < min)    {
                    x = 2*min - x;
                }
                if (x > max)    {
                    x = 2*max - x;
                }
            }
        }
        deltalogprob += MixtureHyperLogPrior() + MixtureHyperSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            x -= m;
        }
        ntot++;
    }
    return nacc/ntot;
}

double MultiGeneCodonM2aModel::MixtureHyperScalingMove(double& x, double tuning, int nrep)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double deltalogprob = - MixtureHyperLogPrior() - MixtureHyperSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        x *= e;
        deltalogprob += MixtureHyperLogPrior() + MixtureHyperSuffStatLogProb();
        deltalogprob += m;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            x /= e;
        }
        ntot++;
    }
    return nacc/ntot;
}

void MultiGeneCodonM2aModel::MasterMoveNucRates()    {

    MoveRR(0.1,1,3);
    MoveRR(0.03,3,3);
    MoveRR(0.01,3,3);

    MoveNucStat(0.1,1,3);
    MoveNucStat(0.01,1,3);
}

double MultiGeneCodonM2aModel::MoveRR(double tuning, int n, int nrep)	{

    vector<double>& nucrelrate = (*nucrelratearray)[0];

    double nacc = 0;
    double ntot = 0;
    double bk[Nrr];
    for (int rep=0; rep<nrep; rep++)	{
        for (int l=0; l<Nrr; l++)	{
            bk[l] = nucrelrate[l];
        }
        double deltalogprob = -GlobalNucRatesLogPrior() -NucRatesSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(nucrelrate,Nrr,tuning,n);
        deltalogprob += loghastings;
        UpdateNucMatrix();
        deltalogprob += GlobalNucRatesLogPrior() + NucRatesSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            for (int l=0; l<Nrr; l++)	{
                nucrelrate[l] = bk[l];
            }
            UpdateNucMatrix();
        }
        ntot++;
    }
    return nacc/ntot;
}

double MultiGeneCodonM2aModel::MoveNucStat(double tuning, int n, int nrep)	{

    vector<double>& nucstat = (*nucstatarray)[0];

    double nacc = 0;
    double ntot = 0;
    double bk[Nnuc];
    for (int rep=0; rep<nrep; rep++)	{
        for (int l=0; l<Nnuc; l++)	{
            bk[l] = nucstat[l];
        }
        double deltalogprob = -GlobalNucRatesLogPrior() -NucRatesSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(nucstat,Nnuc,tuning,n);
        deltalogprob += loghastings;
        UpdateNucMatrix();
        deltalogprob += GlobalNucRatesLogPrior() + NucRatesSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            for (int l=0; l<Nnuc; l++)	{
                nucstat[l] = bk[l];
            }
            UpdateNucMatrix();
        }
        ntot++;
    }
    return nacc/ntot;
}

//-------------------
// MPI Send Receive
// ------------------

void MultiGeneCodonM2aModel::MasterSendGlobalBranchLengths() {

    MasterSendGlobal(*branchlength);
}

void MultiGeneCodonM2aModel::SlaveReceiveGlobalBranchLengths()   {

    SlaveReceiveGlobal(*branchlength);
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetBranchLengths(*branchlength);
    }
}

void MultiGeneCodonM2aModel::MasterSendBranchLengthsHyperParameters() {

    MasterSendGlobal(*branchlength,blhyperinvshape);
}

void MultiGeneCodonM2aModel::SlaveReceiveBranchLengthsHyperParameters()   {

    SlaveReceiveGlobal(*branchlength,blhyperinvshape);
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength,blhyperinvshape);
    }
}

void MultiGeneCodonM2aModel::MasterSendGeneBranchLengths()    {

    MasterSendGeneArray(*branchlengtharray);
}

void MultiGeneCodonM2aModel::SlaveReceiveGeneBranchLengths()   {

    SlaveReceiveGeneArray(*branchlengtharray);
    for (int gene=0; gene<GetLocalNgene(); gene++)    {
        geneprocess[gene]->SetBranchLengths(branchlengtharray->GetVal(gene));
    }
}

void MultiGeneCodonM2aModel::SlaveSendGeneBranchLengths()    {

    SlaveSendGeneArray(*branchlengtharray);
}

void MultiGeneCodonM2aModel::MasterReceiveGeneBranchLengths()    {

    MasterReceiveGeneArray(*branchlengtharray);
}

void MultiGeneCodonM2aModel::MasterSendGlobalNucRates()   {

    MasterSendGlobal(nucrelratearray->GetVal(0),nucstatarray->GetVal(0));
}

void MultiGeneCodonM2aModel::SlaveReceiveGlobalNucRates()   {

    SlaveReceiveGlobal((*nucrelratearray)[0], (*nucstatarray)[0]);

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetNucRates((*nucrelratearray)[0],(*nucstatarray)[0]);
    }
}

void MultiGeneCodonM2aModel::MasterSendGeneNucRates()    {

    MasterSendGeneArray(*nucrelratearray,*nucstatarray);
}

void MultiGeneCodonM2aModel::SlaveReceiveGeneNucRates()   {

    SlaveReceiveGeneArray(*nucrelratearray,*nucstatarray);

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetNucRates((*nucrelratearray)[gene],(*nucstatarray)[gene]);
    }
}

void MultiGeneCodonM2aModel::SlaveSendGeneNucRates()    {

    SlaveSendGeneArray(*nucrelratearray,*nucstatarray);
}

void MultiGeneCodonM2aModel::MasterReceiveGeneNucRates()    {

    MasterReceiveGeneArray(*nucrelratearray,*nucstatarray);
}

void MultiGeneCodonM2aModel::MasterSendMixtureHyperParameters() {

    MasterSendGlobal(mixhyperparam);
}

void MultiGeneCodonM2aModel::SlaveReceiveMixtureHyperParameters()   {

    SlaveReceiveGlobal(mixhyperparam);

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetMixtureHyperParameters(puromhypermean,puromhyperinvconc,dposomhypermean,dposomhyperinvshape,pi,purwhypermean,purwhyperinvconc,poswhypermean,poswhyperinvconc);
    }
}

void MultiGeneCodonM2aModel::MasterSendNucRatesHyperParameters()   {

    MasterSendGlobal(nucrelratehypercenter,nucrelratehyperinvconc);
    MasterSendGlobal(nucstathypercenter,nucstathyperinvconc);
}

void MultiGeneCodonM2aModel::SlaveReceiveNucRatesHyperParameters()   {

    SlaveReceiveGlobal(nucrelratehypercenter,nucrelratehyperinvconc);
    SlaveReceiveGlobal(nucstathypercenter,nucstathyperinvconc);

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter,nucrelratehyperinvconc,nucstathypercenter,nucstathyperinvconc);
    }
}

void MultiGeneCodonM2aModel::SlaveSendBranchLengthsSuffStat()  {

    lengthsuffstatarray->Clear();
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->CollectLengthSuffStat();
        lengthsuffstatarray->Add(*geneprocess[gene]->GetLengthSuffStatArray());
    }
    SlaveSendAdditive(*lengthsuffstatarray);
}

void MultiGeneCodonM2aModel::MasterReceiveBranchLengthsSuffStat()  {

    lengthsuffstatarray->Clear();
    MasterReceiveAdditive(*lengthsuffstatarray);
}

void MultiGeneCodonM2aModel::SlaveSendBranchLengthsHyperSuffStat()   {

    lengthhypersuffstatarray->Clear();
    branchlengtharray->AddSuffStat(*lengthhypersuffstatarray);
    SlaveSendAdditive(*lengthhypersuffstatarray);
}

void MultiGeneCodonM2aModel::MasterReceiveBranchLengthsHyperSuffStat()   {

    lengthhypersuffstatarray->Clear();
    MasterReceiveAdditive(*lengthhypersuffstatarray);
}

void MultiGeneCodonM2aModel::SlaveSendNucRatesHyperSuffStat()   {

    nucrelratesuffstat.Clear();
    nucrelratearray->AddSuffStat(nucrelratesuffstat);
    SlaveSendAdditive(nucrelratesuffstat);

    nucstatsuffstat.Clear();
    nucstatarray->AddSuffStat(nucstatsuffstat);
    SlaveSendAdditive(nucstatsuffstat);
}

void MultiGeneCodonM2aModel::MasterReceiveNucRatesHyperSuffStat()   {

    nucrelratesuffstat.Clear();
    MasterReceiveAdditive(nucrelratesuffstat);

    nucstatsuffstat.Clear();
    MasterReceiveAdditive(nucstatsuffstat);
}

void MultiGeneCodonM2aModel::SlaveSendMixtureHyperSuffStat()  {

    puromsuffstat.Clear();
    puromarray->AddSuffStat(puromsuffstat);
    SlaveSendAdditive(puromsuffstat);

    dposomsuffstat.Clear();
    dposomarray->AddSuffStat(dposomsuffstat,*poswarray);
    SlaveSendAdditive(dposomsuffstat);

    purwsuffstat.Clear();
    purwarray->AddSuffStat(purwsuffstat);
    SlaveSendAdditive(purwsuffstat);

    poswsuffstat.Clear();
    poswarray->AddSuffStat(poswsuffstat);
    SlaveSendAdditive(poswsuffstat);
}

void MultiGeneCodonM2aModel::MasterReceiveMixtureHyperSuffStat()  {

    puromsuffstat.Clear();
    MasterReceiveAdditive(puromsuffstat);
    dposomsuffstat.Clear();
    MasterReceiveAdditive(dposomsuffstat);
    purwsuffstat.Clear();
    MasterReceiveAdditive(purwsuffstat);
    poswsuffstat.Clear();
    MasterReceiveAdditive(poswsuffstat);
}

void MultiGeneCodonM2aModel::SlaveCollectPathSuffStat() {
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->CollectPathSuffStat();
    }
}

void MultiGeneCodonM2aModel::SlaveSendNucPathSuffStat()  {

    nucpathsuffstat.Clear();
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->CollectComponentPathSuffStat();
        geneprocess[gene]->CollectNucPathSuffStat();
        nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
    }

    SlaveSendAdditive(nucpathsuffstat);
}

void MultiGeneCodonM2aModel::MasterReceiveNucPathSuffStat()  {

    nucpathsuffstat.Clear();
    MasterReceiveAdditive(nucpathsuffstat);
}

void MultiGeneCodonM2aModel::SlaveSendMixture()   {

    SlaveSendGeneArray(*puromarray,*dposomarray);
    SlaveSendGeneArray(*purwarray,*poswarray);
}

void MultiGeneCodonM2aModel::MasterReceiveMixture()    {

    MasterReceiveGeneArray(*puromarray,*dposomarray);
    MasterReceiveGeneArray(*purwarray,*poswarray);
}

void MultiGeneCodonM2aModel::MasterSendMixture()    {

    MasterSendGeneArray(*puromarray,*dposomarray);
    MasterSendGeneArray(*purwarray,*poswarray);
}

void MultiGeneCodonM2aModel::SlaveReceiveMixture()   {

    SlaveReceiveGeneArray(*puromarray,*dposomarray);
    SlaveReceiveGeneArray(*purwarray,*poswarray);

    for (int gene=0; gene<GetLocalNgene(); gene++)    {
        geneprocess[gene]->SetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
    }
}

void MultiGeneCodonM2aModel::SlaveSendLogProbs()   {

    GeneLogPrior = 0;
    lnL = 0;
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        GeneLogPrior += geneprocess[gene]->GetLogPrior();
        lnL += geneprocess[gene]->GetLogLikelihood();
    }
    SlaveSendAdditive(GeneLogPrior);
    SlaveSendAdditive(lnL);
}

void MultiGeneCodonM2aModel::MasterReceiveLogProbs()    {

    GeneLogPrior = 0;
    MasterReceiveAdditive(GeneLogPrior);
    lnL = 0;
    MasterReceiveAdditive(lnL);
}

/*
void MultiGeneCodonM2aModel::MasterTraceSitesPostProb(ostream& os)  {
    vector<double> array(GetTotNsite());
    MasterReceiveGene(array);
    int i = 0;
    for (int gene=0; gene<Ngene; gene++)    {
            os << GeneName[gene] << '\t';
            int nsite = GeneNsite[gene];
            for (int k=0; k<nsite; k++) {
                if (array[i] < 0)   {
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
    if (i != totnsite)  {
        cerr << "error in MultiGeneCodonM2aModel::MasterTraceSitesPostProb: non matching number of sites\n";
        exit(1);
    }
}

void MultiGeneCodonM2aModel::SlaveTraceSitesPostProb()  {

}
*/

void MultiGeneCodonM2aModel::MasterTraceSitesPostProb(ostream& os)  {

    for (int proc=1; proc<GetNprocs(); proc++)  {
        int totnsite = GetSlaveTotNsite(proc);
        double* array = new double[totnsite];
        MPI_Status stat;
        MPI_Recv(array,totnsite,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

        int i = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (GeneAlloc[gene] == proc)    {
                os << GeneName[gene] << '\t';
                int nsite = GeneNsite[gene];
                for (int k=0; k<nsite; k++) {
                    if (array[i] < 0)   {
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
        if (i != totnsite)  {
            cerr << "error in MultiGeneCodonM2aModel::MasterTraceSitesPostProb: non matching number of sites\n";
            exit(1);
        }
    }
    os << '\n';
    os.flush();
}

void MultiGeneCodonM2aModel::SlaveTraceSitesPostProb()  {

    int ngene = GetLocalNgene();
    int totnsite = GetLocalTotNsite();
    double* array = new double[totnsite];
    int i = 0;
    for (int gene=0; gene<ngene; gene++)    {
        geneprocess[gene]->GetSitesPostProb(array+i);
        for (int j=0; j<GeneNsite[gene]; j++)   {
            if (array[i+j] < 0) {
                cerr << "error in slave\n";
                cerr << i << '\t' << j << '\t' << GeneName[gene] << '\t' << GeneNsite[gene] << '\t' << geneprocess[gene]->GetNsite() << '\n';
                exit(1);
            }
        }
        i += GetLocalGeneNsite(gene);
    }
    if (i != totnsite)  {
        cerr << "error in MultiGeneCodonM2aModel::SlaveTraceSitesPostProb: non matching number of sites\n";
        exit(1);
    }

    MPI_Send(array,totnsite,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] array;
}
