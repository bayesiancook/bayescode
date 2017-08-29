
#include "MultiGeneCodonM2aModel.hpp"


//-------------------
// Constructing and setting the model
// ------------------

MultiGeneCodonM2aModel::MultiGeneCodonM2aModel(string datafile, string intreefile, double inpihypermean, double inpihyperinvconc, int inmyid, int innprocs) : MultiGeneMPIModule(inmyid,innprocs), nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {

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

    timepercycle.Reset();
    mastersampling.Reset();
    // Allocate();
}

void MultiGeneCodonM2aModel::Allocate() {

    lambda = 10;
    branchlength = new BranchIIDGamma(*tree,1.0,lambda);
    blhyperinvshape = 1.0;
    if (blmode == 2)    {
        lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        lengthhypersuffstatarray = 0;
    }
    else    {
        branchlengtharray.assign(GetLocalNgene(),(GammaWhiteNoise*) 0);
        for (int gene=0; gene<GetLocalNgene(); gene++)  {
            branchlengtharray[gene] = new GammaWhiteNoise(*tree,*branchlength,1.0/blhyperinvshape);
        }
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

    if (myid)   {
        lnL = new double[GetLocalNgene()];
    }
    else    {
        lnL = 0;
    }

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
        MasterReceiveLogLikelihood();
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

        SlaveSendLogLikelihood();
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

/*
void MultiGeneCodonM2aModel::SetBranchLengthsHyperParameters()  {
}

void MultiGeneCodonM2aModel::SetNucRateHyperParameters() {
}
*/

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

void MultiGeneCodonM2aModel::SlaveSetMixtureArrays()    {

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

    os << "#time";
    os << "\tlogprior\tlnL";
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
        os << "\tstdevrr\thyperinvconc";
        os << "\tstdevstat\thyperinvconc";
    }
    os << '\n';
    timepercycle.Start();
}

void MultiGeneCodonM2aModel::Trace(ostream& os)    {
    timepercycle.Stop();
    os << timepercycle.GetTime();
    // cerr << myid << '\t' << timepercycle.GetTime() << '\t' << mastersampling.GetTime() << '\t' << omegachrono.GetTime() << '\t' << hyperchrono.GetTime() << '\n';
    mastersampling.Reset();
    omegachrono.Reset();
    hyperchrono.Reset();
    timepercycle.Reset();
    timepercycle.Start();
    os << '\t' << GetLogPrior();
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
        os << '\t' << sqrt(GetVarNucRelRate()) << '\t' << nucrelratehyperinvconc;
        os << '\t' << sqrt(GetVarNucStat()) << '\t' << nucstathyperinvconc;
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

    double tot = 0;
    for (int j=0; j<Nbranch; j++)   {
        double mean = 0;
        for (int g=0; g<Ngene; g++) {
            double tmp = branchlengtharray[g]->GetVal(j);
            mean += tmp;
        }
        mean /= Ngene;
        tot += mean;
    }
    return tot;
}

double MultiGeneCodonM2aModel::GetVarLength()   {

    if (blmode == 2)    {
        cerr << "error: in getvarlength\n";
        exit(1);
    }

    double tot = 0;
    for (int j=0; j<Nbranch; j++)   {
        double mean = 0;
        double var = 0;
        for (int g=0; g<Ngene; g++) {
            double tmp = branchlengtharray[g]->GetVal(j);
            mean += tmp;
            var += tmp*tmp;
        }
        mean /= Ngene;
        var /= Ngene;
        var -= mean*mean;
        tot += var;
    }
    tot /= Nbranch;
    return tot;
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

    double total = 0;

    // branch lengths
    total += BranchLengthsHyperLogPrior();
    total += BranchLengthsLogPrior();

    // nuc rates
    total += NucRatesHyperLogPrior();
    total += NucRatesLogPrior();

    // mixture
    total += MixtureHyperLogPrior();
    total += MixtureLogPrior();

    return total;
}

double MultiGeneCodonM2aModel::LambdaHyperLogPrior()    {
    return -lambda/10;
}

double MultiGeneCodonM2aModel::BranchLengthsHyperInvShapeLogPrior()    {

    return -blhyperinvshape;
}

double MultiGeneCodonM2aModel::BranchLengthsHyperLogPrior()	{

    // exponential prior on scale parameter
    if (blmode == 0)    {
        return 0;
    }
    double total = LambdaHyperLogPrior();
    if (blmode == 1)    {
        // exponential hyper prior on inverse shape parameter
        total += BranchLengthsHyperInvShapeLogPrior();
        // iid gamma over branch-specific means across genes
        total += branchlength->GetLogProb();
    }
    return total;
}

double MultiGeneCodonM2aModel::BranchLengthsLogPrior()	{

    double total = 0;
    if (blmode == 2)    {
        total += branchlength->GetLogProb();
    }
    else    {
        for (int gene=0; gene<GetLocalNgene(); gene++)  {
            total += branchlengtharray[gene]->GetLogProb();
        }
    }
    return total;
}

double MultiGeneCodonM2aModel::NucRatesHyperLogPrior()  {

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

double MultiGeneCodonM2aModel::NucRatesLogPrior()    {

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

double MultiGeneCodonM2aModel::MixtureLogPrior()    {

    double total = 0;
    total += puromarray->GetLogProb();
    total += dposomarray->GetLogProb();
    total += purwarray->GetLogProb();
    total += poswarray->GetLogProb();
    return total;
}

double MultiGeneCodonM2aModel::GetLogLikelihood()   {
    return totlnL;
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

        if (rep)    {
            mastersampling.Start();
        }

        if (blmode == 2)    {
            MasterReceiveBranchLengthsSuffStat();
            MasterResampleBranchLengths();
            MasterMoveLambda();
            MasterSendGlobalBranchLengths();
        }
        else    {
            if (blmode == 1)    {
                MasterReceiveBranchLengthsHyperSuffStat();
                MasterMoveBranchLengthsHyperParameters();
                MasterSendBranchLengthsHyperParameters();
            }
        }

        omegachrono.Start();
        MasterReceiveMixtureHyperSuffStat();
        omegachrono.Stop();
        hyperchrono.Start();
        MasterMoveMixtureHyperParameters();
        hyperchrono.Stop();
        MasterSendMixtureHyperParameters();

        if (nucmode == 2)   {
            MasterReceiveNucPathSuffStat();
            MasterMoveNucRates();
            MasterSendGlobalNucRates();
        }
        else    {
            if (nucmode == 1)  {
                MasterReceiveNucRatesHyperSuffStat();
                MasterMoveNucRatesHyperParameters();
                MasterSendNucRatesHyperParameters();
            }
        }

        if (rep)    {
            mastersampling.Stop();
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
    MasterReceiveLogLikelihood();
}

// slave move
void MultiGeneCodonM2aModel::SlaveMove() {

    Chrono mapping, sampling;

    mapping.Start();
    SlaveResampleSub();
    mapping.Stop();

    int nrep = 30;

    sampling.Start();
    for (int rep=0; rep<nrep; rep++)	{

        if (blmode == 2)    {
            SlaveSendBranchLengthsSuffStat();
            SlaveReceiveGlobalBranchLengths();
        }
        else    {
            SlaveMoveBranchLengths();

            if (blmode == 1)    {
                SlaveSendBranchLengthsHyperSuffStat();
                SlaveReceiveBranchLengthsHyperParameters();
            }
        }

        SlaveCollectPathSuffStat();
        SlaveMoveOmega();
        SlaveSendMixtureHyperSuffStat();
        SlaveReceiveMixtureHyperParameters();
        SlaveSetMixtureArrays();

        if (nucmode == 2)   {
            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalNucRates();
        }
        else    {
            SlaveMoveNucRates();

            if (nucmode == 1)   {
                SlaveSendNucRatesHyperSuffStat();
                SlaveReceiveNucRatesHyperParameters();
            }
        }
    }
    sampling.Stop();
    burnin++;

    if (blmode != 2)    {
        SlaveSendGeneBranchLengths();
    }
    if (nucmode != 2)   {
        SlaveSendGeneNucRates();
    }
    SlaveSendMixture();
    SlaveSendLogLikelihood();
}

void MultiGeneCodonM2aModel::SlaveResampleSub()  {

    double frac = 1.0;
    /*
    if (burnin > 10)    {
        frac = 0.2;
    }
    */
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->ResampleSub(frac);
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

void MultiGeneCodonM2aModel::SlaveMoveBranchLengths() {

    for (int gene=0; gene<GetLocalNgene(); gene++)  {
        geneprocess[gene]->MoveBranchLengths();
        geneprocess[gene]->GetBranchLengths(*branchlengtharray[gene]);
    }
}

double MultiGeneCodonM2aModel::NucRatesHyperProfileMove(vector<double>& x, double tuning, int n, int nrep)	{

    double nacc = 0;
    double ntot = 0;
    vector<double> bk(x.size(),0);
    for (int rep=0; rep<nrep; rep++)	{
        bk = x;
        double deltalogprob = - NucRatesHyperLogPrior() - NucRatesHyperSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(x,x.size(),tuning,n);
        deltalogprob += NucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
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
        double deltalogprob = - NucRatesHyperLogPrior() - NucRatesHyperSuffStatLogProb();
        double m = tuning * (Random::Uniform() - 0.5);
        double e = exp(m);
        x *= e;
        deltalogprob += NucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
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

    NucRatesHyperProfileMove(nucrelratehypercenter,1.0,1,10);
    NucRatesHyperProfileMove(nucrelratehypercenter,0.3,1,10);
    NucRatesHyperProfileMove(nucrelratehypercenter,0.1,3,10);
    NucRatesHyperScalingMove(nucrelratehyperinvconc,1.0,10);
    NucRatesHyperScalingMove(nucrelratehyperinvconc,0.3,10);

    NucRatesHyperProfileMove(nucstathypercenter,1.0,1,10);
    NucRatesHyperProfileMove(nucstathypercenter,0.3,1,10);
    NucRatesHyperProfileMove(nucstathypercenter,0.1,2,10);
    NucRatesHyperScalingMove(nucstathyperinvconc,1.0,10);
    NucRatesHyperScalingMove(nucstathyperinvconc,0.3,10);
}

void MultiGeneCodonM2aModel::SlaveMoveNucRates()    {

    for (int gene=0; gene<GetLocalNgene(); gene++)  {
        geneprocess[gene]->MoveNucRates();
        geneprocess[gene]->GetNucRates((*nucrelratearray)[gene],(*nucstatarray)[gene]);
    }
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
        double deltalogprob = -NucRatesLogPrior() -NucRatesSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(nucrelrate,Nrr,tuning,n);
        deltalogprob += loghastings;
        UpdateNucMatrix();
        deltalogprob += NucRatesLogPrior() + NucRatesSuffStatLogProb();
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
        double deltalogprob = -NucRatesLogPrior() -NucRatesSuffStatLogProb();
        double loghastings = Random::ProfileProposeMove(nucstat,Nnuc,tuning,n);
        deltalogprob += loghastings;
        UpdateNucMatrix();
        deltalogprob += NucRatesLogPrior() + NucRatesSuffStatLogProb();
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

void MultiGeneCodonM2aModel::SlaveMoveOmega()  {
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->MoveOmega();
        geneprocess[gene]->GetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
    }
}

//-------------------
// MPI Send Receive
// ------------------

void MultiGeneCodonM2aModel::MasterSendGlobalBranchLengths() {

    double* array = new double[Nbranch];
    for (int j=0; j<Nbranch; j++)   {
        array[j] = branchlength->GetVal(j);
    }
    MPI_Bcast(array,Nbranch,MPI_DOUBLE,0,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveReceiveGlobalBranchLengths()   {

    double* array = new double[Nbranch];
    MPI_Bcast(array,Nbranch,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (int j=0; j<Nbranch; j++)   {
        (*branchlength)[j] = array[j];
    }
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetBranchLengths(*branchlength);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterSendBranchLengthsHyperParameters() {

    double* array = new double[Nbranch+1];
    for (int j=0; j<Nbranch; j++)   {
        array[j] = branchlength->GetVal(j);
    }
    array[Nbranch] = blhyperinvshape;

    MPI_Bcast(array,Nbranch+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveReceiveBranchLengthsHyperParameters()   {

    double* array = new double[Nbranch+1];
    MPI_Bcast(array,Nbranch+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (int j=0; j<Nbranch; j++)   {
        (*branchlength)[j] = array[j];
    }
    blhyperinvshape = array[Nbranch];

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength,blhyperinvshape);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterSendGeneBranchLengths()    {

    for (int proc=1; proc<GetNprocs(); proc++)  {
        int ngene = GetSlaveNgene(proc);
        double* array = new double[Nbranch*ngene];
        int index = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (GeneAlloc[gene] == proc)    {
                for (int j=0; j<Nbranch; j++)   {
                    array[index++] = branchlengtharray[gene]->GetVal(j);
                }
            }
        }
        if (index != Nbranch*ngene) {
            cerr << "error when sending gene branch lengths (master): non matching vector size\n";
            exit(1);
        }
        // cerr << "master: " << Nbranch << '\t' << ngene << '\t' << Nbranch*ngene << '\t' << index << '\n';

        MPI_Send(array,Nbranch*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }
}

void MultiGeneCodonM2aModel::SlaveReceiveGeneBranchLengths()   {

    int ngene = GetLocalNgene();
    double* array = new double[Nbranch*ngene];
    MPI_Status stat;
    MPI_Recv(array,Nbranch*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

    int index = 0;
    for (int gene=0; gene<ngene; gene++)    {
        for (int j=0; j<Nbranch; j++)   {
            (*branchlengtharray[gene])[j] = array[index++];
        }
        geneprocess[gene]->SetBranchLengths(*branchlengtharray[gene]);
    }
    if (index != Nbranch*ngene) {
        cerr << "error when receiving gene branch lengths (slave): non matching vector size\n";
        cerr << Nbranch << '\t' << ngene << '\t' << Nbranch*ngene << '\t' << index << '\n';
        exit(1);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveSendGeneBranchLengths()    {

    int ngene = GetLocalNgene();
    double* array = new double[Nbranch*ngene];

    int index = 0;
    for (int gene=0; gene<ngene; gene++)    {
        // geneprocess[gene]->GetBranchLengths(*branchlengtharray[gene]);
        for (int j=0; j<Nbranch; j++)   {
            array[index++] = branchlengtharray[gene]->GetVal(j);
        }
    }
    if (index != Nbranch*ngene) {
        cerr << "error when sending gene branch lengths (slave): non matching vector size\n";
        cerr << Nbranch << '\t' << ngene << '\t' << Nbranch*ngene << '\t' << index << '\n';
        exit(1);
    }
    MPI_Send(array,Nbranch*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterReceiveGeneBranchLengths()    {

    for (int proc=1; proc<GetNprocs(); proc++)  {
        int ngene = GetSlaveNgene(proc);
        double* array = new double[Nbranch*ngene];

        MPI_Status stat;
        MPI_Recv(array,Nbranch*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

        int index = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (GeneAlloc[gene] == proc)    {
                for (int j=0; j<Nbranch; j++)   {
                    (*branchlengtharray[gene])[j] = array[index++];
                }
            }
        }
        if (index != Nbranch*ngene) {
            cerr << "error when receiving gene branch lengths (master): non matching vector size\n";
            exit(1);
        }
        delete[] array;
    }
}

void MultiGeneCodonM2aModel::MasterSendGlobalNucRates()   {

    int N = Nrr + Nnuc;
    double* array = new double[N];
    int i = 0;
    for (int j=0; j<Nrr; j++)   {
        array[i++] = (*nucrelratearray)[0][j];
    }
    for (int j=0; j<Nnuc; j++)  {
        array[i++] = (*nucstatarray)[0][j];
    }
    if (i != N) {
        cerr << "error when sending global nuc rates: non matching vector size\n";
        cerr << i << '\t' << N << '\n';
        exit(1);
    }

    MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveReceiveGlobalNucRates()   {

    int N = Nrr + Nnuc;
    double* array = new double[N];
    MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    int i = 0;
    for (int j=0; j<Nrr; j++)   {
        (*nucrelratearray)[0][j] = array[i++];
    }
    for (int j=0; j<Nnuc; j++)  {
        (*nucstatarray)[0][j] = array[i++];
    }

    if (i != N) {
        cerr << "error when receiving global nuc rates: non matching vector size\n";
        cerr << i << '\t' << N << '\n';
        exit(1);
    }

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetNucRates((*nucrelratearray)[0],(*nucstatarray)[0]);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterSendGeneNucRates()    {

    int N = Nrr + Nnuc;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        int ngene = GetSlaveNgene(proc);
        double* array = new double[N*ngene];
        int index = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (GeneAlloc[gene] == proc)    {
                for (int j=0; j<Nrr; j++)   {
                    array[index++] = (*nucrelratearray)[gene][j];
                }
                for (int j=0; j<Nnuc; j++)  {
                    array[index++] = (*nucstatarray)[gene][j];
                }
            }
        }
        if (index != N*ngene) {
            cerr << "error when sending gene nuc rates (master): non matching vector size\n";
            exit(1);
        }

        MPI_Send(array,N*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }
}

void MultiGeneCodonM2aModel::SlaveReceiveGeneNucRates()   {

    int N = Nrr + Nnuc;
    int ngene = GetLocalNgene();
    double* array = new double[N*ngene];
    MPI_Status stat;
    MPI_Recv(array,N*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

    int index = 0;
    for (int gene=0; gene<ngene; gene++)    {
        for (int j=0; j<Nrr; j++)   {
            (*nucrelratearray)[gene][j] = array[index++];
        }
        for (int j=0; j<Nnuc; j++)  {
            (*nucstatarray)[gene][j] = array[index++];
        }
        geneprocess[gene]->SetNucRates((*nucrelratearray)[gene],(*nucstatarray)[gene]);
    }
    if (index != N*ngene) {
        cerr << "error when receiving gene nuc rates (slave): non matching vector size\n";
        exit(1);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveSendGeneNucRates()    {

    int N = Nrr + Nnuc;
    int ngene = GetLocalNgene();
    double* array = new double[N*ngene];

    int index = 0;
    for (int gene=0; gene<ngene; gene++)    {
        // geneprocess[gene]->GetNucRates((*nucrelratearray)[gene],(*nucstatarray)[gene]);
        for (int j=0; j<Nrr; j++)   {
            array[index++] = (*nucrelratearray)[gene][j];
        }
        for (int j=0; j<Nnuc; j++)  {
            array[index++] = (*nucstatarray)[gene][j];
        }
    }
    if (index != N*ngene) {
        cerr << "error when sending gene nuc rates (slaves): non matching vector size\n";
        exit(1);
    }
    MPI_Send(array,N*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterReceiveGeneNucRates()    {

    int N = Nrr + Nnuc;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        int ngene = GetSlaveNgene(proc);
        double* array = new double[N*ngene];

        MPI_Status stat;
        MPI_Recv(array,N*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

        int index = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (GeneAlloc[gene] == proc)    {
                for (int j=0; j<Nrr; j++)   {
                    (*nucrelratearray)[gene][j] = array[index++];
                }
                for (int j=0; j<Nnuc; j++)  {
                    (*nucstatarray)[gene][j] = array[index++];
                }
            }
        }
        if (index != N*ngene) {
            cerr << "error when receiving gene nuc rates (master): non matching vector size\n";
            exit(1);
        }
        delete[] array;
    }
}

void MultiGeneCodonM2aModel::MasterSendMixtureHyperParameters() {

    int N = 9;
    double* array = new double[N];
    int i = 0;
    array[i++] = puromhypermean;
    array[i++] = puromhyperinvconc;
    array[i++] = dposomhypermean;
    array[i++] = dposomhyperinvshape;
    array[i++] = purwhypermean;
    array[i++] = purwhyperinvconc;
    array[i++] = pi;
    array[i++] = poswhypermean;
    array[i++] = poswhyperinvconc;

    if (i != N) {
        cerr << "error when sending global params: non matching vector size\n";
        cerr << i << '\t' << N << '\n';
        exit(1);
    }

    MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveReceiveMixtureHyperParameters()   {

    int N = 9;
    double* array = new double[N];
    MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    int i = 0;
    puromhypermean = array[i++];
    puromhyperinvconc = array[i++];
    dposomhypermean = array[i++];
    dposomhyperinvshape = array[i++];
    purwhypermean = array[i++];
    purwhyperinvconc = array[i++];
    pi = array[i++];
    poswhypermean = array[i++];
    poswhyperinvconc = array[i++];

    if (i != N) {
        cerr << "error when sending global params: non matching vector size\n";
        cerr << i << '\t' << N << '\n';
        exit(1);
    }

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetMixtureHyperParameters(puromhypermean,puromhyperinvconc,dposomhypermean,dposomhyperinvshape,pi,purwhypermean,purwhyperinvconc,poswhypermean,poswhyperinvconc);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterSendNucRatesHyperParameters()   {

    int N = Nrr + Nnuc + 2;
    double* array = new double[N];
    int i = 0;
    for (int j=0; j<Nrr; j++)   {
        array[i++] = nucrelratehypercenter[j];
    }
    array[i++] = nucrelratehyperinvconc;
    for (int j=0; j<Nnuc; j++)  {
        array[i++] = nucstathypercenter[j];
    }
    array[i++] = nucstathyperinvconc;
    if (i != N) {
        cerr << "error when sending global params: non matching vector size\n";
        cerr << i << '\t' << N << '\n';
        exit(1);
    }

    MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveReceiveNucRatesHyperParameters()   {

    int N = Nrr + Nnuc + 2;
    double* array = new double[N];
    MPI_Bcast(array,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    int i = 0;
    for (int j=0; j<Nrr; j++)   {
        nucrelratehypercenter[j] = array[i++];
    }
    nucrelratehyperinvconc = array[i++];
    for (int j=0; j<Nnuc; j++)  {
        nucstathypercenter[j] = array[i++];
    }
    nucstathyperinvconc = array[i++];

    if (i != N) {
        cerr << "error when sending global params: non matching vector size\n";
        cerr << i << '\t' << N << '\n';
        exit(1);
    }

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter,nucrelratehyperinvconc,nucstathypercenter,nucstathyperinvconc);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::SlaveSendBranchLengthsSuffStat()  {

    lengthsuffstatarray->Clear();
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->CollectLengthSuffStat();
        lengthsuffstatarray->Add(*geneprocess[gene]->GetLengthSuffStatArray());
    }

    int* count = new int[Nbranch];
    double* beta = new double[Nbranch];
    lengthsuffstatarray->Push(count,beta);
    MPI_Send(count,Nbranch,MPI_INT,0,TAG1,MPI_COMM_WORLD);
    MPI_Send(beta,Nbranch,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::MasterReceiveBranchLengthsSuffStat()  {

    int* count = new int[Nbranch];
    double* beta = new double[Nbranch];
    lengthsuffstatarray->Clear();
    MPI_Status stat;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        MPI_Recv(count,Nbranch,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
        MPI_Recv(beta,Nbranch,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
        lengthsuffstatarray->Add(count,beta);
    }
    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::SlaveSendBranchLengthsHyperSuffStat()   {

    int* count = new int[Nbranch];
    double* beta = new double[2*Nbranch];

    lengthhypersuffstatarray->Clear();
    for (int gene=0; gene<GetLocalNgene(); gene++)  {
        branchlengtharray[gene]->AddSuffStat(*lengthhypersuffstatarray);
    }
    for (int j=0; j<Nbranch; j++)   {
        count[j] = lengthhypersuffstatarray->GetVal(j).GetN();
        beta[j] = lengthhypersuffstatarray->GetVal(j).GetSum();
        beta[Nbranch+j] = lengthhypersuffstatarray->GetVal(j).GetSumLog();
    }

    MPI_Send(count,Nbranch,MPI_INT,0,TAG1,MPI_COMM_WORLD);
    MPI_Send(beta,2*Nbranch,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::MasterReceiveBranchLengthsHyperSuffStat()   {

    lengthhypersuffstatarray->Clear();

    int* count = new int[Nbranch];
    double* beta = new double[2*Nbranch];

    MPI_Status stat;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        MPI_Recv(count,Nbranch,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
        MPI_Recv(beta,2*Nbranch,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

        for (int j=0; j<Nbranch; j++)   {
            (*lengthhypersuffstatarray)[j].AddSuffStat(beta[j],beta[Nbranch+j],count[j]);
        }
    }

    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::SlaveSendNucRatesHyperSuffStat()   {

    int Nint = 2;
    int Ndouble = Nrr + Nnuc;
    int* count = new int[Nint];
    double* beta = new double[Ndouble];

    int i = 0;
    int d = 0;

    nucrelratesuffstat.Clear();
    nucrelratearray->AddSuffStat(nucrelratesuffstat);
    count[i++] = nucrelratesuffstat.GetN();
    for (int j=0; j<Nrr; j++)   {
        beta[d++] = nucrelratesuffstat.GetSumLog(j);
    }

    nucstatsuffstat.Clear();
    nucstatarray->AddSuffStat(nucstatsuffstat);
    count[i++] = nucstatsuffstat.GetN();
    for (int j=0; j<Nnuc; j++)   {
        beta[d++] = nucstatsuffstat.GetSumLog(j);
    }

    MPI_Send(count,Nint,MPI_INT,0,TAG1,MPI_COMM_WORLD);
    MPI_Send(beta,Ndouble,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::MasterReceiveNucRatesHyperSuffStat()   {

        nucrelratesuffstat.Clear();
        nucpathsuffstat.Clear();

        int Nint = 2;
        int Ndouble = Nrr + Nnuc;
        int* count = new int[Nint];
        double* beta = new double[Ndouble];

        MPI_Status stat;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPI_Recv(count,Nint,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
            MPI_Recv(beta,Ndouble,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

            int i = 0;
            int d = 0;

            nucrelratesuffstat.AddSuffStat(beta+d,count[i++]);
            d+=Nrr;

            nucstatsuffstat.AddSuffStat(beta+d,count[i++]);
            d+=Nnuc;
        }

        delete[] count;
        delete[] beta;
}

void MultiGeneCodonM2aModel::SlaveSendMixtureHyperSuffStat()  {

    int Nint = 1 + 1 + 1 + 2;
    int Ndouble = 2 + 2 + 2 + 2;
    int* count = new int[Nint];
    double* beta = new double[Ndouble];

    int i = 0;
    int d = 0;

    puromsuffstat.Clear();
    puromarray->AddSuffStat(puromsuffstat);
    count[i++] = puromsuffstat.GetN();
    beta[d++] = puromsuffstat.GetSumLog0();
    beta[d++] = puromsuffstat.GetSumLog1();

    dposomsuffstat.Clear();
    // dposomarray->AddSuffStat(dposomsuffstat);
    dposomarray->AddSuffStat(dposomsuffstat,*poswarray);
    count[i++] = dposomsuffstat.GetN();
    beta[d++] = dposomsuffstat.GetSum();
    beta[d++] = dposomsuffstat.GetSumLog();

    purwsuffstat.Clear();
    purwarray->AddSuffStat(purwsuffstat);
    count[i++] = purwsuffstat.GetN();
    beta[d++] = purwsuffstat.GetSumLog0();
    beta[d++] = purwsuffstat.GetSumLog1();

    poswsuffstat.Clear();
    poswarray->AddSuffStat(poswsuffstat);
    count[i++] = poswsuffstat.GetN0();
    count[i++] = poswsuffstat.GetN1();
    beta[d++] = poswsuffstat.GetSumLog0();
    beta[d++] = poswsuffstat.GetSumLog1();

    MPI_Send(count,Nint,MPI_INT,0,TAG1,MPI_COMM_WORLD);
    MPI_Send(beta,Ndouble,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::MasterReceiveMixtureHyperSuffStat()  {

    puromsuffstat.Clear();
    dposomsuffstat.Clear();
    purwsuffstat.Clear();
    poswsuffstat.Clear();

    int Nint = 1 + 1 + 1 + 2;
    int Ndouble = 2 + 2 + 2 + 2;
    int* count = new int[Nint];
    double* beta = new double[Ndouble];

    MPI_Status stat;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        MPI_Recv(count,Nint,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
        MPI_Recv(beta,Ndouble,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);

        int i = 0;
        int d = 0;

        puromsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
        d+=2;
        dposomsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
        d+=2;
        purwsuffstat.AddSuffStat(beta[d],beta[d+1],count[i++]);
        d+=2;
        poswsuffstat.AddNullSuffStat(count[i++]);
        poswsuffstat.AddPosSuffStat(beta[d],beta[d+1],count[i++]);
        d+=2;
    }

    delete[] count;
    delete[] beta;
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
        nucpathsuffstat.Add(geneprocess[gene]->GetNucPathSuffStat());
    }

    int* count = new int[Nnuc+Nnuc*Nnuc];
    double* beta = new double[Nnuc*Nnuc];
    nucpathsuffstat.Push(count,beta);
    MPI_Send(count,Nnuc+Nnuc*Nnuc,MPI_INT,0,TAG1,MPI_COMM_WORLD);
    MPI_Send(beta,Nnuc*Nnuc,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::MasterReceiveNucPathSuffStat()  {

    int* count = new int[Nnuc+Nnuc*Nnuc];
    double* beta = new double[Nnuc*Nnuc];
    nucpathsuffstat.Clear();
    MPI_Status stat;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        MPI_Recv(count,Nnuc+Nnuc*Nnuc,MPI_INT,proc,TAG1,MPI_COMM_WORLD,&stat);
        MPI_Recv(beta,Nnuc*Nnuc,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
        nucpathsuffstat.Add(count,beta);
    }
    delete[] count;
    delete[] beta;
}

void MultiGeneCodonM2aModel::SlaveSendMixture()   {

    int ngene = GetLocalNgene();
    double* array = new double[4*ngene];
    for (int gene=0; gene<ngene; gene++)    {
        array[gene] = geneprocess[gene]->GetPurOm();
        array[1*ngene+gene] = geneprocess[gene]->GetDPosOm();
        array[2*ngene+gene] = geneprocess[gene]->GetPurW();
        array[3*ngene+gene] = geneprocess[gene]->GetPosW();
    }
    MPI_Send(array,4*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterReceiveMixture()    {

    MPI_Status stat;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        int ngene = GetSlaveNgene(proc);
        double* array = new double[4*ngene];
        MPI_Recv(array,4*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
        int index = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (GeneAlloc[gene] == proc)    {
                (*puromarray)[gene] = array[index];
                (*dposomarray)[gene] = array[1*ngene+index];
                (*purwarray)[gene] = array[2*ngene+index];
                (*poswarray)[gene] = array[3*ngene+index];
                index++;
            }
        }
        delete[] array;
    }
}

void MultiGeneCodonM2aModel::SlaveReceiveMixture()   {

    int ngene = GetLocalNgene();
    double* array = new double[4*ngene];
    MPI_Status stat;
    MPI_Recv(array,4*ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

    for (int gene=0; gene<ngene; gene++)    {
        (*puromarray)[gene] = array[gene];
        (*dposomarray)[gene] = array[1*ngene+gene];
        (*purwarray)[gene] = array[2*ngene+gene];
        (*poswarray)[gene] = array[3*ngene+gene];
        geneprocess[gene]->SetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
    }
    delete[] array;
}

void MultiGeneCodonM2aModel::MasterSendMixture()    {

    for (int proc=1; proc<GetNprocs(); proc++)  {
        int ngene = GetSlaveNgene(proc);
        double* array = new double[4*ngene];
        int index = 0;
        for (int gene=0; gene<Ngene; gene++)    {
            if (GeneAlloc[gene] == proc)    {
                array[index] = (*puromarray)[gene];
                array[1*ngene+index] = (*dposomarray)[gene];
                array[2*ngene+index] = (*purwarray)[gene];
                array[3*ngene+index] = (*poswarray)[gene];
                index++;
            }
        }
        MPI_Send(array,4*ngene,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
        delete[] array;
    }
}

void MultiGeneCodonM2aModel::SlaveSendLogLikelihood()   {

    totlnL = 0;
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        lnL[gene] = geneprocess[gene]->GetLogLikelihood();
        totlnL += lnL[gene];
    }
    MPI_Send(&totlnL,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MultiGeneCodonM2aModel::MasterReceiveLogLikelihood()    {

    MPI_Status stat;
    totlnL = 0;
    for (int proc=1; proc<GetNprocs(); proc++)  {
        double tmp;
        MPI_Recv(&tmp,1,MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
        totlnL += tmp;
    }
}

/*
void MultiGeneCodonM2aModel::SlaveTracePostProbHeader(string name)    {
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        ofstream os((name + "_" + GetLocalGeneName(gene) + ".sitepp").c_str());
    }
}

void MultiGeneCodonM2aModel::SlaveTracePostProb(string name)    {
    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        ofstream os((name + "_" + GetLocalGeneName(gene) + ".sitepp").c_str());
        geneprocess[gene]->TracePostProb(os);
    }
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
