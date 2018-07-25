
#include "MultiGeneAAMutSelM2aModel.hpp"

//-------------------
// Constructing and setting the model
// ------------------

// MultiGeneCodonM2aModel::MultiGeneCodonM2aModel(string datafile, string
// intreefile, double inpihypermean, double inpihyperinvconc, int inmyid, int
// innprocs) :
MultiGeneAAMutSelM2aModel::MultiGeneAAMutSelM2aModel(string datafile, string intreefile,
                                                     double inpihypermean, double inpihyperinvconc,
                                                     int inmyid, int innprocs)
    :

      MultiGeneMPIModule(inmyid, innprocs),
      mixhyperparam(9, 0),
      // puromhypermean(mixhyperparam[0]),
      // puromhyperinvconc(mixhyperparam[1]),
      dposomhypermean(mixhyperparam[2]),
      dposomhyperinvshape(mixhyperparam[3]),
      // purwhypermean(mixhyperparam[4]),
      // purwhyperinvconc(mixhyperparam[5]),
      poswhypermean(mixhyperparam[6]),
      poswhyperinvconc(mixhyperparam[7]),
      pi(mixhyperparam[8]) {
    // nucrelratesuffstat(Nrr), nucstatsuffstat(Nnuc) {

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

    if (!myid) {
        std::cerr << "number of taxa : " << Ntaxa << '\n';
        std::cerr << "number of branches : " << Nbranch << '\n';
        std::cerr << "-- Tree and data fit together\n";
    }
}

void MultiGeneAAMutSelM2aModel::Allocate() {
    lambda = 10;
    branchlength = new BranchIIDGamma(*tree, 1.0, lambda);
    blhyperinvshape = 1.0;
    if (blmode == 2) {
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        lengthhypersuffstatarray = 0;
    } else {
        branchlength->SetAllBranches(1.0 / lambda);
        branchlengtharray =
            new GammaWhiteNoiseArray(GetLocalNgene(), *tree, *branchlength, 1.0 / blhyperinvshape);
        lengthpathsuffstatarray = 0;
        lengthhypersuffstatarray = new GammaSuffStatBranchArray(*tree);
    }

    // nucrelratehypercenter.assign(Nrr,1.0/Nrr);
    // nucrelratehyperinvconc = 1.0 / Nrr;

    // nucstathypercenter.assign(Nnuc,1.0/Nnuc);
    // nucstathyperinvconc = 1.0 / Nnuc;

    /*if (nucmode == 2)   {
        nucrelratearray = new
    IIDDirichlet(1,nucrelratehypercenter,1.0/nucrelratehyperinvconc); nucstatarray
    = new IIDDirichlet(1,nucstathypercenter,1.0/nucstathyperinvconc); nucmatrix =
    new GTRSubMatrix(Nnuc,(*nucrelratearray)[0],(*nucstatarray)[0],true);
    }
    else    {
        nucrelratearray = new
    IIDDirichlet(GetLocalNgene(),nucrelratehypercenter,1.0/nucrelratehyperinvconc);
        nucstatarray = new
    IIDDirichlet(GetLocalNgene(),nucstathypercenter,1.0/nucstathyperinvconc);
        nucmatrix = 0;
    }*/

    // double puromalpha = puromhypermean / puromhyperinvconc;
    // double purombeta = (1-puromhypermean) / puromhyperinvconc;
    // puromarray = new IIDBeta(GetLocalNgene(),puromalpha,purombeta);

    double dposomalpha = 1.0 / dposomhyperinvshape;
    double dposombeta = dposomalpha / dposomhypermean;
    dposomarray = new IIDGamma(GetLocalNgene(), dposomalpha, dposombeta);

    // double purwalpha = purwhypermean / purwhyperinvconc;
    // double purwbeta = (1-purwhypermean) / purwhyperinvconc;
    // purwarray = new IIDBeta(GetLocalNgene(),purwalpha,purwbeta);

    double poswalpha = poswhypermean / poswhyperinvconc;
    double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
    poswarray = new IIDBernoulliBeta(GetLocalNgene(), pi, poswalpha, poswbeta);

    if (!GetMyid()) {
        geneprocess.assign(0, (AAMutSelM2aModel *)0);
    } else {
        geneprocess.assign(GetLocalNgene(), (AAMutSelM2aModel *)0);

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene] = new AAMutSelM2aModel(GetLocalGeneName(gene), treefile, pi);
            // geneprocess[gene]->SetAcrossGenesModes(blmode,nucmode);
            geneprocess[gene]->SetAcrossGenesModes(blmode);
        }
    }
}

void MultiGeneAAMutSelM2aModel::Unfold() {
    if (!GetMyid()) {
        MasterSendBranchLengthsHyperParameters();
        // MasterSendNucRatesHyperParameters();
        MasterSendMixtureHyperParameters();

        if (blmode == 2) {
            MasterSendGlobalBranchLengths();
        } else {
            MasterSendGeneBranchLengths();
        }

        /*if (nucmode == 2)   {
            MasterSendGlobalNucRates();
        }
        else    {
            MasterSendGeneNucRates();
        }*/

        MasterSendMixture();
        MasterReceiveLogProbs();
    } else {
        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->Allocate();
        }

        SlaveReceiveBranchLengthsHyperParameters();
        // SlaveReceiveNucRatesHyperParameters();
        SlaveReceiveMixtureHyperParameters();

        if (blmode == 2) {
            SlaveReceiveGlobalBranchLengths();
        } else {
            SlaveReceiveGeneBranchLengths();
        }

        /*if (nucmode == 2)   {
            SlaveReceiveGlobalNucRates();
        }
        else    {
            SlaveReceiveGeneNucRates();
        }*/

        SlaveReceiveMixture();

        for (int gene = 0; gene < GetLocalNgene(); gene++) {
            geneprocess[gene]->UpdateMatrices();
            geneprocess[gene]->Unfold();
        }

        SlaveSendLogProbs();
    }
}

// void MultiGeneCodonM2aModel::SetAcrossGenesModes(int inblmode, int innucmode,
// int inpurommode, int indposommode, int inpurwmode, int inposwmode)  {
void MultiGeneAAMutSelM2aModel::SetAcrossGenesModes(int inblmode, int indposommode,
                                                    int inposwmode) {
    blmode = inblmode;
    // nucmode = innucmode;
    // purommode = inpurommode;
    dposommode = indposommode;
    // purwmode = inpurwmode;
    poswmode = inposwmode;
}

// void MultiGeneCodonM2aModel::SetMixtureHyperParameters(double
// inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double
// indposomhyperinvshape, double inpurwhypermean, double inpurwhyperinvconc,
// double inposwhypermean, double inposwhyperinvconc)   {
void MultiGeneAAMutSelM2aModel::SetMixtureHyperParameters(double indposomhypermean,
                                                          double indposomhyperinvshape,
                                                          double inposwhypermean,
                                                          double inposwhyperinvconc) {
    // puromhypermean = inpuromhypermean;
    // puromhyperinvconc = inpuromhyperinvconc;
    dposomhypermean = indposomhypermean;
    dposomhyperinvshape = indposomhyperinvshape;
    // purwhypermean = inpurwhypermean;
    // purwhyperinvconc = inpurwhyperinvconc;
    poswhypermean = inposwhypermean;
    poswhyperinvconc = inposwhyperinvconc;
}

// void MultiGeneCodonM2aModel::UpdateNucMatrix()	{
//    nucmatrix->CopyStationary((*nucstatarray)[0]);
//    nucmatrix->CorruptMatrix();
//}

void MultiGeneAAMutSelM2aModel::SetMixtureArrays() {
    // double puromalpha = puromhypermean / puromhyperinvconc;
    // double purombeta = (1-puromhypermean) / puromhyperinvconc;
    // puromarray->SetAlpha(puromalpha);
    // puromarray->SetBeta(purombeta);

    double dposomalpha = 1.0 / dposomhyperinvshape;
    double dposombeta = dposomalpha / dposomhypermean;
    dposomarray->SetShape(dposomalpha);
    dposomarray->SetScale(dposombeta);
    dposomarray->PriorResample(*poswarray);
    // necessary after changing some dposom values
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        // geneprocess[gene]->SetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
        geneprocess[gene]->SetMixtureParameters((*dposomarray)[gene], (*poswarray)[gene]);
    }

    // double purwalpha = purwhypermean / purwhyperinvconc;
    // double purwbeta = (1-purwhypermean) / purwhyperinvconc;
    // purwarray->SetAlpha(purwalpha);
    // purwarray->SetBeta(purwbeta);

    double poswalpha = poswhypermean / poswhyperinvconc;
    double poswbeta = (1 - poswhypermean) / poswhyperinvconc;
    poswarray->SetPi(pi);
    poswarray->SetAlpha(poswalpha);
    poswarray->SetBeta(poswbeta);
}

//-------------------
// Traces and Monitors
// ------------------

void MultiGeneAAMutSelM2aModel::TraceHeader(ostream &os) const {
    os << "#logprior\tlnL";
    if (blmode == 2) {
        os << "\tlength";
    } else {
        os << "\tmeanlength\tstdev";
    }
    os << "\tpi";
    os << "\tnposfrac";
    // os << "\tpurommean\tinvconc";
    os << "\tdposommean\tinvshape";
    // os << "\tpurwmean\tinvconc";
    os << "\tposwmean\tinvconc";
    os << "\tstatent";
    os << "\trrent";
    /*if (nucmode != 2)   {
        os << "\tstdevrr\tcenter\thyperinvconc";
        os << "\tstdevstat\tcenter\thyperinvconc";
        if (nucmode == 1)   {
            //os << "\tnucrracc1\tnucrracc2\tnucrracc3";
            //os << "\tnucstatacc1\tnucstatacc2\tnucstatacc3";
            os << "\tnucrrlogprob\tnucstatlogprob";
            os << "\tnucrrsuffstatlogprob\tnucstatsuffstatlogprob";
        }
    }*/
    os << '\n';
}

void MultiGeneAAMutSelM2aModel::Trace(ostream &os) const {
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
    // os << '\t' << puromhypermean << '\t' << puromhyperinvconc;
    os << '\t' << dposomhypermean << '\t' << dposomhyperinvshape;
    // os << '\t' << purwhypermean << '\t' << purwhyperinvconc;
    os << '\t' << poswhypermean << '\t' << poswhyperinvconc;
    // os << '\t' << nucstatarray->GetMeanEntropy();
    // os << '\t' << nucrelratearray->GetMeanEntropy();
    /*if (nucmode != 2)   {
        os << '\t' << sqrt(GetVarNucRelRate()) << '\t' <<
    Random::GetEntropy(nucrelratehypercenter) << '\t' << nucrelratehyperinvconc;
        os << '\t' << sqrt(GetVarNucStat()) << '\t' <<
    Random::GetEntropy(nucstathypercenter) << '\t' << nucstathyperinvconc; if
    (nucmode == 1)   {
            //os << '\t' << 100*((double) nucrracc1) / nucrrtot1;
            //os << '\t' << 100*((double) nucrracc2) / nucrrtot2;
            //os << '\t' << 100*((double) nucrracc3) / nucrrtot3;
            //os << '\t' << 100*((double) nucstatacc1) / nucstattot1;
            //os << '\t' << 100*((double) nucstatacc2) / nucstattot2;
            //os << '\t' << 100*((double) nucstatacc3) / nucstattot3;
            //nucrracc1 = nucrracc2 = nucrracc3 = nucrrtot1 = nucrrtot2 =
    nucrrtot3 = 0;
            //nucstatacc1 = nucstatacc2 = nucstatacc3 = nucstattot1 = nucstattot2
    = nucstattot3 = 0; os << '\t' << nucrelratearray->GetLogProb(); os << '\t' <<
    nucstatarray->GetLogProb(); os << '\t' <<
    nucrelratesuffstat.GetLogProb(nucrelratehypercenter,1.0/nucrelratehyperinvconc);
            os << '\t' <<
    nucstatsuffstat.GetLogProb(nucstathypercenter,1.0/nucstathyperinvconc);
        }
    }*/
    os << '\n';
    os.flush();
}

void MultiGeneAAMutSelM2aModel::TracePosWeight(ostream &os) const {
    for (int gene = 0; gene < Ngene; gene++) {
        os << poswarray->GetVal(gene) << '\t';
    }
    os << '\n';
    os.flush();
}

void MultiGeneAAMutSelM2aModel::TracePosOm(ostream &os) const {
    for (int gene = 0; gene < Ngene; gene++) {
        os << 1 + dposomarray->GetVal(gene) << '\t';
    }
    os << '\n';
    os.flush();
}

int MultiGeneAAMutSelM2aModel::GetNpos() const { return GetNgene() - poswarray->GetNullSet(); }

double MultiGeneAAMutSelM2aModel::GetMeanTotalLength() const {
    double tot = 0;
    for (int j = 0; j < Nbranch; j++) {
        tot += branchlength->GetVal(j);
    }
    return tot;
}

double MultiGeneAAMutSelM2aModel::GetMeanLength() const {
    if (blmode == 2) {
        cerr << "error: in getvarlength\n";
        exit(1);
    }

    return branchlengtharray->GetMeanLength();
}

double MultiGeneAAMutSelM2aModel::GetVarLength() const {
    if (blmode == 2) {
        cerr << "error: in getvarlength\n";
        exit(1);
    }

    return branchlengtharray->GetVarLength();
}

/*double MultiGeneCodonM2aModel::GetVarNucRelRate() const {

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

double MultiGeneCodonM2aModel::GetVarNucStat() const {

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
}*/

//-------------------
// Log Priors and likelihood
// ------------------

double MultiGeneAAMutSelM2aModel::GetLogPrior() const {
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
    /*if (nucmode == 2)   {
        total += GlobalNucRatesLogPrior();
    }
    else if (nucmode == 1)  {
        total += GeneNucRatesHyperLogPrior();
    }
    else    {
        // nothing: everything accounted for by gene component
    }*/

    // mixture
    total += MixtureHyperLogPrior();
    // already accounted for by gene component
    // total += MixtureLogPrior();

    return total;
}

//-------------------
// Moves
// ------------------

void MultiGeneAAMutSelM2aModel::MasterMove() {
    int nrep = 30;

    for (int rep = 0; rep < nrep; rep++) {
        // mixture hyperparameters
        MasterReceiveMixtureHyperSuffStat();
        MoveMixtureHyperParameters();
        MasterSendMixtureHyperParameters();

        // global branch lengths, or gene branch lengths hyperparameters
        if (blmode == 2) {
            MasterReceiveBranchLengthsSuffStat();
            ResampleBranchLengths();
            MoveLambda();
            MasterSendGlobalBranchLengths();
        } else if (blmode == 1) {
            MasterReceiveBranchLengthsHyperSuffStat();
            MoveBranchLengthsHyperParameters();
            MasterSendBranchLengthsHyperParameters();
        }

        // global nucrates, or gene nucrates hyperparameters
        /*if (nucmode == 2)   {
            MasterReceiveNucPathSuffStat();
            MoveNucRates();
            MasterSendGlobalNucRates();
        }
        else if (nucmode == 1)  {
            MasterReceiveNucRatesHyperSuffStat();
            MoveNucRatesHyperParameters();
            MasterSendNucRatesHyperParameters();
        }*/
    }
    burnin++;
    if (blmode != 2) {
        MasterReceiveGeneBranchLengths();
    }
    // if (nucmode != 2)   {
    //    MasterReceiveGeneNucRates();
    //}
    MasterReceiveMixture();
    MasterReceiveLogProbs();
}

// slave move
void MultiGeneAAMutSelM2aModel::SlaveMove() {
    GeneResampleSub(1.0);

    int nrep = 30;

    for (int rep = 0; rep < nrep; rep++) {
        // gene specific mixture parameters
        // possibly branch lengths and nuc rates (if mode == 1 or 2)
        MoveGeneParameters(1.0);

        // mixture hyperparameters
        SlaveSendMixtureHyperSuffStat();
        SlaveReceiveMixtureHyperParameters();
        SetMixtureArrays();

        // global branch lengths, or gene branch lengths hyperparameters
        if (blmode == 2) {
            SlaveSendBranchLengthsSuffStat();
            SlaveReceiveGlobalBranchLengths();
        } else if (blmode == 1) {
            SlaveSendBranchLengthsHyperSuffStat();
            SlaveReceiveBranchLengthsHyperParameters();
        }

        // global nucrates, or gene nucrates hyperparameters
        /*if (nucmode == 2)   {
            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalNucRates();
        }
        else if (nucmode == 1)  {
            SlaveSendNucRatesHyperSuffStat();
            SlaveReceiveNucRatesHyperParameters();
        }*/
    }
    burnin++;

    // collect current state
    if (blmode != 2) {
        SlaveSendGeneBranchLengths();
    }
    // if (nucmode != 2)   {
    //    SlaveSendGeneNucRates();
    //}
    SlaveSendMixture();
    SlaveSendLogProbs();
}

void MultiGeneAAMutSelM2aModel::GeneResampleSub(double frac) {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->ResampleSub(frac);
    }
}

void MultiGeneAAMutSelM2aModel::MoveGeneParameters(int nrep) {
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->MoveParameters(nrep);
        // geneprocess[gene]->GetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
        geneprocess[gene]->GetMixtureParameters((*dposomarray)[gene], (*poswarray)[gene]);
        if (blmode != 2) {
            geneprocess[gene]->GetBranchLengths((*branchlengtharray)[gene]);
        }
        // if (nucmode != 2)    {
        //    geneprocess[gene]->GetNucRates((*nucrelratearray)[gene],(*nucstatarray)[gene]);
        //}
    }
}

void MultiGeneAAMutSelM2aModel::ResampleBranchLengths() {
    branchlength->GibbsResample(*lengthpathsuffstatarray);
}

void MultiGeneAAMutSelM2aModel::MoveLambda() {
    hyperlengthsuffstat.Clear();
    hyperlengthsuffstat.AddSuffStat(*branchlength);
    ScalingMove(lambda, 1.0, 10, &MultiGeneAAMutSelM2aModel::LambdaHyperLogProb,
                &MultiGeneAAMutSelM2aModel::NoUpdate, this);
    ScalingMove(lambda, 0.3, 10, &MultiGeneAAMutSelM2aModel::LambdaHyperLogProb,
                &MultiGeneAAMutSelM2aModel::NoUpdate, this);
    branchlength->SetScale(lambda);
}

void MultiGeneAAMutSelM2aModel::MoveBranchLengthsHyperParameters() {
    for (int j = 0; j < Nbranch; j++) {
        BranchLengthsHyperScalingMove(1.0, 10);
        BranchLengthsHyperScalingMove(0.3, 10);
    }

    ScalingMove(blhyperinvshape, 1.0, 10, &MultiGeneAAMutSelM2aModel::BranchLengthsHyperLogProb,
                &MultiGeneAAMutSelM2aModel::NoUpdate, this);
    ScalingMove(blhyperinvshape, 0.3, 10, &MultiGeneAAMutSelM2aModel::BranchLengthsHyperLogProb,
                &MultiGeneAAMutSelM2aModel::NoUpdate, this);

    branchlengtharray->SetShape(1.0 / blhyperinvshape);
    MoveLambda();
}

double MultiGeneAAMutSelM2aModel::BranchLengthsHyperScalingMove(double tuning, int nrep) {
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

/*void MultiGeneCodonM2aModel::MoveNucRatesHyperParameters()    {

    ProfileMove(nucrelratehypercenter,1.0,1,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ProfileMove(nucrelratehypercenter,0.3,1,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ProfileMove(nucrelratehypercenter,0.1,3,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ScalingMove(nucrelratehyperinvconc,1.0,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ScalingMove(nucrelratehyperinvconc,0.3,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ScalingMove(nucrelratehyperinvconc,0.03,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);

    ProfileMove(nucstathypercenter,1.0,1,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ProfileMove(nucstathypercenter,0.3,1,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ProfileMove(nucstathypercenter,0.1,2,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ScalingMove(nucstathyperinvconc,1.0,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ScalingMove(nucstathyperinvconc,0.3,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    ScalingMove(nucstathyperinvconc,0.03,10,&MultiGeneCodonM2aModel::NucRatesHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);

    nucrelratearray->SetConcentration(1.0/nucrelratehyperinvconc);
    nucstatarray->SetConcentration(1.0 / nucstathyperinvconc);
}*/

void MultiGeneAAMutSelM2aModel::MoveMixtureHyperParameters() {
    /*if (purommode == 1) {
        SlidingMove(puromhypermean,1.0,10,0,1,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
        SlidingMove(puromhypermean,0.3,10,0,1,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
        ScalingMove(puromhyperinvconc,1.0,10,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
        ScalingMove(puromhyperinvconc,0.3,10,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    }*/

    if (dposommode == 1) {
        SlidingMove(dposomhypermean, 3.0, 10, 0, 1, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
        SlidingMove(dposomhypermean, 1.0, 10, 0, 1, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
        SlidingMove(dposomhypermean, 0.3, 10, 0, 1, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
        ScalingMove(dposomhyperinvshape, 1.0, 10, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
        ScalingMove(dposomhyperinvshape, 0.3, 10, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
    }

    if (poswmode == 1) {
        SlidingMove(poswhypermean, 1.0, 10, 0, 1, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
        SlidingMove(poswhypermean, 0.3, 10, 0, 1, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
        ScalingMove(poswhyperinvconc, 1.0, 10, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
        ScalingMove(poswhyperinvconc, 0.3, 10, &MultiGeneAAMutSelM2aModel::MixtureHyperLogProb,
                    &MultiGeneAAMutSelM2aModel::NoUpdate, this);
    }

    /*if (purwmode == 1)  {
        SlidingMove(purwhypermean,1.0,10,0,1,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
        SlidingMove(purwhypermean,0.3,10,0,1,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
        ScalingMove(purwhyperinvconc,1.0,10,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
        ScalingMove(purwhyperinvconc,0.3,10,&MultiGeneCodonM2aModel::MixtureHyperLogProb,&MultiGeneCodonM2aModel::NoUpdate,this);
    }*/

    if (burnin > 10) {
        if (pihyperinvconc) {
            ResamplePi();
        }
    }
}

void MultiGeneAAMutSelM2aModel::ResamplePi() {
    int n0 = poswsuffstat.GetN0();
    int n1 = poswsuffstat.GetN1();
    if ((n0 + n1) != Ngene) {
        cerr << "error in resample pi\n";
        exit(1);
    }
    double pialpha = pihypermean / pihyperinvconc;
    double pibeta = (1 - pihypermean) / pihyperinvconc;
    double a0 = Random::sGamma(pialpha + n0);
    double a1 = Random::sGamma(pibeta + n1);
    pi = a1 / (a0 + a1);
}

/*void MultiGeneCodonM2aModel::MoveNucRates()    {

    vector<double>& nucrelrate = (*nucrelratearray)[0];
    ProfileMove(nucrelrate,0.1,1,3,&MultiGeneCodonM2aModel::NucRatesLogProb,&MultiGeneCodonM2aModel::UpdateNucMatrix,this);
    ProfileMove(nucrelrate,0.03,3,3,&MultiGeneCodonM2aModel::NucRatesLogProb,&MultiGeneCodonM2aModel::UpdateNucMatrix,this);
    ProfileMove(nucrelrate,0.01,3,3,&MultiGeneCodonM2aModel::NucRatesLogProb,&MultiGeneCodonM2aModel::UpdateNucMatrix,this);

    vector<double>& nucstat = (*nucstatarray)[0];
    ProfileMove(nucstat,0.1,1,3,&MultiGeneCodonM2aModel::NucRatesLogProb,&MultiGeneCodonM2aModel::UpdateNucMatrix,this);
    ProfileMove(nucstat,0.01,1,3,&MultiGeneCodonM2aModel::NucRatesLogProb,&MultiGeneCodonM2aModel::UpdateNucMatrix,this);
}*/

//-------------------
// MPI Send Receive
// ------------------

void MultiGeneAAMutSelM2aModel::MasterSendGlobalBranchLengths() { MasterSendGlobal(*branchlength); }

void MultiGeneAAMutSelM2aModel::SlaveReceiveGlobalBranchLengths() {
    SlaveReceiveGlobal(*branchlength);
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetBranchLengths(*branchlength);
    }
}

void MultiGeneAAMutSelM2aModel::MasterSendBranchLengthsHyperParameters() {
    MasterSendGlobal(*branchlength, blhyperinvshape);
}

void MultiGeneAAMutSelM2aModel::SlaveReceiveBranchLengthsHyperParameters() {
    SlaveReceiveGlobal(*branchlength, blhyperinvshape);
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetBranchLengthsHyperParameters(*branchlength, blhyperinvshape);
    }
}

void MultiGeneAAMutSelM2aModel::MasterSendGeneBranchLengths() {
    MasterSendGeneArray(*branchlengtharray);
}

void MultiGeneAAMutSelM2aModel::SlaveReceiveGeneBranchLengths() {
    SlaveReceiveGeneArray(*branchlengtharray);
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->SetBranchLengths(branchlengtharray->GetVal(gene));
    }
}

void MultiGeneAAMutSelM2aModel::SlaveSendGeneBranchLengths() {
    SlaveSendGeneArray(*branchlengtharray);
}

void MultiGeneAAMutSelM2aModel::MasterReceiveGeneBranchLengths() {
    MasterReceiveGeneArray(*branchlengtharray);
}

/*void MultiGeneCodonM2aModel::MasterSendGlobalNucRates()   {

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
*/

void MultiGeneAAMutSelM2aModel::MasterSendMixtureHyperParameters() {
    MasterSendGlobal(mixhyperparam);
}

void MultiGeneAAMutSelM2aModel::SlaveReceiveMixtureHyperParameters() {
    SlaveReceiveGlobal(mixhyperparam);

    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        // geneprocess[gene]->SetMixtureHyperParameters(puromhypermean,puromhyperinvconc,dposomhypermean,dposomhyperinvshape,pi,purwhypermean,purwhyperinvconc,poswhypermean,poswhyperinvconc);
        geneprocess[gene]->SetMixtureHyperParameters(dposomhypermean, dposomhyperinvshape, pi,
                                                     poswhypermean, poswhyperinvconc);
    }
}

/*void MultiGeneCodonM2aModel::MasterSendNucRatesHyperParameters()   {

    MasterSendGlobal(nucrelratehypercenter,nucrelratehyperinvconc);
    MasterSendGlobal(nucstathypercenter,nucstathyperinvconc);
}

void MultiGeneCodonM2aModel::SlaveReceiveNucRatesHyperParameters()   {

    SlaveReceiveGlobal(nucrelratehypercenter,nucrelratehyperinvconc);
    SlaveReceiveGlobal(nucstathypercenter,nucstathyperinvconc);

    for (int gene=0; gene<GetLocalNgene(); gene++)   {
        geneprocess[gene]->SetNucRatesHyperParameters(nucrelratehypercenter,nucrelratehyperinvconc,nucstathypercenter,nucstathyperinvconc);
    }
}*/

void MultiGeneAAMutSelM2aModel::SlaveSendBranchLengthsSuffStat() {
    lengthpathsuffstatarray->Clear();
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        geneprocess[gene]->CollectLengthSuffStat();
        lengthpathsuffstatarray->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
    }
    SlaveSendAdditive(*lengthpathsuffstatarray);
}

void MultiGeneAAMutSelM2aModel::MasterReceiveBranchLengthsSuffStat() {
    lengthpathsuffstatarray->Clear();
    MasterReceiveAdditive(*lengthpathsuffstatarray);
}

void MultiGeneAAMutSelM2aModel::SlaveSendBranchLengthsHyperSuffStat() {
    lengthhypersuffstatarray->Clear();
    lengthhypersuffstatarray->AddSuffStat(*branchlengtharray);
    SlaveSendAdditive(*lengthhypersuffstatarray);
}

void MultiGeneAAMutSelM2aModel::MasterReceiveBranchLengthsHyperSuffStat() {
    lengthhypersuffstatarray->Clear();
    MasterReceiveAdditive(*lengthhypersuffstatarray);
}

void MultiGeneAAMutSelM2aModel::SlaveSendMixtureHyperSuffStat() {
    // puromsuffstat.Clear();
    // puromarray->AddSuffStat(puromsuffstat);
    // SlaveSendAdditive(puromsuffstat);

    dposomsuffstat.Clear();
    dposomsuffstat.AddSuffStat(*dposomarray, *poswarray);
    SlaveSendAdditive(dposomsuffstat);

    // purwsuffstat.Clear();
    // purwarray->AddSuffStat(purwsuffstat);
    // SlaveSendAdditive(purwsuffstat);

    poswsuffstat.Clear();
    poswarray->AddSuffStat(poswsuffstat);
    SlaveSendAdditive(poswsuffstat);
}

void MultiGeneAAMutSelM2aModel::MasterReceiveMixtureHyperSuffStat() {
    // puromsuffstat.Clear();
    // MasterReceiveAdditive(puromsuffstat);
    dposomsuffstat.Clear();
    MasterReceiveAdditive(dposomsuffstat);
    // purwsuffstat.Clear();
    // MasterReceiveAdditive(purwsuffstat);
    poswsuffstat.Clear();
    MasterReceiveAdditive(poswsuffstat);
}

/*void MultiGeneCodonM2aModel::SlaveSendNucPathSuffStat()  {

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
}*/

void MultiGeneAAMutSelM2aModel::SlaveSendMixture() {
    // SlaveSendGeneArray(*puromarray,*dposomarray);
    // SlaveSendGeneArray(*purwarray,*poswarray);
    SlaveSendGeneArray(*dposomarray);
    SlaveSendGeneArray(*poswarray);
}

void MultiGeneAAMutSelM2aModel::MasterReceiveMixture() {
    // MasterReceiveGeneArray(*puromarray,*dposomarray);
    // MasterReceiveGeneArray(*purwarray,*poswarray);
    MasterReceiveGeneArray(*dposomarray);
    MasterReceiveGeneArray(*poswarray);
}

void MultiGeneAAMutSelM2aModel::MasterSendMixture() {
    // MasterSendGeneArray(*puromarray,*dposomarray);
    // MasterSendGeneArray(*purwarray,*poswarray);
    MasterSendGeneArray(*dposomarray);
    MasterSendGeneArray(*poswarray);
}

void MultiGeneAAMutSelM2aModel::SlaveReceiveMixture() {
    SlaveReceiveGeneArray(*dposomarray);
    SlaveReceiveGeneArray(*poswarray);

    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        // geneprocess[gene]->SetMixtureParameters((*puromarray)[gene],(*dposomarray)[gene],(*purwarray)[gene],(*poswarray)[gene]);
        geneprocess[gene]->SetMixtureParameters((*dposomarray)[gene], (*poswarray)[gene]);
    }
}

void MultiGeneAAMutSelM2aModel::SlaveSendLogProbs() {
    GeneLogPrior = 0;
    lnL = 0;
    for (int gene = 0; gene < GetLocalNgene(); gene++) {
        GeneLogPrior += geneprocess[gene]->GetLogPrior();
        lnL += geneprocess[gene]->GetLogLikelihood();
    }
    SlaveSendAdditive(GeneLogPrior);
    SlaveSendAdditive(lnL);
}

void MultiGeneAAMutSelM2aModel::MasterReceiveLogProbs() {
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
        cerr << "error in MultiGeneCodonM2aModel::MasterTraceSitesPostProb: non
matching number of sites\n"; exit(1);
    }
}

void MultiGeneCodonM2aModel::SlaveTraceSitesPostProb()  {

}
*/

void MultiGeneAAMutSelM2aModel::MasterTraceSitesPostProb(ostream &os) {
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
    }
    os << '\n';
    os.flush();
}

void MultiGeneAAMutSelM2aModel::SlaveTraceSitesPostProb() {
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
