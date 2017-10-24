#include "AAMutSelM2aModel.hpp"

// constuction allocation
//

AAMutSelM2aModel::AAMutSelM2aModel(string datafile, string treefile, double inpi) : aahypersuffstat(20)	{

    //std::cerr << "in principal constructor" << std::endl;

    blmode = 0;
    nucmode = 0;
    data = new FileSequenceAlignment(datafile);
    codondata = new CodonSequenceAlignment(data, true);
    pi = inpi;

    Nsite = codondata->GetNsite();    // # columns
    Ntaxa = codondata->GetNtaxa();

    //std::cerr << "-- Number of sites: " << Nsite << std::endl;

    taxonset = codondata->GetTaxonSet();

    // get tree from file (newick format)
    tree = new Tree(treefile);

    // check whether tree and data fits together
    tree->RegisterWith(taxonset);

    tree->SetIndices();
    Nbranch = tree->GetNbranch();
    //std::cerr << "end of  principal constructor" << std::endl;
}


void AAMutSelM2aModel::Unfold()   {

    phyloprocess->Unfold();
    phyloprocess->ResampleSub();
}

void AAMutSelM2aModel::Allocate()	{

    lambda = 10.0;
    blhypermean = new BranchIIDGamma(*tree,1.0,lambda);
    blhypermean->SetAllBranches(1.0 / lambda);
    blhyperinvshape = 1.0;
    branchlength = new GammaWhiteNoise(*tree,*blhypermean,1.0/blhyperinvshape);

    //purom = 0.5;
    //puromhypermean = 0.5;
    //puromhyperinvconc = 0.5;

    dposom = 1.0;
    dposomhypermean = 0.5;
    dposomhyperinvshape = 0.5;

    //purw = 0.1;
    //purwhypermean = 0.5;
    //purwhyperinvconc = 0.5;

    if (! pi)   {
        posw = 0;
        poswhypermean = 0;
        poswhyperinvconc = 0;
    }
    else    {
        posw = 0.1;
        poswhypermean = 0.5;
        poswhyperinvconc = 0.5;
    }

    //componentomegaarray = new M2aMix(purom,dposom+1,purw,posw);
    componentomegaarray = new MechM2aMix(dposom+1,posw);
    sitealloc = new MultinomialAllocationVector(GetNsite(),componentomegaarray->GetWeights());
    siteomegaarray = new ConstMixtureArray<double> (componentomegaarray,sitealloc);
    sitepostprobarray.assign(GetNsite(),vector<double>(2,0));

    nucrelratehypercenter.assign(Nrr,1.0/Nrr);
    nucrelratehyperinvconc = 1.0 / Nrr;

    nucstathypercenter.assign(Nnuc,1.0/Nnuc);
    nucstathyperinvconc = 1.0 / Nnuc;

    nucrelrate.assign(Nrr,0);
    Random::DirichletSample(nucrelrate,nucrelratehypercenter,1.0/nucrelratehyperinvconc);

    nucstat.assign(Nnuc,0);
    Random::DirichletSample(nucstat,nucstathypercenter,1.0/nucstathyperinvconc);

    nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

    aacenter.assign(Naa,1.0/Naa);
    aainvconc = 1.0/Naa;
    aafitnessarray = new IIDDirichlet(Nsite,aacenter,1.0/aainvconc);

    sitecodonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,aafitnessarray,siteomegaarray);

    //sitesubmatrixarray = new ConstMixtureArray<SubMatrix>(componentcodonmatrixarray,sitealloc);
    //sitecodonmatrixarray = new ConstMixtureArray<MGOmegaCodonSubMatrix>(componentcodonmatrixarray,sitealloc);

    phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitecodonmatrixarray);

    lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);
    sitepathsuffstatarray = new PathSuffStatArray(GetNsite());
    //componentpathsuffstatarray = new PathSuffStatArray(2);
    siteomegasuffstatarray = new OmegaSuffStatArray(GetNsite());
}

void AAMutSelM2aModel::Update()    {

    blhypermean->SetAllBranches(1.0/lambda);
    cerr << "check that\n";
    componentomegaarray->SetParameters(dposom+1,posw);
    UpdateMatrices();
    /*
    UpdateMatrices();
    componentomegaarray->SetParameters(purom,dposom+1,purw,posw);
    */
    ResampleAlloc();
}

// setting model features and (hyper)parameters
//

void AAMutSelM2aModel::SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
    branchlength->Copy(inbranchlength);
}

void AAMutSelM2aModel::GetBranchLengths(BranchArray<double>& inbranchlength) const   {
    inbranchlength.Copy(*branchlength);
}

void AAMutSelM2aModel::SetBranchLengthsHyperParameters(const ConstBranchArray<double>& inblmean, double inblinvshape) {
    blhypermean->Copy(inblmean);
    branchlength->SetShape(1.0 / blhyperinvshape);
}

void AAMutSelM2aModel::SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
    nucrelrate = innucrelrate;
    nucstat = innucstat;
    UpdateMatrices();
}

void AAMutSelM2aModel::GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const {
    innucrelrate = nucrelrate;
    innucstat = nucstat;
}

void AAMutSelM2aModel::SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc, const std::vector<double>& innucstathypercenter, double innucstathyperinvconc) {
    nucrelratehypercenter = innucrelratehypercenter;
    nucrelratehyperinvconc = innucrelratehyperinvconc;
    nucstathypercenter = innucstathypercenter;
    nucstathyperinvconc = innucstathyperinvconc;
}

void AAMutSelM2aModel::SetMixtureParameters(double indposom, double inposw)    {
    dposom = indposom;
    posw = inposw;
    componentomegaarray->SetParameters(dposom+1,posw);
}

void AAMutSelM2aModel::GetMixtureParameters(double& indposom, double& inposw)  const {
    indposom = dposom;
    inposw = posw;
}

void AAMutSelM2aModel::SetMixtureHyperParameters(double indposomhypermean, double indposomhyperinvshape, double inpi, double inposwhypermean, double inposwhyperinvconc)  {
    dposomhypermean = indposomhypermean;
    dposomhyperinvshape = indposomhyperinvshape;
    pi = inpi;
    poswhypermean = inposwhypermean;
    poswhyperinvconc = inposwhyperinvconc;
}

// 
// Matrices
//

void AAMutSelM2aModel::UpdateNucMatrix()	{
    nucmatrix->CopyStationary(nucstat);
    nucmatrix->CorruptMatrix();
}

void AAMutSelM2aModel::UpdateCodonMatrices()	{
    sitecodonmatrixarray->UpdateCodonMatrices();
}

void AAMutSelM2aModel::UpdateCodonMatrix(int site)    {
    (*sitecodonmatrixarray)[site].CorruptMatrix();
}
    
void AAMutSelM2aModel::UpdateMatrices()   {
    UpdateNucMatrix();
    UpdateCodonMatrices();
}

//
// Likelihood
//

double AAMutSelM2aModel::GetLogLikelihood() const {
    // return GetIntegratedLogLikelihood();
    return phyloprocess->GetLogProb();
}

//
// Suff Stat and suffstatlogprobs
//

const PoissonSuffStatBranchArray* AAMutSelM2aModel::GetLengthSuffStatArray() const {
    return lengthsuffstatarray;
}

double AAMutSelM2aModel::LambdaHyperSuffStatLogProb() const {
    return lambdasuffstat.GetLogProb(1.0,lambda);
}

//const NucPathSuffStat& AAMutSelM2aModel::GetNucPathSuffStat() const {
//    return nucpathsuffstat;
//}

//double AAMutSelM2aModel::NucRatesSuffStatLogProb() const {
//   return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
//}

double AAMutSelM2aModel::PathSuffStatLogProb() const {
    //return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
    return sitepathsuffstatarray->GetLogProb(*sitecodonmatrixarray);
}

double AAMutSelM2aModel::PathSuffStatLogProb(int site) const {
    return sitepathsuffstatarray->GetVal(site).GetLogProb(sitecodonmatrixarray->GetVal(site));
}

double AAMutSelM2aModel::OmegaSuffStatLogProb() const {
    componentomegaarray->SetParameters(dposom+1,posw);
    return componentomegaarray->GetPostProbArray(*siteomegasuffstatarray,sitepostprobarray);
}
    
double AAMutSelM2aModel::AAHyperSuffStatLogProb() const   {
    return aahypersuffstat.GetLogProb(aacenter,1.0/aainvconc);
}

//
// Priors
//

double AAMutSelM2aModel::GetLogPrior() const {

    double total = 0;

    if (! FixedBranchLengths()) {
        total += BranchLengthsLogPrior();
    }
    if (! FixedNucRates())  {
        total += NucRatesLogPrior();
    }
    total += OmegaLogPrior();
    total += AAHyperLogPrior();
    total += AALogPrior();
    return total;
}

double AAMutSelM2aModel::BranchLengthsLogPrior() const  {

    double total = 0;
    if (blmode == 0)    {
        total += LambdaHyperLogPrior();
    }
    total += branchlength->GetLogProb();
    return total;
}

double AAMutSelM2aModel::LambdaHyperLogPrior()	 const {
    return -lambda / 10;
}

double AAMutSelM2aModel::NucRatesLogPrior() const {

    double total = 0;
    total += Random::logDirichletDensity(nucrelrate,nucrelratehypercenter,1.0/nucrelratehyperinvconc);
    total += Random::logDirichletDensity(nucstat,nucstathypercenter,1.0/nucstathyperinvconc);
    return total;
}

//
// Hyper priors for omega mixture
//

double AAMutSelM2aModel::OmegaLogPrior() const {
    double total = 0;
    //total += PurOmegaLogProb();
    total += PosOmegaLogProb();
    //total += PurWeightLogProb();
    total += PosWeightLogProb();
    return total;
}

// Beta prior for purifmean
//double CodonM2aModel::PurOmegaLogProb()  const {
//    double alpha = puromhypermean / puromhyperinvconc;
//    double beta = (1-puromhypermean) / puromhyperinvconc;
//    return Random::logBetaDensity(purom,alpha,beta);
//}

// Gamma prior for dposom
double AAMutSelM2aModel::PosOmegaLogProb() const {
    double alpha = 1.0 / dposomhyperinvshape;
    double beta = alpha / dposomhypermean;
    return Random::logGammaDensity(dposom,alpha,beta);
}

// Beta prior for purw
//double CodonM2aModel::PurWeightLogProb() const {
//    double alpha = purwhypermean / purwhyperinvconc;
//    double beta = (1 - purwhypermean) / purwhyperinvconc;
//    return Random::logBetaDensity(purw,alpha,beta);
//}

// mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1 - pi) for posw
double AAMutSelM2aModel::PosWeightLogProb() const {
    if (posw)   {
        if (! pi)   {
            cerr << "in PosWeightLogProb: pi == 0 and posw > 0\n";
            exit(1);
        }

        double alpha = poswhypermean / poswhyperinvconc;
        double beta = (1 - poswhypermean) / poswhyperinvconc;
        return log(pi) + Random::logBetaDensity(posw,alpha,beta);
    }
    else    {
        return log(1-pi);
    }
}

// Bernoulli for whether posw == 0 or > 0
double AAMutSelM2aModel::PosSwitchLogProb() const {
    if (posw)   {
        return log(pi);
    }
    return log(1-pi);
}

// AA priors

    double AAMutSelM2aModel::AAHyperLogPrior() const {
        return -aainvconc;
    }

    double AAMutSelM2aModel::AALogPrior() const {
        return aafitnessarray->GetLogProb();
    }

    double AAMutSelM2aModel::AALogPrior(int i) const {
        return aafitnessarray->GetLogProb(i);
    }


//
//  Moves 
//

double AAMutSelM2aModel::Move()	{

    ResampleSub(1.0);
    MoveParameters(30);
    return 1;
}

void AAMutSelM2aModel::MoveParameters(int nrep)    {

    for (int rep=0; rep<nrep; rep++)	{

        if (! FixedBranchLengths()) {
            MoveBranchLengths();
        }

        CollectPathSuffStat();

        MoveAA();
        MoveAAHyperParameters();
        MoveOmega();

        if (! FixedNucRates())  {
            UpdateMatrices();
            MoveNucRates();
        }
    }
}

void AAMutSelM2aModel::ResampleSub(double frac)  {
    UpdateMatrices();
    phyloprocess->Move(frac);
}

//
// Branch Lengths and hyperparam lambda
//

void AAMutSelM2aModel::MoveBranchLengths()    {
        ResampleBranchLengths();
        if (blmode == 0)    {
            MoveLambda();
        }
}

void AAMutSelM2aModel::ResampleBranchLengths()	{

    CollectLengthSuffStat();
    branchlength->GibbsResample(*lengthsuffstatarray);
}

void AAMutSelM2aModel::CollectLengthSuffStat()    {

    lengthsuffstatarray->Clear();
    phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
}

void AAMutSelM2aModel::MoveLambda()	{

    lambdasuffstat.Clear();
    branchlength->AddSuffStat(lambdasuffstat);

    ScalingMove(lambda,1.0,10,&AAMutSelM2aModel::LambdaHyperLogProb,&AAMutSelM2aModel::NoUpdate,this);
    ScalingMove(lambda,0.3,10,&AAMutSelM2aModel::LambdaHyperLogProb,&AAMutSelM2aModel::NoUpdate,this);

    blhypermean->SetAllBranches(1.0/lambda);
}

//
// Omega mixture 
//

void AAMutSelM2aModel::CollectPathSuffStat()	{

    sitepathsuffstatarray->Clear();
    phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
}

//void AAMutSelM2aModel::CollectComponentPathSuffStat()	{
//
//    componentpathsuffstatarray->Clear();
//    sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
//}

void AAMutSelM2aModel::MoveOmega() 	{

    CollectOmegaSuffStat();

    //SlidingMove(purom,0.1,10,0,1,&AAMutSelM2aModel::OmegaLogProb,&AAMutSelM2aModel::NoUpdate,this);
    //SlidingMove(purw,1.0,10,0,1,&AAMutSelM2aModel::OmegaLogProb,&AAMutSelM2aModel::NoUpdate,this);
    if (pi != 0)    {
        ScalingMove(dposom,1.0,10,&AAMutSelM2aModel::OmegaLogProb,&AAMutSelM2aModel::NoUpdate,this);
        SlidingMove(posw,1.0,10,0,1,&AAMutSelM2aModel::OmegaLogProb,&AAMutSelM2aModel::NoUpdate,this);
    }
    if ((pi != 0) && (pi != 1))    {
        SwitchPosWeight(10);
    }
    ResampleAlloc();
}

void AAMutSelM2aModel::CollectOmegaSuffStat()	{

    siteomegasuffstatarray->Clear();
    siteomegasuffstatarray->AddSuffStat(*sitecodonmatrixarray,*sitepathsuffstatarray);
}

void AAMutSelM2aModel::ResampleAlloc()	{
    OmegaSuffStatLogProb();
    sitealloc->GibbsResample(sitepostprobarray);
    UpdateMatrices();
}

double AAMutSelM2aModel::DrawBetaPosWeight()    {
    double alpha = poswhypermean / poswhyperinvconc;
    double beta = (1-poswhypermean) / poswhyperinvconc;
    return Random::BetaSample(alpha,beta);
}

double AAMutSelM2aModel::SwitchPosWeight(int nrep)	{

    double nacc = 0;
    double ntot = 0;
    for (int rep=0; rep<nrep; rep++)	{
        double bkposw = posw;
        double deltalogprob = - PosSwitchLogProb() - OmegaSuffStatLogProb();
        if (posw)   {
            posw = 0;
        }
        else    {
            posw = DrawBetaPosWeight();
        }
        deltalogprob += PosSwitchLogProb() + OmegaSuffStatLogProb();
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (accepted)	{
            nacc ++;
        }
        else	{
            posw = bkposw;
        }
        ntot++;
    }
    return nacc/ntot;
}

    void AAMutSelM2aModel::MoveAAHyperParameters()    {
        aahypersuffstat.Clear();
        aafitnessarray->AddSuffStat(aahypersuffstat);
        for (int rep=0; rep<10; rep++)  {
            ProfileMove(aacenter,0.1,1,10,&AAMutSelM2aModel::AAHyperLogProb,&AAMutSelM2aModel::NoUpdate,this);
            ProfileMove(aacenter,0.03,3,10,&AAMutSelM2aModel::AAHyperLogProb,&AAMutSelM2aModel::NoUpdate,this);
            ProfileMove(aacenter,0.01,3,10,&AAMutSelM2aModel::AAHyperLogProb,&AAMutSelM2aModel::NoUpdate,this);
            ScalingMove(aainvconc,1.0,10,&AAMutSelM2aModel::AAHyperLogProb,&AAMutSelM2aModel::NoUpdate,this);
            ScalingMove(aainvconc,0.3,10,&AAMutSelM2aModel::AAHyperLogProb,&AAMutSelM2aModel::NoUpdate,this);
        }
        aafitnessarray->SetCenter(aacenter);
        aafitnessarray->SetConcentration(1.0/aainvconc);
    }

    double AAMutSelM2aModel::MoveAA() {
        MoveAA(1.0,1,3);
        MoveAA(0.3,1,3);
        MoveAA(0.1,3,3);
        MoveAA(0.1,5,3);
        return 1.0;
    }

	double AAMutSelM2aModel::MoveAA(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Naa];
        for (int i=0; i<Nsite; i++) {
            vector<double>& aa = (*aafitnessarray)[i];
            for (int rep=0; rep<nrep; rep++)	{
                for (int l=0; l<Naa; l++)	{
                    bk[l] = aa[l];
                }
                double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                double loghastings = Random::ProfileProposeMove(aa,Naa,tuning,n);
                deltalogprob += loghastings;
                UpdateCodonMatrix(i);
                deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted)	{
                    nacc ++;
                }
                else	{
                    for (int l=0; l<Naa; l++)	{
                        aa[l] = bk[l];
                    }
                    UpdateCodonMatrix(i);
                }
                ntot++;
            }
        }
        return 1.0;
    }
//
// nucleotide parameters
//

/*void AAMutSelM2aModel::CollectNucPathSuffStat()   {
    UpdateMatrices();
    nucpathsuffstat.Clear();
    nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
}*/

void AAMutSelM2aModel::MoveNucRates()	{

    //CollectComponentPathSuffStat();
    //CollectNucPathSuffStat();

    ProfileMove(nucrelrate,0.1,1,3,&AAMutSelM2aModel::NucRatesLogProb,&AAMutSelM2aModel::UpdateMatrices,this);
    ProfileMove(nucrelrate,0.03,3,3,&AAMutSelM2aModel::NucRatesLogProb,&AAMutSelM2aModel::UpdateMatrices,this);
    ProfileMove(nucrelrate,0.01,3,3,&AAMutSelM2aModel::NucRatesLogProb,&AAMutSelM2aModel::UpdateMatrices,this);

    ProfileMove(nucstat,0.1,1,3,&AAMutSelM2aModel::NucRatesLogProb,&AAMutSelM2aModel::UpdateMatrices,this);
    ProfileMove(nucstat,0.01,1,3,&AAMutSelM2aModel::NucRatesLogProb,&AAMutSelM2aModel::UpdateMatrices,this);

    //UpdateMatrices();
}

// summary statistics

double AAMutSelM2aModel::GetMeanOmega() const {
    //return posw*(1 + dposom) + (1-posw)*(purw*purom + (1-purw));
    return posw*(1 + dposom) + (1-posw);
}

void AAMutSelM2aModel::TraceHeader(std::ostream& os) const {
    os << "#logprior\tlnL\tlength\t";
    os << "posom\tposw\t";
    os << "statent\t";
    os << "rrent\n";
}

void AAMutSelM2aModel::Trace(ostream& os) const {	
    os << GetLogPrior() << '\t';
    os << GetLogLikelihood() << '\t';
    os << branchlength->GetTotalLength() << '\t';
    os << dposom+1 << '\t' << posw << '\t';
    os << Random::GetEntropy(nucstat) << '\t';
    os << Random::GetEntropy(nucrelrate) << '\n';
    SubMatrix::diagerr = 0;
}

void AAMutSelM2aModel::TracePostProb(ostream& os) const {
    for (int i=0; i<GetNsite(); i++)    {
        os << sitepostprobarray[i][2] << '\t';
    }
    os << '\n';
}

void AAMutSelM2aModel::GetSitesPostProb(double* array) const {
    for (int i=0; i<GetNsite(); i++)    {
        array[i] = sitepostprobarray[i][2];
        if (sitepostprobarray[i][2] < 0)    {
            cerr << "error in AAMutSelM2aModel::GetSitesPostProb: negative prob\n";
            exit(1);
        }
    }
}

void AAMutSelM2aModel::ToStream(ostream& os) const {

    os << lambda << '\n';
    os << *branchlength << '\n';
    //os << purom << '\t' << dposom << '\t' << purw << '\t' << posw << '\n';
    os << dposom << '\t' << posw << '\n';
    os << nucrelrate << '\n';
    os << nucstat << '\n';
}

void AAMutSelM2aModel::FromStream(istream& is) {

    is >> lambda;
    is >> *branchlength;
    //is >> purom >> dposom >> purw >> posw;
    is >> dposom >> posw;
    is >> nucrelrate;
    is >> nucstat;
}
