#include "CodonM2aModel.hpp"

// constuction allocation
//

CodonM2aModel::CodonM2aModel(string datafile, string treefile, double inpi)	{

    blmode = 0;
    nucmode = 0;
    data = new FileSequenceAlignment(datafile);
    codondata = new CodonSequenceAlignment(data, true);
    pi = inpi;

    Nsite = codondata->GetNsite();    // # columns
    Ntaxa = codondata->GetNtaxa();

    // std::cerr << "-- Number of sites: " << Nsite << std::endl;

    taxonset = codondata->GetTaxonSet();

    // get tree from file (newick format)
    tree = new Tree(treefile);

    // check whether tree and data fits together
    tree->RegisterWith(taxonset);

    tree->SetIndices();
    Nbranch = tree->GetNbranch();
}


void CodonM2aModel::Unfold()   {

    phyloprocess->Unfold();
    phyloprocess->ResampleSub();
}

void CodonM2aModel::Allocate()	{

    lambda = 10.0;
    blhypermean = new BranchIIDGamma(*tree,1.0,lambda);
    blhypermean->SetAllBranches(1.0 / lambda);
    blhyperinvshape = 1.0;
    branchlength = new GammaWhiteNoise(*tree,*blhypermean,1.0/blhyperinvshape);

    purom = 0.5;
    puromhypermean = 0.5;
    puromhyperinvconc = 0.5;

    dposom = 1.0;
    dposomhypermean = 0.5;
    dposomhyperinvshape = 0.5;

    purw = 0.1;
    purwhypermean = 0.5;
    purwhyperinvconc = 0.5;

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

    componentomegaarray = new M2aMix(purom,dposom+1,purw,posw);
    sitealloc = new MultinomialAllocationVector(GetNsite(),componentomegaarray->GetWeights());
    sitepostprobarray.assign(GetNsite(),vector<double>(3,0));

    nucrelratehypercenter.assign(Nrr,1.0/Nrr);
    nucrelratehyperinvconc = 1.0 / Nrr;

    nucstathypercenter.assign(Nnuc,1.0/Nnuc);
    nucstathyperinvconc = 1.0 / Nnuc;

    nucrelrate.assign(Nrr,0);
    Random::DirichletSample(nucrelrate,nucrelratehypercenter,1.0/nucrelratehyperinvconc);

    nucstat.assign(Nnuc,0);
    Random::DirichletSample(nucstat,nucstathypercenter,1.0/nucstathyperinvconc);

    nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

    componentcodonmatrixarray = new MGOmegaCodonSubMatrixArray((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,componentomegaarray);

    sitesubmatrixarray = new ConstMixtureArray<SubMatrix>(componentcodonmatrixarray,sitealloc);
    sitecodonmatrixarray = new ConstMixtureArray<MGOmegaCodonSubMatrix>(componentcodonmatrixarray,sitealloc);

    phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,sitesubmatrixarray);

    lengthsuffstatarray = new PoissonSuffStatBranchArray(*tree);
    sitepathsuffstatarray = new PathSuffStatArray(GetNsite());
    componentpathsuffstatarray = new PathSuffStatArray(3);
    siteomegasuffstatarray = new OmegaSuffStatArray(GetNsite());
}

// setting model features and (hyper)parameters
//

void CodonM2aModel::SetBranchLengths(const ConstBranchArray<double>& inbranchlength)    {
    branchlength->Copy(inbranchlength);
}

void CodonM2aModel::GetBranchLengths(BranchArray<double>& inbranchlength) const   {
    inbranchlength.Copy(*branchlength);
}

void CodonM2aModel::SetBranchLengthsHyperParameters(const ConstBranchArray<double>& inblmean, double inblinvshape) {
    blhypermean->Copy(inblmean);
    branchlength->SetShape(1.0 / blhyperinvshape);
}

void CodonM2aModel::SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
    nucrelrate = innucrelrate;
    nucstat = innucstat;
    UpdateMatrices();
}

void CodonM2aModel::GetNucRates(std::vector<double>& innucrelrate, std::vector<double>& innucstat) const {
    innucrelrate = nucrelrate;
    innucstat = nucstat;
}

void CodonM2aModel::SetNucRatesHyperParameters(const std::vector<double>& innucrelratehypercenter, double innucrelratehyperinvconc, const std::vector<double>& innucstathypercenter, double innucstathyperinvconc) {
    nucrelratehypercenter = innucrelratehypercenter;
    nucrelratehyperinvconc = innucrelratehyperinvconc;
    nucstathypercenter = innucstathypercenter;
    nucstathyperinvconc = innucstathyperinvconc;
}

void CodonM2aModel::SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw)    {
    purom = inpurom;
    dposom = indposom;
    purw = inpurw;
    posw = inposw;
    componentomegaarray->SetParameters(purom,dposom+1,purw,posw);
}

void CodonM2aModel::GetMixtureParameters(double& inpurom, double& indposom, double& inpurw, double& inposw)  const {
    inpurom = purom;
    indposom = dposom;
    inpurw = purw;
    inposw = posw;
}

void CodonM2aModel::SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape, double inpi, double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc)  {
    puromhypermean = inpuromhypermean;
    puromhyperinvconc = inpuromhyperinvconc;
    dposomhypermean = indposomhypermean;
    dposomhyperinvshape = indposomhyperinvshape;
    pi = inpi;
    purwhypermean = inpurwhypermean;
    purwhyperinvconc = inpurwhyperinvconc;
    poswhypermean = inposwhypermean;
    poswhyperinvconc = inposwhyperinvconc;
}

// 
// Matrices
//

void CodonM2aModel::UpdateNucMatrix()	{
    nucmatrix->CopyStationary(nucstat);
    nucmatrix->CorruptMatrix();
}

void CodonM2aModel::UpdateCodonMatrices()	{
    componentcodonmatrixarray->UpdateCodonMatrices();
}
    
void CodonM2aModel::UpdateMatrices()   {
    UpdateNucMatrix();
    UpdateCodonMatrices();
}

//
// Likelihood
//

double CodonM2aModel::GetLogLikelihood()	{
    // return GetIntegratedLogLikelihood();
    return phyloprocess->GetLogProb();
}

double CodonM2aModel::GetIntegratedLogLikelihood() {

    int ncat = 3;

    double total = 0;
    double logp[ncat];
    const vector<double>& w = componentomegaarray->GetWeights();
    double max = 0;
    for (int i=0; i<GetNsite(); i++) {
        int bkalloc = sitealloc->GetVal(i);

        for (int k=0; k<ncat; k++) {
            (*sitealloc)[i] = k;
            logp[k] = phyloprocess->SiteLogLikelihood(i);
            if ((!k) || (max<logp[k]))  {
                max = logp[k];
            }
        }

        double p = 0;
        for (int k=0; k<ncat; k++) {
            p += w[k] * exp(logp[k]-max);
        }
        double logl = log(p) + max;
        total += logl;

        (*sitealloc)[i] = bkalloc;
    }
    return total;
}

//
// Suff Stat and suffstatlogprobs
//

const PoissonSuffStatBranchArray* CodonM2aModel::GetLengthSuffStatArray()  {
    return lengthsuffstatarray;
}

double CodonM2aModel::LambdaHyperSuffStatLogProb() {
    return lambdasuffstat.GetLogProb(1.0,lambda);
}

const NucPathSuffStat& CodonM2aModel::GetNucPathSuffStat() {
    return nucpathsuffstat;
}

double CodonM2aModel::NucRatesSuffStatLogProb() {
    return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
}

double CodonM2aModel::PathSuffStatLogProb()	{
    return componentpathsuffstatarray->GetLogProb(*componentcodonmatrixarray);
}

double CodonM2aModel::OmegaSuffStatLogProb()    {
    componentomegaarray->SetParameters(purom,dposom+1,purw,posw);
    return componentomegaarray->GetPostProbArray(*siteomegasuffstatarray,sitepostprobarray);
}

//
// Priors
//

double CodonM2aModel::GetLogPrior() const {

    double total = 0;

    if (! FixedBranchLengths()) {
        total += BranchLengthsLogPrior();
    }
    if (! FixedNucRates())  {
        total += NucRatesLogPrior();
    }
    total += OmegaLogPrior();
    return total;
}

double CodonM2aModel::BranchLengthsLogPrior() const  {

    double total = 0;
    if (blmode == 0)    {
        total += LambdaHyperLogPrior();
    }
    total += branchlength->GetLogProb();
    return total;
}

double CodonM2aModel::LambdaHyperLogPrior()	 const {
    return -lambda / 10;
}

double CodonM2aModel::NucRatesLogPrior() const {

    double total = 0;
    total += Random::logDirichletDensity(nucrelrate,nucrelratehypercenter,1.0/nucrelratehyperinvconc);
    total += Random::logDirichletDensity(nucstat,nucstathypercenter,1.0/nucstathyperinvconc);
    return total;
}

//
// Hyper priors for omega mixture
//

double CodonM2aModel::OmegaLogPrior() const {
    double total = 0;
    total += PurOmegaLogProb();
    total += PosOmegaLogProb();
    total += PurWeightLogProb();
    total += PosWeightLogProb();
    return total;
}

// Beta prior for purifmean
double CodonM2aModel::PurOmegaLogProb()  const {
    double alpha = puromhypermean / puromhyperinvconc;
    double beta = (1-puromhypermean) / puromhyperinvconc;
    return Random::logBetaDensity(purom,alpha,beta);
}

// Gamma prior for dposom
double CodonM2aModel::PosOmegaLogProb() const {
    double alpha = 1.0 / dposomhyperinvshape;
    double beta = alpha / dposomhypermean;
    return Random::logGammaDensity(dposom,alpha,beta);
}

// Beta prior for purw
double CodonM2aModel::PurWeightLogProb() const {
    double alpha = purwhypermean / purwhyperinvconc;
    double beta = (1 - purwhypermean) / purwhyperinvconc;
    return Random::logBetaDensity(purw,alpha,beta);
}

// mixture of point mass at 0 (with prob pi) and Beta distribution (with prob 1 - pi) for posw
double CodonM2aModel::PosWeightLogProb() const {
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
double CodonM2aModel::PosSwitchLogProb() const {
    if (posw)   {
        return log(pi);
    }
    return log(1-pi);
}

//
//  Moves 
//

double CodonM2aModel::Move()	{

    ResampleSub(1.0);
    MoveParameters(30);
    return 1;
}

void CodonM2aModel::MoveParameters(int nrep)    {

    for (int rep=0; rep<nrep; rep++)	{

        if (! FixedBranchLengths()) {
            MoveBranchLengths();
        }

        CollectPathSuffStat();

        MoveOmega();

        if (! FixedNucRates())  {
            UpdateMatrices();
            MoveNucRates();
        }
    }
}

void CodonM2aModel::ResampleSub(double frac)  {
    UpdateMatrices();
    phyloprocess->Move(frac);
}

//
// Branch Lengths and hyperparam lambda
//

void CodonM2aModel::MoveBranchLengths()    {
        ResampleBranchLengths();
        if (blmode == 0)    {
            MoveLambda();
        }
}

void CodonM2aModel::ResampleBranchLengths()	{

    CollectLengthSuffStat();
    branchlength->GibbsResample(*lengthsuffstatarray);
}

void CodonM2aModel::CollectLengthSuffStat()    {

    lengthsuffstatarray->Clear();
    phyloprocess->AddLengthSuffStat(*lengthsuffstatarray);
}

void CodonM2aModel::MoveLambda()	{

    lambdasuffstat.Clear();
    branchlength->AddSuffStat(lambdasuffstat);

    ScalingMove(lambda,1.0,10,&CodonM2aModel::LambdaHyperLogProb,&CodonM2aModel::NoUpdate,this);
    ScalingMove(lambda,0.3,10,&CodonM2aModel::LambdaHyperLogProb,&CodonM2aModel::NoUpdate,this);

    blhypermean->SetAllBranches(1.0/lambda);
}

//
// Omega mixture 
//

void CodonM2aModel::CollectPathSuffStat()	{

    sitepathsuffstatarray->Clear();
    phyloprocess->AddPathSuffStat(*sitepathsuffstatarray);
}

void CodonM2aModel::CollectComponentPathSuffStat()	{

    componentpathsuffstatarray->Clear();
    sitepathsuffstatarray->AddToComponents(*componentpathsuffstatarray,*sitealloc);
}

void CodonM2aModel::MoveOmega() 	{

    CollectOmegaSuffStat();

    SlidingMove(purom,0.1,10,0,1,&CodonM2aModel::OmegaLogProb,&CodonM2aModel::NoUpdate,this);
    SlidingMove(purw,1.0,10,0,1,&CodonM2aModel::OmegaLogProb,&CodonM2aModel::NoUpdate,this);
    if (pi != 0)    {
        ScalingMove(dposom,1.0,10,&CodonM2aModel::OmegaLogProb,&CodonM2aModel::NoUpdate,this);
        SlidingMove(posw,1.0,10,0,1,&CodonM2aModel::OmegaLogProb,&CodonM2aModel::NoUpdate,this);
    }
    if ((pi != 0) && (pi != 1))    {
        SwitchPosWeight(10);
    }
    ResampleAlloc();
}

void CodonM2aModel::CollectOmegaSuffStat()	{

    siteomegasuffstatarray->Clear();
    siteomegasuffstatarray->AddSuffStat(*sitecodonmatrixarray,*sitepathsuffstatarray);
}

void CodonM2aModel::ResampleAlloc()	{
    OmegaSuffStatLogProb();
    sitealloc->GibbsResample(sitepostprobarray);
}

double CodonM2aModel::DrawBetaPosWeight()    {
    double alpha = poswhypermean / poswhyperinvconc;
    double beta = (1-poswhypermean) / poswhyperinvconc;
    return Random::BetaSample(alpha,beta);
}

double CodonM2aModel::SwitchPosWeight(int nrep)	{

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

//
// nucleotide parameters
//

void CodonM2aModel::CollectNucPathSuffStat()   {
    UpdateMatrices();
    nucpathsuffstat.Clear();
    nucpathsuffstat.AddSuffStat(*componentcodonmatrixarray,*componentpathsuffstatarray);
}

void CodonM2aModel::MoveNucRates()	{

    CollectComponentPathSuffStat();
    CollectNucPathSuffStat();

    ProfileMove(nucrelrate,0.1,1,3,&CodonM2aModel::NucRatesLogProb,&CodonM2aModel::UpdateNucMatrix,this);
    ProfileMove(nucrelrate,0.03,3,3,&CodonM2aModel::NucRatesLogProb,&CodonM2aModel::UpdateNucMatrix,this);
    ProfileMove(nucrelrate,0.01,3,3,&CodonM2aModel::NucRatesLogProb,&CodonM2aModel::UpdateNucMatrix,this);

    ProfileMove(nucstat,0.1,1,3,&CodonM2aModel::NucRatesLogProb,&CodonM2aModel::UpdateNucMatrix,this);
    ProfileMove(nucstat,0.01,1,3,&CodonM2aModel::NucRatesLogProb,&CodonM2aModel::UpdateNucMatrix,this);

    UpdateMatrices();
}

// summary statistics

double CodonM2aModel::GetMeanOmega() const {
    return posw*(1 + dposom) + (1-posw)*(purw*purom + (1-purw));
}

void CodonM2aModel::TraceHeader(std::ostream& os) {
    os << "#logprior\tlnL\tlength\t";
    os << "purom\tposom\tpurw\tposw\t";
    os << "statent\t";
    os << "rrent\n";
}

void CodonM2aModel::Trace(ostream& os)  {	
    os << GetLogPrior() << '\t';
    os << GetLogLikelihood() << '\t';
    os << branchlength->GetTotalLength() << '\t';
    os << purom << '\t' << dposom+1 << '\t' << purw << '\t' << posw << '\t';
    os << Random::GetEntropy(nucstat) << '\t';
    os << Random::GetEntropy(nucrelrate) << '\n';
    SubMatrix::diagerr = 0;
}

void CodonM2aModel::TracePostProb(ostream& os) {
    for (int i=0; i<GetNsite(); i++)    {
        os << sitepostprobarray[i][2] << '\t';
    }
    os << '\n';
}

void CodonM2aModel::GetSitesPostProb(double* array) const {
    for (int i=0; i<GetNsite(); i++)    {
        array[i] = sitepostprobarray[i][2];
        if (sitepostprobarray[i][2] < 0)    {
            cerr << "error in CodonM2aModel::GetSitesPostProb: negative prob\n";
            exit(1);
        }
    }
}

