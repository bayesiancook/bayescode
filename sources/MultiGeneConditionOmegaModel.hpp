
// this is a multigene version of singleomegamodel
//
// - branch lengths are shared across genes, and are iid Exponential of rate lambda
// - nucleotide relative exchangeabilities and stationaries are also shared across genes (uniform Dirichlet)
// - the array of gene-specific omega's are iid gamma with hyperparameters omegahypermean and omegahyperinvshape
//
// the sequence of MCMC moves is as follows:
// - genes resample substitution histories, gather path suff stats and move their omega's
// - master receives the array of omega's across genes, moves their hyperparameters and then broadcast the new value of these hyperparams
// - master collects branch length suff stats across genes, moves branch lengths and broadcasts their new value
// - master collects nuc path suffstats across genes, moves nuc rates and broadcasts their new value

#include "ConditionOmegaModel.hpp"
#include "Parallel.hpp"
#include "MultiGeneProbModel.hpp"

/**
 * \brief An Array of BranchArray's of OmegaPathSuffStat
 *
 * used in MultiGeneConditionOmegaModel, where each gene has a BranchArray of OmegaPathSuffStat.
 */

class OmegaPathSuffStatBidimArray : public Array<OmegaPathSuffStatArray>   {

    public:

    //! constructor, parameterized by underlying tree and size (number of genes)
    OmegaPathSuffStatBidimArray(int inncond, int insize): ncond(inncond), array(insize,(OmegaPathSuffStatArray*)0) {
        for (int i=0; i<GetSize(); i++) {
            array[i] = new OmegaPathSuffStatArray(ncond);
        }
    }

    ~OmegaPathSuffStatBidimArray() {
        for (int i=0; i<GetSize(); i++) {
            delete[] array[i];
        }
    }

    int GetSize() const override    {
        return array.size();
    }

    const OmegaPathSuffStatArray& GetVal(int i) const override {
        return *array[i];
    }

    OmegaPathSuffStatArray& operator[](int i) override {
        return *array[i];
    }

    //! clear all suff stats
    void Clear()    {
        for (int i=0; i<GetSize(); i++) {
            array[i]->Clear();
        }
    }

    private:

    int ncond;
    vector<OmegaPathSuffStatArray*> array;
};

class MultiGeneConditionOmegaModel : public MultiGeneProbModel {

    private:

	Tree* tree;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;
    int Ncond;
    int Nlevel;

    // branch lengths are shared across genes
    // iid expo (gamma of shape 1 and scale lambda)
    // where lambda is a hyperparameter
	double lambda;
	BranchIIDGamma* branchlength;
	
    // nucleotide rates are shared across genes
	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    double condvhypermean;
    double condvhyperinvshape;
    IIDGamma* condvarray;
    GammaSuffStat hypercondvsuffstat;

	double genewhypermean;
	double genewhyperinvshape;
	IIDGamma* genewarray;
    GammaSuffStat hypergenewsuffstat;

    double omegainvshape;

    ProductArray* meanomegabidimarray;
    ConditionSpecificMeanGammaBidimArray* condomegabidimarray;

    OmegaPathSuffStatBidimArray* omegapathsuffstatbidimarray;

    // suffstats for paths, as a function of branch lengths
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
    // suff stats for branch lengths, as a function of lambda
	GammaSuffStat hyperlengthsuffstat;

    // suffstats for paths, as a function of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // each gene defines its own ConditionOmegaModel
    std::vector<ConditionOmegaModel*> geneprocess;

    // total log likelihood (summed across all genes)
    double lnL;
    // total logprior for gene-specific variables (here, omega only)
    // summed over all genes
    double GeneLogPrior;

    public:

    //-------------------
    // Construction and allocation
    //-------------------

    MultiGeneConditionOmegaModel(string datafile, string intreefile, int inNcond, int inNlevel, int inmyid, int innprocs) : MultiGeneProbModel(inmyid,innprocs) {

        AllocateAlignments(datafile);
        treefile = intreefile;
        Ncond = inNcond;
        Nlevel = inNlevel;

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

    void Allocate() {

        lambda = 10;
        branchlength = new BranchIIDGamma(*tree,1.0,lambda);
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

        nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));

        nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        condvhypermean = 1.0;
        condvhyperinvshape = 1.0;
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condvarray = new IIDGamma(Ncond,alpha,beta);
        for (int j=0; j<Ncond; j++)  {
            (*condvarray)[j] = 1.0;
        }

        genewhypermean = 1.0;
        genewhyperinvshape = 1.0;
        double genealpha = 1.0 / genewhyperinvshape;
        double genebeta = genealpha / genewhypermean;
        genewarray = new IIDGamma(Ngene,genealpha,genebeta);

        meanomegabidimarray = new ProductArray(*condvarray,*genewarray);
        omegainvshape = 1.0;
        condomegabidimarray = new ConditionSpecificMeanGammaBidimArray(*meanomegabidimarray,omegainvshape);

        // should be a branch site structure
        omegapathsuffstatbidimarray = new OmegaPathSuffStatBidimArray(Ncond,Ngene);

        lnL = 0;
        GeneLogPrior = 0;

        if (! GetMyid())    {
            geneprocess.assign(0,(ConditionOmegaModel*) 0);
        }
        else    {
            geneprocess.assign(GetLocalNgene(),(ConditionOmegaModel*) 0);

            for (int gene=0; gene<GetLocalNgene(); gene++)   {
                geneprocess[gene] = new ConditionOmegaModel(GetLocalGeneName(gene),treefile,Ncond,Nlevel);
                geneprocess[gene]->Allocate();
            }
        }
    }

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    int GetNbranch() const {
        return tree->GetNbranch();
    }

    const Tree* GetTree() const {
        return tree;
    }

    int GetNcond() const    {
        return Ncond;
    }

    //-------------------
    // Traces and Monitors
    //-------------------

    void TraceHeader(ostream& os) const {

        os << "#logprior\tlnL\tlength\t";
        os << "genemean\tinvshape\t";
        os << "condinvshape\t";
        os << "omegainvshape\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream& os) const {
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << genewhypermean << '\t' << genewhyperinvshape << '\t';
        os << condvhyperinvshape << '\t';
        os << omegainvshape << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';
		os.flush();
    }

    void PrintGeneEffects(ostream& os) const    {
        os << *genewarray;
        os << '\n';
        os.flush();
    }

    void PrintCondEffects(ostream& os) const  {
        os << *condvarray;
        os << '\n';
        os.flush();
    }

    void PrintDeviations(ostream& os) const {
        for (int j=0; j<GetNcond(); j++)  {
            for (int i=0; i<GetNgene(); i++)    {
                os << GetOmega(i,j) / GetMeanOmega(i,j) << '\t';
            }
        }
        os << '\n';
        os.flush();
    }

	void Monitor(ostream& os) const {}

	void MasterToStream(ostream& os) const {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << condvhypermean << '\t' << condvhyperinvshape << '\n';
        os << *condvarray << '\n';
        os << genewhypermean << '\t' << genewhyperinvshape << '\n';
        os << *genewarray << '\n';
        os << omegainvshape << '\n';
        os << *condomegabidimarray << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }

	void MasterFromStream(istream& is) {
        is >> lambda;
        is >> *branchlength;
        is >> condvhypermean >> condvhyperinvshape;
        is >> *condvarray;
        is >> genewhypermean >> genewhyperinvshape;
        is >> *genewarray;
        is >> omegainvshape;
        is >> *condomegabidimarray;
        is >> nucrelrate;
        is >> nucstat;
    }

    //-------------------
    // Updates
    //-------------------

    void FastUpdate()   {
        branchlength->SetScale(lambda);
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condvarray->SetShape(alpha);
        condvarray->SetScale(beta);
        double genealpha = 1.0 / genewhyperinvshape;
        double genebeta = genealpha / genewhypermean;
        genewarray->SetShape(genealpha);
        genewarray->SetScale(genebeta);
        condomegabidimarray->SetInvShape(omegainvshape);
    }

    void MasterUpdate() override {

        FastUpdate();

        if (nprocs > 1) {
            MasterSendGlobalBranchLengths();
            MasterSendGlobalNucRates();
            MasterSendOmegaHyperParameters();
            MasterSendOmega();
            MasterReceiveLogProbs();
        }
    }

    void SlaveUpdate() override {
        SlaveReceiveGlobalBranchLengths();
        SlaveReceiveGlobalNucRates();
        SlaveReceiveOmegaHyperParameters();
        SlaveReceiveOmega();
        GeneUpdate();
        SlaveSendLogProbs();
    }

    void GeneUpdate()	{
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->Update();
        }
    }

	void TouchNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

    void NoUpdate() {}

    //-------------------
    // Log Prior and Likelihood
    //-------------------
    

    double GetOmega(int gene, int cond) const {
        return condomegabidimarray->GetVal(gene).GetVal(cond);
    }

    double GetMeanOmega(int gene, int cond) const {
        return meanomegabidimarray->GetVal(gene).GetVal(cond);
    }

    double GetLogPrior() const {
		double total = 0;
		total += BranchLengthsHyperLogPrior();
		total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += CondVHyperLogPrior();
        total += CondVLogPrior();
        total += GeneWHyperLogPrior();
        total += GeneWLogPrior();
        total += OmegaInvShapeLogPrior();
		total += OmegaLogPrior();
		return total;
    }

    double OmegaLogPrior() const    {
        meanomegabidimarray->Update();
        condomegabidimarray->SetInvShape(omegainvshape);
        return condomegabidimarray->GetLogProb();
    }

	double BranchLengthsHyperLogPrior() const {
		return -lambda / 10;
	}

	double BranchLengthsLogPrior() const {
		return branchlength->GetLogProb();
	}

    double NucRatesLogPrior() const {
        return 0;
    }

    double CondVHyperLogPrior() const {
        return -condvhypermean - condvhyperinvshape;
    }

    double GeneWHyperLogPrior() const {
        return -genewhypermean - genewhyperinvshape;
    }

    double CondVLogPrior() const {
        return condvarray->GetLogProb();
    }

    double GeneWLogPrior() const {
        return genewarray->GetLogProb();
    }

    double OmegaInvShapeLogPrior() const   {
        return omegainvshape;
    }

    double GetLogLikelihood() const {
        return lnL;
    }

    //-------------------
    // Suff Stat Log Probs
    //-------------------

    // suff stat for moving branch lengths hyperparameter (lambda)
	double BranchLengthsHyperSuffStatLogProb() const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    // suff stats for moving nuc rates
    double NucRatesSuffStatLogProb() const {
        return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    }

    double CondVHyperSuffStatLogProb() const  {
        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        return hypercondvsuffstat.GetLogProb(alpha,beta);
    }

    double GeneWHyperSuffStatLogProb() const    {
        double alpha = 1.0 / genewhyperinvshape;
        double beta = alpha / genewhypermean;
        return hypergenewsuffstat.GetLogProb(alpha,beta);
    }

    double OmegaSuffStatLogProb() const {
        double total = 0;
        for (int i=0; i<GetNgene(); i++)    {
            total += GeneOmegaSuffStatLogProb(i);
        }
        return total;
    }

    double GeneOmegaSuffStatLogProb(int gene) const {
        double total = 0;
        for (int j=0; j<GetNcond(); j++)  {
            total += GeneConditionOmegaSuffStatLogProb(gene,j);
        }
        return total;
    }

    double ConditionOmegaSuffStatLogProb(int cond) const {
        double total = 0;
        for (int i=0; i<GetNgene(); i++) {
            total += GeneConditionOmegaSuffStatLogProb(i,cond);
        }
        return total;
    }

    double GeneConditionOmegaSuffStatLogProb(int gene, int cond) const {

        const OmegaPathSuffStat& suffstat = omegapathsuffstatbidimarray->GetVal(gene).GetVal(cond);
        int count = suffstat.GetCount();
        double b = suffstat.GetBeta();

        double mean = genewarray->GetVal(gene) * condvarray->GetVal(cond);
        double alpha = 1.0 / omegainvshape;
        double beta = alpha / mean;

        return alpha*log(beta) - Random::logGamma(alpha) + Random::logGamma(alpha + count) - (alpha+count)*log(beta+b);
    }

    //-------------------
    // Log Probs for MH moves
    //-------------------

    // log prob for moving branch lengths hyperparameter (lambda)
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + NucRatesSuffStatLogProb();
    }

    // log prob for moving branch v hyperparameters 
    double CondVHyperLogProb() const {
        return CondVHyperLogPrior() + CondVHyperSuffStatLogProb();
    }

    // log prob for moving gene w hyperparameters
    double GeneWHyperLogProb() const {
        return GeneWHyperLogPrior() + GeneWHyperSuffStatLogProb();
    }

    //-------------------
    // Moves
    //-------------------

    // all methods starting with Master are called only by master
    // for each such method, there is a corresponding method called by slave, and starting with Slave
    //
    // all methods starting with Gene are called only be slaves, and do some work across all genes allocated to that slave

    void MasterMove() override {

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            MasterReceiveOmegaSuffStat();
            MoveOmegaHyperParameters(3);
            MasterSendOmegaHyperParameters();

            MasterReceiveLengthSuffStat();
            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();
            MasterSendGlobalBranchLengths();

            MasterReceiveNucPathSuffStat();
            MoveNucRates();
            MasterSendGlobalNucRates();
        }

        MasterReceiveOmega();
        MasterReceiveLogProbs();
    }

    // slave move
    void SlaveMove() override {

        GeneResampleSub(1.0);

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

            GeneCollectPathSuffStat();
            SlaveSendOmegaSuffStat();
            SlaveReceiveOmegaHyperParameters();
            GeneResampleOmega();

            SlaveSendLengthSuffStat();
            SlaveReceiveGlobalBranchLengths();

            SlaveSendNucPathSuffStat();
            SlaveReceiveGlobalNucRates();
        }

        SlaveSendOmega();
        SlaveSendLogProbs();
    }

    void GeneResampleSub(double frac)  {

        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleSub(frac);
        }
    }

    void GeneCollectPathSuffStat()  {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectPathSuffStat();
        }
    }

    void GeneResampleOmega()  {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->ResampleOmega();
        }
    }

    void ResampleBranchLengths()    {
		branchlength->GibbsResample(*lengthpathsuffstatarray);
    }

	void MoveBranchLengthsHyperParameter()	{

		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);

        ScalingMove(lambda,1.0,10,&MultiGeneConditionOmegaModel::BranchLengthsHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&MultiGeneConditionOmegaModel::BranchLengthsHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);

		branchlength->SetScale(lambda);
    }

    void MoveNucRates()    {

        ProfileMove(nucrelrate,0.1,1,3,&MultiGeneConditionOmegaModel::NucRatesLogProb,&MultiGeneConditionOmegaModel::TouchNucMatrix,this);
        ProfileMove(nucrelrate,0.03,3,3,&MultiGeneConditionOmegaModel::NucRatesLogProb,&MultiGeneConditionOmegaModel::TouchNucMatrix,this);
        ProfileMove(nucrelrate,0.01,3,3,&MultiGeneConditionOmegaModel::NucRatesLogProb,&MultiGeneConditionOmegaModel::TouchNucMatrix,this);

        ProfileMove(nucstat,0.1,1,3,&MultiGeneConditionOmegaModel::NucRatesLogProb,&MultiGeneConditionOmegaModel::TouchNucMatrix,this);
        ProfileMove(nucstat,0.01,1,3,&MultiGeneConditionOmegaModel::NucRatesLogProb,&MultiGeneConditionOmegaModel::TouchNucMatrix,this);
    }

    /*
    double BranchGeneCompMove(double tuning, int nrep)   {

        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            double deltalogprob = -CondVLogPrior() - GeneWLogPrior();

            double m = tuning*(Random::Uniform() - 0.5);
            double e = exp(m);
            for (int i=0; i<GetNgene(); i++) {
                (*genewarray)[i] *= e;
            }
            for (int j=0; j<GetNbranch(); j++)  {
                (*condvarray)[j] /= e;
            }
            
            deltalogprob += CondVLogPrior() + GeneWLogPrior();
            deltalogprob += (GetNgene() - GetNbranch())*m;

            int acc = (log(Random::Uniform()) < deltalogprob);
            if (acc)  {
                nacc++;
            }
            else    {
                for (int i=0; i<GetNgene(); i++) {
                    (*genewarray)[i] /= e;
                }
                for (int j=0; j<GetNbranch(); j++)  {
                    (*condvarray)[j] *= e;
                }
            }
        }
        return nacc/nrep;
    }
    */

    void MoveOmegaHyperParameters(int nrep)  {
        for (int rep=0; rep<nrep; rep++)    {
            MoveGeneW(1.0,1);
            MoveCondV(1.0,1);
            MoveGeneW(0.3,1);
            MoveCondV(0.3,1);
            /*
            BranchGeneCompMove(1.0,1);
            BranchGeneCompMove(0.3,1);
            */
            MoveOmegaInvShape(1.0,1);
            MoveOmegaInvShape(0.3,1);
        }
        // MoveCondVHyperParams(1.0,100);
        MoveGeneWHyperParams(1.0,100);
    }

    double MoveCondV(double tuning, int nrep) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int j=1; j<GetNcond(); j++)  {
                double deltalogprob = - condvarray->GetLogProb(j);
                deltalogprob -= ConditionOmegaSuffStatLogProb(j);
                double m = tuning*(Random::Uniform() - 0.5);
                double e = exp(m);
                (*condvarray)[j] *= e;
                deltalogprob += condvarray->GetLogProb(j);
                deltalogprob += ConditionOmegaSuffStatLogProb(j);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc)    {
                    nacc++;
                }
                else    {
                    (*condvarray)[j] /= e;
                }
            }
        }
        return nacc / GetNbranch() / nrep;
    }

    double MoveGeneW(double tuning, int nrep)   {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<GetNgene(); i++)    {
                double deltalogprior = - genewarray->GetLogProb(i);
                double deltasuffstatlogprob = - GeneOmegaSuffStatLogProb(i);
                double m = tuning*(Random::Uniform() - 0.5);
                double e = exp(m);
                (*genewarray)[i] *= e;
                deltalogprior += genewarray->GetLogProb(i);
                deltasuffstatlogprob += GeneOmegaSuffStatLogProb(i);
                double deltalogprob = deltalogprior + deltasuffstatlogprob + m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc)    {
                    nacc++;
                }
                else    {
                    (*genewarray)[i] /= e;
                }
            }
        }
        return nacc / GetNgene() / nrep;
    }

    double MoveOmegaInvShape(double tuning, int nrep)   {

        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            double deltalogprob = -OmegaInvShapeLogPrior() - OmegaSuffStatLogProb();
            double m = tuning*(Random::Uniform() - 0.5);
            double e = exp(m);
            omegainvshape *= e;
            deltalogprob += OmegaInvShapeLogPrior() + OmegaSuffStatLogProb();
            deltalogprob += m;
            int acc = (log(Random::Uniform()) < deltalogprob);
            if (acc)    {
                nacc++;
            }
            else    {
                omegainvshape /= e;
            }
        }
        condomegabidimarray->SetInvShape(omegainvshape);
        return nacc/nrep;
    }

    void MoveCondVHyperParams(double tuning, int nrep)  {

        hypercondvsuffstat.Clear();
        hypercondvsuffstat.AddSuffStat(*condvarray);

        ScalingMove(condvhyperinvshape,1.0,10,&MultiGeneConditionOmegaModel::CondVHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);
        ScalingMove(condvhyperinvshape,0.3,10,&MultiGeneConditionOmegaModel::CondVHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);

        double alpha = 1.0 / condvhyperinvshape;
        double beta = alpha / condvhypermean;
        condvarray->SetShape(alpha);
        condvarray->SetScale(beta);
    }

    void MoveGeneWHyperParams(double tuning, int nrep)  {

        hypergenewsuffstat.Clear();
        hypergenewsuffstat.AddSuffStat(*genewarray);

        ScalingMove(genewhypermean,1.0,10,&MultiGeneConditionOmegaModel::GeneWHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);
        ScalingMove(genewhypermean,0.3,10,&MultiGeneConditionOmegaModel::GeneWHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);
        ScalingMove(genewhyperinvshape,1.0,10,&MultiGeneConditionOmegaModel::GeneWHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);
        ScalingMove(genewhyperinvshape,0.3,10,&MultiGeneConditionOmegaModel::GeneWHyperLogProb,&MultiGeneConditionOmegaModel::NoUpdate,this);

        double alpha = 1.0 / genewhyperinvshape;
        double beta = alpha / genewhypermean;
        genewarray->SetShape(alpha);
        genewarray->SetScale(beta);
    }

    //-------------------
    // MPI send / receive
    //-------------------

    // global branch lengths
    
    void MasterSendGlobalBranchLengths() {
        MasterSendGlobal(*branchlength);
    }

    void SlaveReceiveGlobalBranchLengths()   {
        SlaveReceiveGlobal(*branchlength);
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetBranchLengths(*branchlength);
        }
    }

    // global nuc rates

    void MasterSendGlobalNucRates()   {
        MasterSendGlobal(nucrelrate,nucstat);
    }

    void SlaveReceiveGlobalNucRates()   {

        SlaveReceiveGlobal(nucrelrate,nucstat);
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetNucRates(nucrelrate,nucstat);
        }
    }

    // omega arrays

    void SlaveSendOmega()   {
        for (int gene=0; gene<GetLocalNgene(); gene++)  {
            (*condomegabidimarray)[gene].Copy(*geneprocess[gene]->GetOmegaArray());
        }
        SlaveSendGeneArray(*condomegabidimarray);
    }

    void MasterReceiveOmega()    {
        MasterReceiveGeneArray(*condomegabidimarray);
    }

    void MasterSendOmega()  {
        MasterSendGeneArray(*condomegabidimarray);
    }
    
    void SlaveReceiveOmega()    {
        SlaveReceiveGeneArray(*condomegabidimarray);
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->SetOmegaTree(condomegabidimarray->GetVal(gene));
        }
    }

    // omega hyperparameters

    void MasterSendOmegaHyperParameters()   {
        MasterSendGlobal(condvhypermean,condvhyperinvshape);
        MasterSendGlobal(genewhypermean,genewhyperinvshape);
        MasterSendGlobal(omegainvshape);
        MasterSendGlobal(*condvarray);
        MasterSendGlobal(*genewarray);
    }

    void SlaveReceiveOmegaHyperParameters() {
        SlaveReceiveGlobal(condvhypermean,condvhyperinvshape);
        SlaveReceiveGlobal(genewhypermean,genewhyperinvshape);
        SlaveReceiveGlobal(omegainvshape);
        SlaveReceiveGlobal(*condvarray);
        SlaveReceiveGlobal(*genewarray);
        for (int gene=0; gene<GetLocalNgene(); gene++)    {
            geneprocess[gene]->SetCondVHyperParams(condvhypermean,condvhyperinvshape);
            geneprocess[gene]->SetCondV(*condvarray);
            geneprocess[gene]->SetGeneW(genewarray->GetVal(gene));
            geneprocess[gene]->SetOmegaHyperInvShape(omegainvshape);
        }
    }

    // branch length suff stat

    void SlaveSendLengthSuffStat()  {
        lengthpathsuffstatarray->Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectLengthSuffStat();
            lengthpathsuffstatarray->Add(*geneprocess[gene]->GetLengthPathSuffStatArray());
        }
        SlaveSendAdditive(*lengthpathsuffstatarray);
    }

    void MasterReceiveLengthSuffStat()  {

        lengthpathsuffstatarray->Clear();
        MasterReceiveAdditive(*lengthpathsuffstatarray);
    }

    // nuc path suff stat

    void SlaveSendNucPathSuffStat()  {
        nucpathsuffstat.Clear();
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectNucPathSuffStat();
            nucpathsuffstat += geneprocess[gene]->GetNucPathSuffStat();
        }
        SlaveSendAdditive(nucpathsuffstat);
    }

    void MasterReceiveNucPathSuffStat()  {
        nucpathsuffstat.Clear();
        MasterReceiveAdditive(nucpathsuffstat);
    }

    // omega path suff stat

    void SlaveSendOmegaSuffStat()   {
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            geneprocess[gene]->CollectOmegaSuffStat();
            (*omegapathsuffstatbidimarray)[gene].Array<OmegaPathSuffStat>::Copy(*geneprocess[gene]->GetOmegaPathSuffStatArray());
        }
        SlaveSendGeneArray(*omegapathsuffstatbidimarray);
    }

    void MasterReceiveOmegaSuffStat()   {
        MasterReceiveGeneArray(*omegapathsuffstatbidimarray);
    }

    // log probs

    void SlaveSendLogProbs()   {

        GeneLogPrior = 0;
        lnL = 0;
        for (int gene=0; gene<GetLocalNgene(); gene++)   {
            GeneLogPrior += geneprocess[gene]->GetLogPrior();
            lnL += geneprocess[gene]->GetLogLikelihood();
        }
        SlaveSendAdditive(GeneLogPrior);
        SlaveSendAdditive(lnL);
    }

    void MasterReceiveLogProbs()    {

        GeneLogPrior = 0;
        MasterReceiveAdditive(GeneLogPrior);
        lnL = 0;
        MasterReceiveAdditive(lnL);
    }
};

