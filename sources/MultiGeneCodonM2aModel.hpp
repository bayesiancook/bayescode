
#include "CodonM2aModel.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "IIDDirichlet.hpp"

class MultiGeneCodonM2aModel : public MultiGeneMPIModule	{

    public:

    //-------------------
    // Constructors
    // ------------------

    MultiGeneCodonM2aModel(string datafile, string intreefile, double inpihypermean, double inpihyperinvconc, int inmyid, int innprocs);
    void Allocate();
    void Unfold();

    //-------------------
    // Accessors
    // ------------------

	CodonStateSpace* GetCodonStateSpace()   {
		return (CodonStateSpace*) refcodondata->GetStateSpace();
	}

    //-------------------
    // Setting and updating
    // ------------------

    void SetAcrossGenesModes(int inblmode, int innucmode, int inpurommode, int indposommode, int inpurwmode, int inposwmode);

    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape, double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc);

	void UpdateNucMatrix();
    void SetMixtureArrays();

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream& os);
    void Trace(ostream& os);

    void TracePosWeight(ostream& os);
    void TracePosOm(ostream& os);

    void MasterTraceSitesPostProb(ostream& os);
    void SlaveTraceSitesPostProb();

	void Monitor(ostream& os) {}
	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

    int GetNpos();
	double GetMeanTotalLength();
    double GetMeanLength();
    double GetVarLength();
    double GetVarNucRelRate();
    double GetVarNucStat();

    //-------------------
    // Log Probs
    // ------------------

    // priors
    double GetLogPrior();

    double GlobalBranchLengthsLogPrior();
    double GeneBranchLengthsHyperLogPrior();

    double LambdaHyperLogPrior();
    double BranchLengthsHyperInvShapeLogPrior();

    double GlobalNucRatesLogPrior();
    double GeneNucRatesHyperLogPrior();

    double MixtureHyperLogPrior();

    // likelihood
    double GetLogLikelihood()   {
        return lnL;
    }

    // hyper suffstatlogprob
	double LambdaHyperSuffStatLogProb();
    double BranchLengthsHyperSuffStatLogProb();
    double NucRatesHyperSuffStatLogProb();
    double MixtureHyperSuffStatLogProb();

    // suffstatlogprob
    double NucRatesSuffStatLogProb();

    //-------------------
    // Moves
    // ------------------

    // general move schedule
    void MasterMove();
    void SlaveMove();

    void GeneResampleSub(double frac);
    void MoveGeneParameters(int nrep);

    void ResampleBranchLengths();
	void MoveLambda();
	double MoveLambda(double tuning, int nrep);

    void MoveBranchLengthsHyperParameters();
    double BranchLengthsHyperScalingMove(double tuning, int nrep);
    double BranchLengthsHyperInvShapeMove(double tuning, int nrep);

    /*
    void SetNucRelRateCenterToMean();
    void SetNucStatCenterToMean();
    */
    void MoveNucRatesHyperParameters();

    void MoveMixtureHyperParameters() ;

	double NucRatesHyperProfileMove(vector<double>& x, double tuning, int n, int nrep);
	double NucRatesHyperScalingMove(double& x, double tuning, int nrep);

    void ResamplePi();

	double MixtureHyperSlidingMove(double& x, double tuning, int nrep, double min = 0, double max = 0);
	double MixtureHyperScalingMove(double& x, double tuning, int nrep);

    void MoveNucRates();

	double MoveRR(double tuning, int n, int nrep);
	double MoveNucStat(double tuning, int n, int nrep);

    //-------------------
    // MPI send/receive
    // ------------------

    // global parameters

    void MasterSendGlobalBranchLengths();
    void SlaveReceiveGlobalBranchLengths();

    void MasterSendGlobalNucRates();
    void SlaveReceiveGlobalNucRates();

    // gene-specific parameters

    void MasterSendGeneBranchLengths();
    void SlaveReceiveGeneBranchLengths();
    void SlaveSendGeneBranchLengths();
    void MasterReceiveGeneBranchLengths();

    void MasterSendGeneNucRates();
    void SlaveReceiveGeneNucRates();
    void SlaveSendGeneNucRates();
    void MasterReceiveGeneNucRates();

    void SlaveSendMixture();
    void MasterReceiveMixture();
    void SlaveReceiveMixture();
    void MasterSendMixture();

    // hyper parameters

    void MasterSendBranchLengthsHyperParameters();
    void SlaveReceiveBranchLengthsHyperParameters();

    void MasterSendNucRatesHyperParameters();
    void SlaveReceiveNucRatesHyperParameters();

    void MasterSendMixtureHyperParameters();
    void SlaveReceiveMixtureHyperParameters();

    // suff stats

    // branch lengths
    void SlaveSendBranchLengthsSuffStat();
    void MasterReceiveBranchLengthsSuffStat();
    
    void SlaveSendNucPathSuffStat();
    void MasterReceiveNucPathSuffStat();

    // hyper suffstats

    void SlaveSendBranchLengthsHyperSuffStat();
    void MasterReceiveBranchLengthsHyperSuffStat();

    void SlaveSendNucRatesHyperSuffStat();
    void MasterReceiveNucRatesHyperSuffStat();

    void SlaveSendMixtureHyperSuffStat();
    void MasterReceiveMixtureHyperSuffStat();

    // log likelihoods
    
    void SlaveSendLogProbs();
    void MasterReceiveLogProbs();

    private:

	Tree* tree;
	CodonSequenceAlignment* refcodondata;
	const TaxonSet* taxonset;

    string treefile;

	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	GammaSuffStat lambdasuffstat;
	
    double blhyperinvshape;
    GammaWhiteNoiseArray* branchlengtharray;
	PoissonSuffStatBranchArray* lengthsuffstatarray;
    GammaSuffStatBranchArray* lengthhypersuffstatarray;

    vector<double> mixhyperparam;

    double& puromhypermean;
	double& puromhyperinvconc;
	IIDBeta* puromarray;
	BetaSuffStat puromsuffstat;

    double& dposomhypermean;
    double& dposomhyperinvshape;
    IIDGamma* dposomarray;
    GammaSuffStat dposomsuffstat;

    double& purwhypermean;
    double& purwhyperinvconc;
    IIDBeta* purwarray;
    BetaSuffStat purwsuffstat;

    double& poswhypermean;
    double& poswhyperinvconc;
    IIDBernoulliBeta* poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double& pi;

    // shared nuc rates
	GTRSubMatrix* nucmatrix;
    NucPathSuffStat nucpathsuffstat;

    // gene-specific nuc rates
    vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    IIDDirichlet* nucrelratearray;
    DirichletSuffStat nucrelratesuffstat;

    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    IIDDirichlet* nucstatarray;
    DirichletSuffStat nucstatsuffstat;

    std::vector<CodonM2aModel*> geneprocess;

    double lnL;
    double GeneLogPrior;

    int burnin;

    // 0: free (fixed hyper parameters)
    // 1: free and shrinkage (free hyper parameters)
    // 2: shared across genes
    int blmode;
    int nucmode;
    int purommode;
    int dposommode;
    int purwmode;
    int poswmode;

    /*
    double nucrracc1, nucrracc2, nucrracc3, nucrrtot1, nucrrtot2, nucrrtot3;
    double nucstatacc1, nucstatacc2, nucstatacc3, nucstattot1, nucstattot2, nucstattot3;
    */
};

