
#include "CodonM2aModel.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "IIDDirichlet.hpp"

#include "Chrono.hpp"


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
    /*
    void SetBranchLengthsHyperParameters();
    void SetNucRateHyperParameters();
    */

    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape, double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean, double inposwhyperinvconc);

	void UpdateNucMatrix();
    void SlaveSetMixtureArrays();

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream& os);
    void Trace(ostream& os);

    void TracePosWeight(ostream& os);
    void TracePosOm(ostream& os);
    /*
    void SlaveTracePostProbHeader(string name);
    void SlaveTracePostProb(string name);
    */

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

    double LambdaHyperLogPrior();
    double BranchLengthsHyperInvShapeLogPrior();

    double BranchLengthsHyperLogPrior();
    double BranchLengthsLogPrior();

    double NucRatesHyperLogPrior();
    double NucRatesLogPrior();

    double MixtureHyperLogPrior();
    double MixtureLogPrior();

    // likelihood
    double GetLogLikelihood();

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

    void SlaveResampleSub();

    void MasterResampleBranchLengths();
	void MasterMoveLambda();
	double MoveLambda(double tuning, int nrep);

    void MasterMoveBranchLengthsHyperParameters();
    double BranchLengthsHyperScalingMove(double tuning, int nrep);
    double BranchLengthsHyperInvShapeMove(double tuning, int nrep);

    void SlaveMoveBranchLengths();

    void MasterMoveNucRatesHyperParameters();
    void SlaveMoveNucRates();

    void MasterMoveMixtureHyperParameters() ;

	double NucRatesHyperProfileMove(vector<double>& x, double tuning, int n, int nrep);
	double NucRatesHyperScalingMove(double& x, double tuning, int nrep);

    void ResamplePi();

	double MixtureHyperSlidingMove(double& x, double tuning, int nrep, double min = 0, double max = 0);
	double MixtureHyperScalingMove(double& x, double tuning, int nrep);

    void MasterMoveNucRates();

	double MoveRR(double tuning, int n, int nrep);
	double MoveNucStat(double tuning, int n, int nrep);

    void SlaveMoveOmega();

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
    
    void SlaveCollectPathSuffStat();

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
    
    void SlaveSendLogLikelihood();
    void MasterReceiveLogLikelihood();

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
    vector<GammaWhiteNoise*> branchlengtharray;
	PoissonSuffStatBranchArray* lengthsuffstatarray;
    GammaSuffStatBranchArray* lengthhypersuffstatarray;

    double puromhypermean;
	double puromhyperinvconc;
	IIDBeta* puromarray;
	BetaSuffStat puromsuffstat;

    double dposomhypermean;
    double dposomhyperinvshape;
    IIDGamma* dposomarray;
    GammaSuffStat dposomsuffstat;

    double purwhypermean;
    double purwhyperinvconc;
    IIDBeta* purwarray;
    BetaSuffStat purwsuffstat;

    double poswhypermean;
    double poswhyperinvconc;
    IIDBernoulliBeta* poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double pi;

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

    double totlnL;
    double* lnL;

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

    Chrono timepercycle;
    Chrono omegachrono,hyperchrono,mastersampling;
};

