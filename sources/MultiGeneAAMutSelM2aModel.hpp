
#include "AAMutSelM2aModel.hpp"
#include "IIDBernoulliBeta.hpp"
#include "IIDBeta.hpp"
#include "IIDDirichlet.hpp"
#include "MultiGeneMPIModule.hpp"
#include "Parallel.hpp"
#include "ProbModel.hpp"

class MultiGeneAAMutSelM2aModel : public MultiGeneMPIModule, public ProbModel {
  public:
    //-------------------
    // Constructors
    // ------------------

    MultiGeneAAMutSelM2aModel(string datafile, string intreefile, double inpihypermean,
                              double inpihyperinvconc, int inmyid, int innprocs);
    void Allocate();
    void Unfold();

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
    // void SetAcrossGenesModes(int inblmode, int innucmode, int inpurommode, int
    // indposommode, int inpurwmode, int inposwmode);
    void SetAcrossGenesModes(int inblmode, int indposommode, int inposwmode);

    // void SetMixtureHyperParameters(double inpuromhypermean, double
    // inpuromhyperinvconc, double indposomhypermean, double indposomhyperinvshape,
    // double inpurwhypermean, double inpurwhyperinvconc, double inposwhypermean,
    // double inposwhyperinvconc);
    void SetMixtureHyperParameters(double indposomhypermean, double indposomhyperinvshape,
                                   double inposwhypermean, double inposwhyperinvconc);

    // void UpdateNucMatrix();
    void SetMixtureArrays();

    void NoUpdate() {}

    //-------------------
    // Traces and Monitors
    // ------------------

    void TraceHeader(ostream &os) const;
    void Trace(ostream &os) const;

    void TracePosWeight(ostream &os) const;
    void TracePosOm(ostream &os) const;

    void MasterTraceSitesPostProb(ostream &os);
    void SlaveTraceSitesPostProb();

    void Monitor(ostream &os) const {}
    void FromStream(istream &is) {}
    void ToStream(ostream &os) const {}

    // summary statistics for tracing MCMC
    int GetNpos() const;
    double GetMeanTotalLength() const;
    double GetMeanLength() const;
    double GetVarLength() const;
    // double GetVarNucRelRate() const;
    // double GetVarNucStat() const;

    //-------------------
    // Log Probs
    // ------------------

    //-------------------
    // Log Priors
    // ------------------

    // total log prior
    double GetLogPrior() const;

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
    // double GlobalNucRatesLogPrior() const {
    //    return nucrelratearray->GetLogProb() + nucstatarray->GetLogProb();
    //}

    // exponential of mean 1 for nucrelrate and nucstat hyper inverse
    // concentration
    /*double GeneNucRatesHyperLogPrior() const {
        double total = 0;
        if (nucmode == 1)   {
            total -= nucrelratehyperinvconc;
            total -= nucstathyperinvconc;
        }
        return total;
    }*/

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
        // total -= puromhyperinvconc;
        // total -= purwhyperinvconc;
        // exponential of mean 10 for poswhyperinvconc
        total -= 10 * poswhyperinvconc;
        // exponential of mean 1 for dposomhypermean
        total -= dposomhypermean;
        // exponential of mean 10 for dposomhyperinvshape
        total -= 10 * dposomhyperinvshape;
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
    // double NucRatesSuffStatLogProb() const {
    //    return nucpathsuffstat.GetLogProb(*nucmatrix,*GetCodonStateSpace());
    //}

    // suff stat for gene-specific nuc rates, as a function of nucrate
    // hyperparameters
    // double NucRatesHyperSuffStatLogProb() const {
    //    double total = 0;
    //    total +=
    //    nucrelratesuffstat.GetLogProb(nucrelratehypercenter,1.0/nucrelratehyperinvconc);
    //    total +=
    //    nucstatsuffstat.GetLogProb(nucstathypercenter,1.0/nucstathyperinvconc);
    //    return total;
    //}

    // suff stat for gene-specific mixture parameters, as a function of mixture
    // hyperparameters
    double MixtureHyperSuffStatLogProb() const {
        double total = 0;

        // double puromalpha = puromhypermean / puromhyperinvconc;
        // double purombeta = (1-puromhypermean) / puromhyperinvconc;
        // total += puromsuffstat.GetLogProb(puromalpha,purombeta);

        double dposomalpha = 1.0 / dposomhyperinvshape;
        double dposombeta = dposomalpha / dposomhypermean;
        total += dposomsuffstat.GetLogProb(dposomalpha, dposombeta);

        // double purwalpha = purwhypermean / purwhyperinvconc;
        // double purwbeta = (1-purwhypermean) / purwhyperinvconc;
        // total += purwsuffstat.GetLogProb(purwalpha,purwbeta);

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

    // log prob for moving mixture hyper params
    double MixtureHyperLogProb() const {
        return MixtureHyperLogPrior() + MixtureHyperSuffStatLogProb();
    }

    // log prob for moving nuc rates hyper params
    // double NucRatesHyperLogProb() const {
    //    return GeneNucRatesHyperLogPrior() + NucRatesHyperSuffStatLogProb();
    //}

    // log prob for moving nuc rates
    // double NucRatesLogProb() const {
    //    return GlobalNucRatesLogPrior()  + NucRatesSuffStatLogProb();
    //}

    //-------------------
    // Moves
    // ------------------

    // general move schedule
    void MasterMove();
    void SlaveMove();

    // moving gene-specific parameters
    void GeneResampleSub(double frac);
    void MoveGeneParameters(int nrep);

    // moving global branch lengths and lambda hyperparam
    void ResampleBranchLengths();
    void MoveLambda();

    // moving hyperparams for gene-specific branch lengths
    void MoveBranchLengthsHyperParameters();
    double BranchLengthsHyperScalingMove(double tuning, int nrep);

    // moving global nuc rates
    // void MoveNucRates();

    // moving gene-specific nuc rates hyper params
    // void MoveNucRatesHyperParameters();

    // moving mixture hyper params
    void MoveMixtureHyperParameters();
    // special function for moving pi
    void ResamplePi();

    //-------------------
    // MPI send/receive
    // ------------------

    // global parameters

    void MasterSendGlobalBranchLengths();
    void SlaveReceiveGlobalBranchLengths();

    // void MasterSendGlobalNucRates();
    // void SlaveReceiveGlobalNucRates();

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

    // void SlaveSendNucPathSuffStat();
    // void MasterReceiveNucPathSuffStat();

    // hyper suffstats

    void SlaveSendBranchLengthsHyperSuffStat();
    void MasterReceiveBranchLengthsHyperSuffStat();

    // void SlaveSendNucRatesHyperSuffStat();
    // void MasterReceiveNucRatesHyperSuffStat();

    void SlaveSendMixtureHyperSuffStat();
    void MasterReceiveMixtureHyperSuffStat();

    // log likelihoods

    void SlaveSendLogProbs();
    void MasterReceiveLogProbs();

    //-------------------
    // Data structures
    // ------------------

  private:
    Tree *tree;
    CodonSequenceAlignment *refcodondata;
    const TaxonSet *taxonset;

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

    vector<double> mixhyperparam;

    // double& puromhypermean;
    // double& puromhyperinvconc;
    // IIDBeta* puromarray;
    // BetaSuffStat puromsuffstat;

    double &dposomhypermean;
    double &dposomhyperinvshape;
    IIDGamma *dposomarray;
    GammaSuffStat dposomsuffstat;

    // double& purwhypermean;
    // double& purwhyperinvconc;
    // IIDBeta* purwarray;
    // BetaSuffStat purwsuffstat;

    double &poswhypermean;
    double &poswhyperinvconc;
    IIDBernoulliBeta *poswarray;
    BernoulliBetaSuffStat poswsuffstat;

    double pihypermean;
    double pihyperinvconc;
    double &pi;

    // shared nuc rates
    // GTRSubMatrix* nucmatrix;
    // NucPathSuffStat nucpathsuffstat;

    // gene-specific nuc rates
    /*vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    IIDDirichlet* nucrelratearray;
    DirichletSuffStat nucrelratesuffstat;

    vector<double> nucstathypercenter;
    double nucstathyperinvconc;
    IIDDirichlet* nucstatarray;
    DirichletSuffStat nucstatsuffstat;
    */

    std::vector<AAMutSelM2aModel *> geneprocess;

    double lnL;
    double GeneLogPrior;

    int burnin;

    // 0: free (fixed hyper parameters)
    // 1: free and shrinkage (free hyper parameters)
    // 2: shared across genes
    int blmode;
    // int nucmode;
    // int purommode;
    int dposommode;
    // int purwmode;
    int poswmode;

    /*
    double nucrracc1, nucrracc2, nucrracc3, nucrrtot1, nucrrtot2, nucrrtot3;
    double nucstatacc1, nucstatacc2, nucstatacc3, nucstattot1, nucstattot2,
    nucstattot3;
    */
};
