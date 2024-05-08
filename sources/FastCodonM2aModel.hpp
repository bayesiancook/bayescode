
#include "IIDGamma.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include "CodonSuffStat.hpp"
#include "M2aMix.hpp"
#include "MultinomialAllocationVector.hpp"
#include "ProbModel.hpp"

class CodonM2aModel : public ProbModel {

    string dsomsuffstatfile;

    // synonymous rate shared across sites
    double ds;

    // parameters of the distribution of omega across sites
    // omega0 < 1: weight purw * (1 - posw)
    // omega1 = 1: weight (1-purw) * (1-posw)
    // omega2 > 1: weight posw

    // omega0 = purom
    double purom;

    // omega2 = 1 + dposom
    double dposom;

    double posw;
    double purw;

    M2aMix *componentomegaarray;
    MultinomialAllocationVector *sitealloc;
    mutable vector<vector<double>> sitepostprobarray;
    dSOmegaPathSuffStatArray *dsomss;

    //
    // hyperparameters of the priors over the mixture parameters
    //

    // prior probability for the gene to be under positive selection (i.e. prior
    // prob that posw > 0)
    double pi;

    // Beta prior for purom (with hypermean and hyper inverse concentration)
    double puromhypermean;
    double puromhyperinvconc;

    // Gamma prior for dposom = omega_pos - 1 (with hyper mean and inverse shape parameter)
    double dposomhypermean;
    double dposomhyperinvshape;

    // Beta prior for purw
    double purwhypermean;
    double purwhyperinvconc;

    // Beta prior for posw (assuming posw>0)
    double poswhypermean;
    double poswhyperinvconc;

  public:

    int GetNsite() const    {
        return dsomss->GetSize();
        // return siteomegapathsuffstatarray->GetSize();
    }

    CodonM2aModel(string datapath, string indsomsuffstatfile, double inpi)  {

        pi = inpi;
        dsomsuffstatfile = datapath + indsomsuffstatfile;

        ds = 1.0;

        // default values
        puromhypermean = 0.5;
        puromhyperinvconc = 0.5;
        purwhypermean = 0.5;
        purwhyperinvconc = 0.5;
        poswhypermean = 0.5;
        poswhyperinvconc = 0.1;
        dposomhypermean = 1.0;
        dposomhyperinvshape = 0.5;

    }

    void Allocate() {
        purom = puromhypermean;
        dposom = dposomhypermean;
        purw = purwhypermean;
        posw = poswhypermean;

        ifstream is(dsomsuffstatfile.c_str());
        int n;
        is >> n;
        dsomss = new dSOmegaPathSuffStatArray(n);
        is >> (*dsomss);
        componentomegaarray = new M2aMix(purom, dposom + 1, purw, posw);
        sitealloc = new MultinomialAllocationVector(GetNsite(), componentomegaarray->GetWeights());
        sitepostprobarray.assign(GetNsite(), vector<double>(3, 0));
    }

    //! return value of omega_0 < 1
    double GetPurOm() const { return purom; }

    //! return value of omega_2 > 1
    double GetPosOm() const { return 1.0 + dposom; }

    //! return value of dposom  = omega_2 - 1 > 0
    double GetDPosOm() const { return dposom; }

    //! return proportion of sites under strictly purifying selection (under
    //! omega_0 < 1)
    double GetPurW() const { return purw; }

    //! return proportion of sites under positive selection
    double GetPosW() const { return posw; }

    double GetdS() const { return ds; }

    //! set omega mixture parameters to a new value
    void SetMixtureParameters(double inpurom, double indposom, double inpurw, double inposw)    {
        purom = inpurom;
        dposom = indposom;
        purw = inpurw;
        posw = inposw;
        componentomegaarray->SetParameters(purom, dposom + 1, purw, posw);
    }

    //! set omega mixture hyperparameters to a new value
    void SetMixtureHyperParameters(double inpuromhypermean, double inpuromhyperinvconc,
                                   double indposomhypermean, double indposomhyperinvshape,
                                   double inpi, double inpurwhypermean, double inpurwhyperinvconc,
                                   double inposwhypermean, double inposwhyperinvconc)   {
        puromhypermean = inpuromhypermean;
        puromhyperinvconc = inpuromhyperinvconc;
        dposomhypermean = indposomhypermean;
        dposomhyperinvshape = indposomhyperinvshape;
        pi = inpi;
        purwhypermean = inpurwhypermean;
        purwhyperinvconc = inpurwhyperinvconc;
        poswhypermean = inposwhypermean;
        poswhyperinvconc = inposwhyperinvconc;

        if (!pi) {
            poswhypermean = 0;
            poswhyperinvconc = 0;
        }

        if (!puromhyperinvconc) {
            purom = puromhypermean;
        }
        if (! dposomhyperinvshape)  {
            dposom = dposomhypermean;
        }
        if (! purwhyperinvconc) {
            purw = purwhypermean;
        }
        if (! poswhyperinvconc) {
            posw = poswhypermean;
        }
    }

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() override  {
        UpdateOmegaMixture();
    }

    void UpdateOmegaMixture()   {
        componentomegaarray->SetParameters(purom, dposom + 1, purw, posw);
    }

    void NoUpdate() {}

    void TraceHeader(ostream &os) const override    {
        os << "#logprior\tlnL";
        os << "\tds";
        os << "\tpurom\tposom\tpurw\tposw";
        os << '\n';
    }

    void Trace(ostream &os) const override  {
        os << GetLogPrior() << '\t' << GetLogLikelihood();
        os << '\t' << ds;
        os << '\t' << purom << '\t' << dposom + 1 << '\t' << purw << '\t' << posw;
        os << '\n';
    }

    //! \brief write current site post probs (of being under positive selection)
    //! on one line
    void TracePostProb(ostream &os) const   {
        for (int i = 0; i < GetNsite(); i++) {
            os << sitepostprobarray[i][2] << '\t';
        }
        os << '\n';
    }

    //! \brief write current site omega values implied by mixture model
    //! on one line
    void TraceSiteOmega(ostream &os) const  {
        for (int i = 0; i < GetNsite(); i++) {
            os << componentomegaarray->GetVal(sitealloc->GetVal(i)) << '\t';
        }
        os << '\n';
    }

    //! \brief get a copy of current site post probs (of being under positive
    //! selection) into array
    void GetSitesPostProb(double *array) const  {
        for (int i = 0; i < GetNsite(); i++) {
            array[i] = sitepostprobarray[i][2];
            if (sitepostprobarray[i][2] < 0) {
                cerr << "error in CodonM2aModel::GetSitesPostProb: negative prob\n";
                exit(1);
            }
        }
    }

    void FromStream(istream &is) override   {
        is >> ds >> purom >> dposom >> purw >> posw;
    }

    void ToStream(ostream &os) const override   {
        os << ds << '\t' << purom << '\t' << dposom << '\t' << purw << '\t' << posw << '\n';
    }

    void ToStreamHeader(ostream &os) const override {
        os << "ds" << '\t' << "purom" << '\t' << "dposom" << '\t' << "purw" << '\t' << "posw" << '\n';
    }

    //-------------------
    // Likelihood
    //-------------------

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const override {
        return GetLogPrior() + GetLogLikelihood();
        // return GetLogPrior() + GetIntegratedLogLikelihood();
    }

    //! return current value of likelihood, averaged over omega mixture
    //! allocations
    double GetLogLikelihood() const {
        int ncat = 3;

        double total = 0;
        double logp[ncat];
        const vector<double> &w = componentomegaarray->GetWeights();
        double max = 0;
        for (int i = 0; i < GetNsite(); i++) {
            for (int k = 0; k < ncat; k++) {
                logp[k] = dsomss->GetVal(i).GetLogProb(ds, componentomegaarray->GetVal(k));
                if ((!k) || (max < logp[k])) {
                    max = logp[k];
                }
            }

            double p = 0;
            for (int k = 0; k < ncat; k++) {
                double tmp = w[k] * exp(logp[k] - max);
                p += tmp;
                sitepostprobarray[i][k] = tmp;
            }
            double logl = log(p) + max;
            total += logl;
            for (int k = 0; k < ncat; k++) {
                sitepostprobarray[i][k] /= p;
            }
        }
        sitealloc->GibbsResample(sitepostprobarray);
        return total;
    }

    //! log prior over omega mixture
    double GetLogPrior() const    {
        double total = 0;
        total += dSLogPrior();
        total += OmegaLogPrior();
        return total;
    }

    double OmegaLogPrior() const    {
        double total = 0;
        total += PurOmegaLogPrior();
        total += PosOmegaLogPrior();
        total += PurWeightLogPrior();
        total += PosWeightLogPrior();
        return total;
    }

    double dSLogPrior() const   {
        return -ds/10;
    }

    //! beta prior for purom
    double PurOmegaLogPrior() const {
        if (! puromhyperinvconc)    {
            return 0;
        }
        double alpha = puromhypermean / puromhyperinvconc;
        double beta = (1 - puromhypermean) / puromhyperinvconc;
        return Random::logBetaDensity(purom, alpha, beta);
    }

    //! gamma prior for dposom
    double PosOmegaLogPrior() const {
        if (! dposomhyperinvshape)  {
            return 0;
        }
        double alpha = 1.0 / dposomhyperinvshape;
        double beta = alpha / dposomhypermean;
        return Random::logGammaDensity(dposom, alpha, beta);
    }

    //! beta prior for purw
    double PurWeightLogPrior() const    {
        if (! purwhyperinvconc) {
            return 0;
        }
        double alpha = purwhypermean / purwhyperinvconc;
        double beta = (1 - purwhypermean) / purwhyperinvconc;
        return Random::logBetaDensity(purw, alpha, beta);
    }

    //! mixture of point mass at 0 (with prob pi) and Beta distribution (with prob
    //! 1 - pi) for posw
    double PosWeightLogPrior() const    {
        if (posw) {
            if (!pi) {
                cerr << "in PosWeightLogProb: pi == 0 and posw > 0\n";
                exit(1);
            }

            if (! poswhyperinvconc)   {
                return 0;
            }
            double alpha = poswhypermean / poswhyperinvconc;
            double beta = (1 - poswhypermean) / poswhyperinvconc;
            return log(pi) + Random::logBetaDensity(posw, alpha, beta);
        } else {
            return log(1 - pi);
        }
    }

    //! Bernoulli for whether posw == 0 or > 0
    double PosSwitchLogPrior() const    {
        if (posw) {
            return log(pi);
        }
        return log(1 - pi);
    }

    //! \brief return log prob of current substitution mapping, as a function of
    //! omega mixture configuration
    //!
    //! Calculated using siteomegapathsuffstatarray
    double dSOmegaPathSuffStatLogProb() const {
        return componentomegaarray->GetPostProbArray(*dsomss, ds, sitepostprobarray);
    }

    //! \brief log prob factor to be recomputed when moving parameters of omega
    //! mixture
    double OmegaLogProb() const { return OmegaLogPrior() + dSOmegaPathSuffStatLogProb(); }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() override  {
        MoveParameters(30);
        return 1;
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep)   {
        ResampledS();
        for (int rep = 0; rep < nrep; rep++) {
            MoveOmega();
        }
    }

    void ResampledS()   {
        double count = 0;
        double beta = 0;
        for (int i=0; i<GetNsite(); i++)    {
            double om = componentomegaarray->GetVal(sitealloc->GetVal(i));
            count += dsomss->GetVal(i).GetSynCount() + dsomss->GetVal(i).GetNonSynCount();
            beta += dsomss->GetVal(i).GetSynBeta() + om * dsomss->GetVal(i).GetNonSynBeta();
        }
        ds = Random::GammaSample(1.0 + count, 10.0 + beta);
    }

    //! complete move schedule for omega mixture parameters
    void MoveOmega()    {
        if (puromhyperinvconc)  {
            SlidingMove(purom, 0.1, 10, 0, 1, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture, this);
        }
        if (purwhyperinvconc)   {
            SlidingMove(purw, 1.0, 10, 0, 1, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture, this);
        }
        if (pi != 0) {
            if (dposomhyperinvshape)    {
                ScalingMove(dposom, 1.0, 10, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture, this);
            }
            if (poswhyperinvconc)   {
                SlidingMove(posw, 1.0, 10, 0, 1, &CodonM2aModel::OmegaLogProb, &CodonM2aModel::UpdateOmegaMixture,
                        this);
            }
        }
        if ((pi != 0) && (pi != 1)) {
            SwitchPosWeight(10);
        }
        ResampleAlloc();
    }

    //! resample site allocations of omega mixture
    void ResampleAlloc()    {
        dSOmegaPathSuffStatLogProb();
        sitealloc->GibbsResample(sitepostprobarray);
    }

    //! reversible jump move on posw
    double SwitchPosWeight(int nrep)    {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double bkposw = posw;
            double deltalogprob = -PosSwitchLogPrior() - dSOmegaPathSuffStatLogProb();
            if (posw) {
                posw = 0;
            } else {
                double alpha = poswhypermean / poswhyperinvconc;
                double beta = (1 - poswhypermean) / poswhyperinvconc;
                posw = Random::BetaSample(alpha, beta);
            }
            UpdateOmegaMixture();
            deltalogprob += PosSwitchLogPrior() + dSOmegaPathSuffStatLogProb();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                posw = bkposw;
                UpdateOmegaMixture();
            }
            ntot++;
        }
        return nacc / ntot;
    }
};
