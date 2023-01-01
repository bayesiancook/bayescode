#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "ConditionSpecificMeanGammaMixArray.hpp"
#include "MeanPoissonSuffStat.hpp"

class GeneBranchGammaEffects    {

    public:

    GeneBranchGammaEffects(int inNgene, int inNbranch, int indevmode,
            int infixgene_hypermean, int infixbranch_hypermean) :
        Ngene(inNgene), Nbranch(inNbranch), devmode(indevmode),
        fixgene_hypermean(infixgene_hypermean), fixbranch_hypermean(infixbranch_hypermean)  {

        branch_hypermean = 1.0;
        branch_hyperinvshape = 1.0;
        double branch_alpha = 1.0 / branch_hyperinvshape;
        double branch_beta = branch_alpha / branch_hypermean;
        branch_array = new IIDGamma(Nbranch, branch_alpha, branch_beta);
        if (fixbranch_hypermean)    {
            for (int j=0; j<Nbranch; j++) {
                (*branch_array)[j] = 1.0;
            }
        }

        gene_hypermean = 1.0;
        gene_hyperinvshape = 1.0;
        double gene_alpha = 1.0 / gene_hyperinvshape;
        double gene_beta = gene_alpha / gene_hypermean;
        gene_array = new IIDGamma(Ngene, gene_alpha, gene_beta);
        if (fixgene_hypermean)    {
            for (int i=0; i<Ngene; i++) {
                (*gene_array)[i] = 1.0;
            }
        }

        dev_invshape = 1.0;

        dev_invshape2 = 1.0;
        dev_mean2 = 2.0;
        if (devmode == 2)    {
            dev_pi = 0.01;
        }
        else    {
            dev_pi = 0;
        }

        mean_bidimarray = new BidimProduct(*gene_array, *branch_array);

        dev_bidimarray = 0;
        if (devmode)    {
            dev_bidimarray =
                new ConditionSpecificMeanGammaMixBidimArray(*mean_bidimarray,
                        dev_invshape, dev_mean2, dev_invshape2, dev_pi);
        }
    }


    double GetMeanVal(int gene, int branch) const   {
        return mean_bidimarray->GetVal(gene, branch);
    }

    double GetVal(int gene, int branch) const   {
        if (devmode)    {
           return dev_bidimarray->GetVal(gene, branch);
        }
        return mean_bidimarray->GetVal(gene, branch);
    }

    double GetMeanBranchTotal() const	{
        return branch_hypermean * Nbranch;
    }

    double GetDevInvShape() const   {
        return dev_invshape;
    }

    double GetDevPi() const {
        return dev_pi;
    }

    void TraceHeader(ostream &os, string prefix) const {
        if (! fixgene_hypermean)    {
            os << "\t" << prefix << "_genemean";
        }
        os << "\t" << prefix << "_geneinvshape";
        if (! fixbranch_hypermean)  {
            os << "\t" << prefix << "_branchmean";
        }
        os << "\t" << prefix << "_branchinvshape";

        if (devmode)    {
            os << "\t" << prefix << "_dev";
            if (devmode == 2)   {
                os << '\t' << prefix << "_pi";
                os << '\t' << prefix << "_mean2";
                // os << '\t' << prefix << "_dev2";
            }
        }
    }

    void Trace(ostream& os) const   {
        if (! fixgene_hypermean)    {
            os << '\t' << gene_hypermean;
        }
        os << '\t' << gene_hyperinvshape;
        if (! fixbranch_hypermean)  {
            os << '\t' << branch_hypermean;
        }
        os << '\t' << branch_hyperinvshape;
        if (devmode)    {
            os << '\t' << dev_invshape;
            if (devmode == 2)   {
                os << '\t' << dev_pi;
                os << '\t' << dev_mean2;
                // os << '\t' << dev_invshape2;
            }
        }
    }

    void ToStream(ostream& os) const    {
        if (! fixbranch_hypermean)  {
            os << '\t' << branch_hypermean;
        }
        os << '\t' << branch_hyperinvshape;
        os << '\t' << *branch_array;

        if (! fixgene_hypermean)    {
            os << '\t' << gene_hypermean;
        }
        os << '\t' << gene_hyperinvshape;
        os << '\t' << *gene_array;

        if (devmode)    {
            os << '\t' << dev_invshape;
            if (devmode == 2)   {
                os << '\t' << dev_pi << '\t' << dev_mean2 << '\t' << dev_invshape2;
            }
            os << '\t' << *dev_bidimarray;
        }
    }

    void FromStream(istream& is)    {
        if (! fixbranch_hypermean)  {
            is >> branch_hypermean;
        }
        is >> branch_hyperinvshape;
        is >> *branch_array;

        if (! fixgene_hypermean)    {
            is >> gene_hypermean;
        }
        is >> gene_hyperinvshape;
        is >> *gene_array;

        if (devmode)    {
            is >> dev_invshape;
            if (devmode == 2)   {
                is >> dev_pi >> dev_mean2 >> dev_invshape2;
            }
            is >> *dev_bidimarray;
        }
    }

    double GetGeneMean() const   {
        double tot = 0;
        for (int i=0; i<Ngene; i++) {
            tot += gene_array->GetVal(i);
        }
        tot /= Ngene;
        return tot;
    }

    double GetBranchMean() const   {
        double tot = 0;
        for (int j=0; j<Nbranch; j++)   {
            tot += branch_array->GetVal(j);
        }
        tot /= Nbranch;
        return tot;
    }

    void AddStats(double& mean, double& gene, double& branch, double& dev)  {
        mean += branch_hypermean * gene_hypermean;
        gene += gene_hyperinvshape;
        branch += branch_hyperinvshape;
        dev += dev_invshape;
    }

    void AddGeneArrayTo(vector<double>& array) const {
        for (int i=0; i<Ngene; i++)   {
            array[i] += gene_array->GetVal(i);
        }
    }

    void AddSquaredGeneArrayTo(vector<double>& array) const {
        for (int i=0; i<Ngene; i++)   {
            array[i] += gene_array->GetVal(i) * gene_array->GetVal(i);
        }
    }

    void AddBranchArrayTo(vector<double>& array) const {
        for (int j=0; j<Nbranch; j++)   {
            array[j] += branch_array->GetVal(j);
        }
    }

    void AddSquaredBranchArrayTo(vector<double>& array) const {
        for (int j=0; j<Nbranch; j++)   {
            array[j] += branch_array->GetVal(j) * branch_array->GetVal(j);
        }
    }

    void AddRelVarTo(vector<vector<double>>& array) const   {
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                double mean = GetMeanVal(i,j);
                double tmp = (GetVal(i,j) - mean) / mean;
                array[i][j] += tmp*tmp / dev_invshape;
            }
        }
    }

    void AddZscoreTo(vector<vector<double>>& array) const   {
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                double mean = GetMeanVal(i,j);
                double tmp = (GetVal(i,j) - mean) / mean;
                array[i][j] += tmp / sqrt(dev_invshape);
                // array[i][j] += tmp*tmp / dev_invshape;
            }
        }
    }

    void AddDevToHist(vector<double>& post, vector<double>& ppred, int offset) const {
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                double mean = GetMeanVal(i,j);
                post[offset] = (GetVal(i,j) - mean) / mean / sqrt(dev_invshape);
                double alpha = 1.0 / dev_invshape;
                double beta = alpha / mean;
                double tmp = Random::Gamma(alpha, beta);
                ppred[offset] = (tmp - mean) / mean / sqrt(dev_invshape);
                offset++;
            }
        }
    }

    void AddDevZscoreTo(vector<vector<double>>& array) const  {
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                double mean = GetMeanVal(i,j);
                double z = (GetVal(i,j) - mean) / mean / sqrt(dev_invshape);
                array[i][j] += z;
            }
        }
    }

    void Update(IIDGamma& array, double hypermean, double hyperinvshape)    {
        double alpha = 1.0 / hyperinvshape;
        double beta = alpha / hypermean;
        array.SetShape(alpha);
        array.SetScale(beta);
    }

    void Update()   {
        if (devmode)    {
            dev_bidimarray->SetParams(dev_invshape, dev_mean2, dev_invshape2, dev_pi);
        }
        Update(*branch_array, branch_hypermean, branch_hyperinvshape);
        Update(*gene_array, gene_hypermean, gene_hyperinvshape);
    }

    double GetLogPrior() const {
        double total = 0;
        total += BranchHyperLogPrior();
        total += BranchLogPrior();
        total += GeneHyperLogPrior();
        total += GeneLogPrior();
        if (devmode)    {
            total += DevHyperLogPrior();
            total += DevLogPrior();
        }
        return total;
    }

    double BranchHyperLogPrior() const { 
        double ret = 0;
        if (! fixbranch_hypermean)  {
            ret -= branch_hypermean;
        }
        ret -= branch_hyperinvshape;
        return ret;
    }

    double BranchLogPrior() const { 
        return branch_array->GetLogProb(); 
    }

    double GeneHyperLogPrior() const {
        double ret = 0;
        if (! fixgene_hypermean)    {
            ret -= gene_hypermean;
        }
        ret -= gene_hyperinvshape;
        return ret;
    }

    double GeneLogPrior() const {
        return gene_array->GetLogProb(); 
    }

    double DevHyperLogPrior() const {
        double ret = -dev_invshape;
        if (devmode == 2)   {
            ret -= dev_mean2;
            ret -= dev_invshape2;
        }
        return ret;
    }

    double DevLogPrior() const  {
        return dev_bidimarray->GetLogProb();
    }

    double GeneDevLogPrior(int gene) const  {
        return dev_bidimarray->GetRowLogProb(gene);
    }

    double BranchDevLogPrior(int branch) const  {
        return dev_bidimarray->GetColLogProb(branch);
    }


    double HyperSuffStatLogProb(const GammaSuffStat& suffstat, 
        double mean, double invshape) const    {
        double alpha = 1.0 / invshape;
        double beta = alpha / mean;
        return suffstat.GetLogProb(alpha, beta);
    }

    template<class SS>
    double SuffStatLogProb(SS get_ss) const {
        double total = 0;
        for (int i=0; i<Ngene; i++) {
            total += GeneSuffStatLogProb(i, get_ss);
        }
        return total;
    }

    template<class SS>
    double GeneSuffStatLogProb(int gene, SS get_ss) const {
        double total = 0;
        double pp = 1.0;
        for (int j=0; j<Nbranch; j++) {
            total += GeneBranchSuffStatLogProb(gene, j, get_ss, 0, pp);
        }
        return total;
    }

    template<class SS>
    double BranchSuffStatLogProb(int branch, SS get_ss) const {
        double total = 0;
        double pp = 1.0;
        for (int i=0; i<Ngene; i++) {
            total += GeneBranchSuffStatLogProb(i, branch, get_ss, 0, pp);
        }
        return total;
    }

    template<class SS>
    double GeneBranchSuffStatLogProb(int gene, int branch, SS get_ss, int gibbs_resample, double& postprob) const {

        const MeanPoissonSuffStat& ss = get_ss(gene, branch);

        postprob = 1.0;

        double mean = gene_array->GetVal(gene) * branch_array->GetVal(branch);

        if (! devmode)  {
            return ss.count * log(mean)  - ss.beta * mean;
        }

        if (devmode == 2)   {
            double alpha1 = 1.0 / dev_invshape;
            double beta1 = alpha1 / mean;
            double postalpha1 = alpha1 + ss.count;
            double postbeta1 = beta1 + ss.beta;

            double logl1 = alpha1 * log(beta1) - Random::logGamma(alpha1) 
                + Random::logGamma(postalpha1) - postalpha1 * log(postbeta1);

            double mean2 = mean * dev_mean2;
            double alpha2 = alpha1 / dev_invshape2;
            double beta2 = alpha2 / mean2;
            double postalpha2 = alpha2 + ss.count;
            double postbeta2 = beta2 + ss.beta;

            double logl2 = alpha2 * log(beta2) - Random::logGamma(alpha2) 
                + Random::logGamma(postalpha2) - postalpha2 * log(postbeta2);

            double max = (logl1 > logl2) ? logl1 : logl2;
            double l1 = exp(logl1-max);
            double l2 = exp(logl2-max);
            double p1 = (1-dev_pi) * l1;
            double p2 = dev_pi * l2;
            double tot = p1 + p2;
            p1 /= tot;
            p2 /= tot;
            postprob = p1;
            double logl = log(tot) + max;

            if (gibbs_resample) {
                if (Random::Uniform() < p2) {
                    (*dev_bidimarray)(gene,branch) = Random::Gamma(postalpha2, postbeta2);
                }
                else    {
                    (*dev_bidimarray)(gene,branch) = Random::Gamma(postalpha1, postbeta1);
                }
            }
            return  logl;
        }

        double alpha = 1.0 / dev_invshape;
        double beta = alpha / mean;

        double postalpha = alpha + ss.count;
        double postbeta = beta + ss.beta;

        double logl = alpha * log(beta) - Random::logGamma(alpha) 
            + Random::logGamma(postalpha) - postalpha * log(postbeta);

        if (gibbs_resample) {
            (*dev_bidimarray)(gene,branch) = Random::Gamma(postalpha, postbeta);
        }
        return logl;
    }

    template<class SS>
    void GibbsResample(SS get_ss) const {
        double pp = 1.0;
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                GeneBranchSuffStatLogProb(i, j, get_ss, 1, pp);
            }
        }
    }

    template<class SS>
    double GeneBranchDevLogLikelihood(int gene, int branch, SS get_ss) const   {
        double val = GetVal(gene, branch);
        const MeanPoissonSuffStat& ss = get_ss(gene, branch);
        return ss.count * log(val)  - ss.beta * val;
    }

    template<class SS>
    void AddDevPostProbsTo(vector<vector<double>>& pp, SS get_ss)   {
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                double p = 0;
                GeneBranchSuffStatLogProb(i,j, get_ss, 0, p);
                pp[i][j] += p;
            }
        }
    }

    // log prob for moving branch omega hyperparameters

    double BranchHyperLogProb() const { 
        return BranchHyperLogPrior() + 
            HyperSuffStatLogProb(branch_hypersuffstat, branch_hypermean, branch_hyperinvshape); 
    }

    // log prob for moving gene omega hyperparameters
    double GeneHyperLogProb() const { 
        return GeneHyperLogPrior() + 
            HyperSuffStatLogProb(gene_hypersuffstat, gene_hypermean, gene_hyperinvshape); 
    }

    // needs a lambda returning suffstatcount and beta for a given gene-branch pair
    template<class SS>
    double IntegratedMove(double tuning, int nrep, SS get_ss)  {

        for (int rep=0; rep<nrep; rep++)    {

            IntegratedMoveGeneArray(1.0, 1, get_ss);
            IntegratedMoveBranchArray(1.0, 1, get_ss);
            CompensatoryMove(1.0, 3);

            if (devmode)    {
                ScalingMove(dev_invshape, 0.3, 1, 
                        [this, get_ss] () {
                        return this->DevHyperLogPrior() + this->SuffStatLogProb(get_ss);},
                        [] () {});

                if (devmode == 2)   {
                    ScalingMove(dev_mean2, 0.3, 1, 1.0, 50.0, 
                        [this, get_ss] () {
                        return this->DevHyperLogPrior() + this->SuffStatLogProb(get_ss);},
                        [] () {});

                    /*
                    ScalingMove(dev_invshape2, 0.3, 1, 1.0, 10.0, 
                        [this, get_ss] () {
                        return this->DevHyperLogPrior() + this->SuffStatLogProb(get_ss);},
                        [] () {});
                    */

                    ScalingMove(dev_pi, 0.3, 1, 0, 0.3, 
                        [this, get_ss] () {
                        return this->DevHyperLogPrior() + this->SuffStatLogProb(get_ss);},
                        [] () {});
                }

                dev_bidimarray->SetParams(dev_invshape, dev_mean2, dev_invshape2, dev_pi);
            }
        }

        if (devmode)    {
            GibbsResample(get_ss);
        }

        return 1.0;
    }

    template<class SS>
    double NonIntegratedMove(double tuning, int nrep, SS get_ss)   {

        for (int rep=0; rep<nrep; rep++)    {

            if (devmode)    {
                MoveDev(1.0, 3, get_ss);
            }

            NonIntegratedMoveGeneArray(1.0, 1);
            NonIntegratedMoveBranchArray(1.0, 1);
            CompensatoryMove(1.0, 3);

            auto update = [this] () {
                this->dev_bidimarray->SetParams(dev_invshape, dev_mean2, dev_invshape2, dev_pi);};
            auto logprob = [this] () {
                return this->DevHyperLogPrior() + this->DevLogPrior();};

            if (devmode)    {
                ScalingMove(dev_invshape, 0.3, 1, logprob, update);
                if (devmode == 2)   {
                    ScalingMove(dev_mean2, 0.3, 1, 1.0, 50.0, logprob, update);
                    // ScalingMove(dev_invshape2, 0.3, 1, 1.0, 10.0, logprob, update);
                    ScalingMove(dev_pi, 0.3, 1, 0, 0.3, logprob, update);
                }
            }
        }
        return 1.0;
    }

    double MoveHyper(double tuning, int nrep)   {

        if (! fixbranch_hypermean)  {
            MoveHyperParam(branch_hypersuffstat, *branch_array,
                    branch_hypermean, branch_hypermean, branch_hyperinvshape,
                    [this] () {return this->BranchHyperLogProb();},
                    [] () {},
                    1.0,100);
        }

        MoveHyperParam(branch_hypersuffstat, *branch_array,
                branch_hyperinvshape, branch_hypermean, branch_hyperinvshape,
                [this] () {return this->BranchHyperLogProb();},
                [] () {},
                1.0,100);

        if (! fixgene_hypermean)    {
            MoveHyperParam(gene_hypersuffstat, *gene_array,
                    gene_hypermean, gene_hypermean, gene_hyperinvshape,
                    [this] () {return this->GeneHyperLogProb();},
                    [] () {},
                    1.0,100);
        }

        MoveHyperParam(gene_hypersuffstat, *gene_array,
                gene_hyperinvshape, gene_hypermean, gene_hyperinvshape,
                [this] () {return this->GeneHyperLogProb();},
                [] () {},
                1.0,100);

        return 1.0;
    }

    double MoveBranchHyperMean(double tuning, int nrep)   {
        return MoveHyperParam(branch_hypersuffstat, *branch_array,
                branch_hypermean, branch_hypermean, branch_hyperinvshape,
                [this] () {return this->BranchHyperLogProb();},
                [] () {}, tuning, nrep);
    }

    double MoveBranchHyperInvShape(double tuning, int nrep)   {
        return MoveHyperParam(branch_hypersuffstat, *branch_array,
                branch_hyperinvshape, branch_hypermean, branch_hyperinvshape,
                [this] () {return this->BranchHyperLogProb();},
                [] () {}, tuning, nrep);
    }

    double MoveGeneHyperMean(double tuning, int nrep)   {
        return MoveHyperParam(gene_hypersuffstat, *gene_array,
                gene_hypermean, gene_hypermean, gene_hyperinvshape,
                [this] () {return this->GeneHyperLogProb();},
                [] () {}, tuning, nrep);
    }

    double MoveGeneHyperInvShape(double tuning, int nrep)   {
        return MoveHyperParam(gene_hypersuffstat, *gene_array,
                gene_hyperinvshape, gene_hypermean, gene_hyperinvshape,
                [this] () {return this->GeneHyperLogProb();},
                [] () {}, tuning, nrep);
    }

    template <class SuffStat>
    double IntegratedMoveBranchArray(double tuning, int nrep, SuffStat get_suffstat) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++) {
            for (int j=0; j<Nbranch; j++) {
                double deltalogprob = -branch_array->GetLogProb(j);
                deltalogprob -= BranchSuffStatLogProb(j, get_suffstat);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*branch_array)[j] *= e;
                // update(j);
                deltalogprob += branch_array->GetLogProb(j);
                deltalogprob += BranchSuffStatLogProb(j, get_suffstat);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*branch_array)[j] /= e;
                    // update(j);
                }
            }
        }
        double ret = nacc / Nbranch / nrep;
        return ret;
    }

    template <class SuffStat>
    double IntegratedMoveGeneArray(double tuning, int nrep, SuffStat get_suffstat)    {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++) {
            for (int i=0; i<Ngene; i++) {
                double deltalogprob = -gene_array->GetLogProb(i);
                deltalogprob -= GeneSuffStatLogProb(i, get_suffstat);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*gene_array)[i] *= e;
                // update(i);
                deltalogprob += gene_array->GetLogProb(i);
                deltalogprob += GeneSuffStatLogProb(i, get_suffstat);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*gene_array)[i] /= e;
                    // update(i);
                }
            }
        }
        double ret = nacc / Ngene / nrep;
        return ret;
    }

    double NonIntegratedMoveBranchArray(double tuning, int nrep) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++) {
            for (int j=0; j<Nbranch; j++) {
                double deltalogprob = -branch_array->GetLogProb(j);
                deltalogprob -= BranchDevLogPrior(j);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*branch_array)[j] *= e;
                deltalogprob += branch_array->GetLogProb(j);
                deltalogprob += BranchDevLogPrior(j);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*branch_array)[j] /= e;
                }
            }
        }
        double ret = nacc / Nbranch / nrep;
        return ret;
    }

    double NonIntegratedMoveGeneArray(double tuning, int nrep) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++) {
            for (int i=0; i<Ngene; i++) {
                double deltalogprob = -gene_array->GetLogProb(i);
                deltalogprob -= GeneDevLogPrior(i);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*gene_array)[i] *= e;
                deltalogprob += gene_array->GetLogProb(i);
                deltalogprob += GeneDevLogPrior(i);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*gene_array)[i] /= e;
                }
            }
        }
        double ret = nacc / Ngene / nrep;
        return ret;
    }

    template<class SuffStat>
    double MoveDev(double tuning, int nrep, SuffStat get_ss) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<Ngene; i++) {
                for (int j=0; j<Nbranch; j++)   {
                    double deltalogprob = -dev_bidimarray->GetLogProb(i,j);
                    deltalogprob -= GeneBranchDevLogLikelihood(i, j, get_ss);
                    double m = tuning * (Random::Uniform() - 0.5);
                    double e = exp(m);
                    (*dev_bidimarray)(i,j) *= e;
                    deltalogprob += dev_bidimarray->GetLogProb(i,j);
                    deltalogprob += GeneBranchDevLogLikelihood(i, j, get_ss);
                    deltalogprob += m;
                    int acc = (log(Random::Uniform()) < deltalogprob);
                    if (acc) {
                        nacc++;
                    } else {
                        (*dev_bidimarray)(i,j) /= e;
                    }
                }
            }
        }
        return nacc / Ngene / Nbranch / nrep;
    }

    double CompensatoryMove(double tuning, int nrep)    {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {

            double deltalogprob = - gene_array->GetLogProb() - branch_array->GetLogProb();
            
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            for (int i=0; i<Ngene; i++) {
                (*gene_array)[i] *= e;
            }
            for (int j=0; j<Nbranch; j++)   {
                (*branch_array)[j] /= e;
            }

            deltalogprob += gene_array->GetLogProb() + branch_array->GetLogProb();
            deltalogprob += (Ngene - Nbranch) * m;

            int acc = (log(Random::Uniform()) < deltalogprob);
            if (acc) {
                nacc++;
            } else {
                for (int i=0; i<Ngene; i++) {
                    (*gene_array)[i] /= e;
                }
                for (int j=0; j<Nbranch; j++)   {
                    (*branch_array)[j] *= e;
                }
            }
        }
        double ret = nacc / nrep;
        return ret;
    }

    template <class LogProb, class UpdateF>
    double MoveHyperParam(GammaSuffStat& hypersuffstat, IIDGamma& array, 
            double& hyperparam, double& hypermean, double& hyperinvshape, 
            LogProb logprob, UpdateF update, double tuning, int nrep)    {

        hypersuffstat.Clear();
        hypersuffstat.AddSuffStat(array);

        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -logprob();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            hyperparam *= e;
            update();
            deltalogprob += logprob();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                hyperparam /= e;
                update();
            }
            ntot++;
        }

        Update(array, hypermean, hyperinvshape);

        return nacc / ntot;
    }

    template <class LogProbF, class UpdateF>
    double SlidingMove(double &x, double tuning, int nrep, double min, double max,
                       LogProbF logprobf, UpdateF updatef)  {

        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double bk = x;
            double deltalogprob = -logprobf();
            double m = tuning * (Random::Uniform() - 0.5);
            x += m;
            if (max > min) {
                while ((x < min) || (x > max)) {
                    if (x < min) {
                        x = 2 * min - x;
                    }
                    if (x > max) {
                        x = 2 * max - x;
                    }
                }
            }
            updatef();
            deltalogprob += logprobf();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x = bk;
                updatef();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    template <class LogProbF, class UpdateF>
    double ScalingMove(double &x, double tuning, int nrep, LogProbF logprobf, UpdateF updatef)  {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -logprobf();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            x *= e;
            updatef();
            deltalogprob += logprobf();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                x /= e;
                updatef();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    template <class LogProbF, class UpdateF>
    double ScalingMove(double &x, double tuning, int nrep, double min, double max, LogProbF logprobf, UpdateF updatef)  {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -logprobf();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            x *= e;
            if ((x <= min) || (x >= max)) {
                x /=e;
            }
            else    {
                updatef();
                deltalogprob += logprobf();
                deltalogprob += m;
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x /= e;
                    updatef();
                }
            }
            ntot++;
        }
        return nacc / ntot;
    }

    private:

    int Ngene;
    int Nbranch;

    int devmode;

    int fixgene_hypermean;
    int fixbranch_hypermean;

    double gene_hypermean;
    double gene_hyperinvshape;
    double branch_hypermean;
    double branch_hyperinvshape;
    double dev_invshape;
    double dev_mean2;
    double dev_invshape2;
    double dev_pi;

    IIDGamma *gene_array;
    IIDGamma *branch_array;
    BidimProduct* mean_bidimarray;
    ConditionSpecificMeanGammaMixBidimArray *dev_bidimarray;
    GammaSuffStat gene_hypersuffstat;
    GammaSuffStat branch_hypersuffstat;
};

