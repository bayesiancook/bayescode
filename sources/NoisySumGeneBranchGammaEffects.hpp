#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "GammaBidimArray.hpp"
#include "MeanPoissonSuffStat.hpp"
#include "BidimWeightedSum.hpp"

class NoisyGeneBranchGammaEffects    {

    public:

    NoisyGeneBranchGammaEffects(int inNgene, int inNbranch, const Selector<double>* intimescale,
            int infixgene_hypermean, int infixbranch_hypermean) :
        Ngene(inNgene), Nbranch(inNbranch), timescale(intimescale),
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

        relative = 1;
        if (! timescale)    {
            relative = 0;
            timescale = branch_array;
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
        dev_mean2 = 1.0;
        refmean = 1.0;

        mean_bidimarray = new BidimProduct(*gene_array, *branch_array);
        mean2_bidimarray = new BidimHomogeneousSelector<double>(Ngene, Nbranch, refmean);
        // mean2_bidimarray = new BidimHomogeneousSelector<double>(Ngene, Nbranch, dev_mean2);

        dev1_bidimarray = new GammaBidimArray<BidimProduct>(*mean_bidimarray, dev_invshape);
        dev2_bidimarray = new GammaBidimArray<BidimHomogeneousSelector<double>>(
                *mean2_bidimarray, dev_invshape2);
        dev_bidimarray = new BidimWeightedSum(*dev1_bidimarray, *dev2_bidimarray, 
                *timescale, dev_mean2, relative);

    }

    double GetMeanVal(int gene, int branch) const   {
        return mean_bidimarray->GetVal(gene, branch);
    }

    double GetVal(int gene, int branch) const   {
        return dev_bidimarray->GetVal(gene, branch);
    }

    double GetDevInvShape() const   {
        return dev_invshape;
    }

    double GetDevMean2() const   {
        return dev_mean2;
    }

    double GetDevInvShape2() const   {
        return dev_invshape2;
    }

    const Selector<double>* GetBranchArray() const  {
        return branch_array;
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

        os << "\t" << prefix << "_dev";
        os << "\t" << prefix << "_mean2";
        os << "\t" << prefix << "_dev2";
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

        os << '\t' << dev_invshape;
        os << '\t' << dev_mean2;
        os << '\t' << dev_invshape2;
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

        os << '\t' << dev_invshape;
        os << '\t' << dev_mean2 << '\t' << dev_invshape2;
        os << '\t' << *dev1_bidimarray;
        os << '\t' << *dev2_bidimarray;
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

        is >> dev_invshape;
        is >> dev_mean2 >> dev_invshape2;
        is >> *dev1_bidimarray;
        is >> *dev2_bidimarray;
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

    void Update(IIDGamma& array, double hypermean, double hyperinvshape)    {
        double alpha = 1.0 / hyperinvshape;
        double beta = alpha / hypermean;
        array.SetShape(alpha);
        array.SetScale(beta);
    }

    void Update()   {
        dev1_bidimarray->SetInvShape(dev_invshape);
        dev2_bidimarray->SetInvShape(dev_invshape2);
        Update(*branch_array, branch_hypermean, branch_hyperinvshape);
        Update(*gene_array, gene_hypermean, gene_hyperinvshape);
    }

    double GetLogPrior() const {
        double total = 0;
        total += BranchHyperLogPrior();
        total += BranchLogPrior();
        total += GeneHyperLogPrior();
        total += GeneLogPrior();
        total += DevHyperLogPrior();
        total += DevLogPrior();
        if (std::isinf(total))  {
            cerr << "rescaled gene branch gamma: get log prior: inf\n";
            cerr << BranchLogPrior() << '\t' << GeneLogPrior() << '\t' << DevLogPrior() << '\n';
            exit(1);
        }
        if (std::isnan(total))  {
            cerr << "rescaled gene branch gamma: get log prior: nan\n";
            cerr << BranchLogPrior() << '\t' << GeneLogPrior() << '\t' << DevLogPrior() << '\n';
            exit(1);
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
        ret -= dev_mean2;
        ret -= dev_invshape2;
        return ret;
    }

    double DevLogPrior() const  {
        double ret = dev1_bidimarray->GetLogProb();
        ret += dev2_bidimarray->GetLogProb();
        if (std::isinf(ret))    {
            cerr << "dev log prior: inf\n";
            exit(1);
        }
        if (std::isnan(ret))    {
            cerr << "dev log prior: nan\n";
            exit(1);
        }
        return ret;
    }

    double GeneDevLogPrior(int gene) const  {
        return dev1_bidimarray->GetRowLogProb(gene) + dev2_bidimarray->GetRowLogProb(gene);
    }

    double BranchDevLogPrior(int branch) const  {
        return dev1_bidimarray->GetColLogProb(branch) + dev2_bidimarray->GetColLogProb(branch);
    }

    double HyperSuffStatLogProb(const GammaSuffStat& suffstat, 
            double mean, double invshape) const    {
        double alpha = 1.0 / invshape;
        double beta = alpha / mean;
        return suffstat.GetLogProb(alpha, beta);
    }

    template<class SS>
    double DevLogLikelihood(SS get_ss) const    {
        double tot = 0;
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                tot += GeneBranchDevLogLikelihood(i, j, get_ss);
            }
        }
        return tot;
    }

    template<class SS>
    double GeneDevLogLikelihood(int gene, SS get_ss) const  {
        double tot = 0;
        for (int j=0; j<Nbranch; j++)   {
            tot += GeneBranchDevLogLikelihood(gene, j, get_ss);
        }
        return tot;
    }

    template<class SS>
    double BranchDevLogLikelihood(int branch, SS get_ss) const  {
        double tot = 0;
        for (int i=0; i<Ngene; i++) {
            tot += GeneBranchDevLogLikelihood(i, branch, get_ss);
        }
        return tot;
    }

    template<class SS>
    double GeneBranchDevLogLikelihood(int gene, int branch, SS get_ss) const   {
        double val = GetVal(gene, branch);
        const MeanPoissonSuffStat& ss = get_ss(gene, branch);
        return ss.count * log(val)  - ss.beta * val;
    }

    template<class SS>
    void AddDevPostProbsTo(vector<vector<double>>& pp, SS get_ss)   {
        /*
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                double p = 0;
                GeneBranchSuffStatLogProb(i,j, get_ss, 0, p);
                pp[i][j] += p;
            }
        }
        */
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

    template<class SS, class LogProbF, class UpdateF, class GlobalLogProbF, class GlobalUpdateF>
    double NonIntegratedMove(double tuning, int nrep, SS get_ss,
            LogProbF logprob, UpdateF update, 
            GlobalLogProbF global_logprob, GlobalUpdateF global_update)   {

        for (int rep=0; rep<nrep; rep++)    {

            MoveDev1(1.0, 3, get_ss);
            MoveDev2(1.0, 3, get_ss);
            // MoveDev12(1.0, 3, get_ss);

            NonIntegratedMoveGeneArray(1.0, 1);
            NonIntegratedMoveBranchArray(1.0, 1, logprob, update);
            CompensatoryMove(1.0, 3, global_logprob, global_update);

            auto update = [this] () {
                this->Update();
            };

            auto logprob = [this] () {
                return this->DevHyperLogPrior() + this->DevLogPrior();
            };

            ScalingMove(dev_invshape, 0.3, 1, logprob, update);
            // ScalingMove(dev_mean2, 0.3, 1, logprob, update);
            ScalingMove(dev_invshape2, 0.3, 1, logprob, update);

            auto update2 = [this] () {
                this->Update();
            };

            auto logprob2 = [this, get_ss] () {
                return this->DevHyperLogPrior() + this->DevLogPrior() 
                    + this->DevLogLikelihood(get_ss);
            };

            ScalingMove(dev_mean2, 0.3, 1, logprob2, update2);

            /*
            if (fixbranch_hypermean)    {
                Mean2ArrayCompensatoryMove(*gene_array, 0.3, 1, logprob2, update2);
            }
            else    {
                Mean2ArrayCompensatoryMove(*branch_array, 0.3, 1, logprob2, update2);
            }
            */
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

    template<class LogProbF, class UpdateF>
    double NonIntegratedMoveBranchArray(double tuning, int nrep,
            LogProbF logprob, UpdateF update)   {

        double nacc = 0;
        for (int rep=0; rep<nrep; rep++) {
            for (int j=0; j<Nbranch; j++) {
                double deltalogprob = -branch_array->GetLogProb(j);
                deltalogprob -= BranchDevLogPrior(j);
                deltalogprob -= logprob(j);
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*branch_array)[j] *= e;
                update(j);
                deltalogprob += branch_array->GetLogProb(j);
                deltalogprob += BranchDevLogPrior(j);
                deltalogprob += logprob(j);
                deltalogprob += m;
                int acc = (log(Random::Uniform()) < deltalogprob);
                if (acc) {
                    nacc++;
                } else {
                    (*branch_array)[j] /= e;
                    update(j);
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
    double MoveDev1(double tuning, int nrep, SuffStat get_ss) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<Ngene; i++) {
                for (int j=0; j<Nbranch; j++)   {
                    double deltalogprob = -dev1_bidimarray->GetLogProb(i,j);
                    deltalogprob -= GeneBranchDevLogLikelihood(i, j, get_ss);
                    double m = tuning * (Random::Uniform() - 0.5);
                    double e = exp(m);
                    (*dev1_bidimarray)(i,j) *= e;
                    deltalogprob += dev1_bidimarray->GetLogProb(i,j);
                    deltalogprob += GeneBranchDevLogLikelihood(i, j, get_ss);
                    deltalogprob += m;
                    int acc = (log(Random::Uniform()) < deltalogprob);
                    if (acc) {
                        nacc++;
                    } else {
                        (*dev1_bidimarray)(i,j) /= e;
                    }
                }
            }
        }
        return nacc / Ngene / Nbranch / nrep;
    }

    template<class SuffStat>
    double MoveDev2(double tuning, int nrep, SuffStat get_ss) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<Ngene; i++) {
                for (int j=0; j<Nbranch; j++)   {
                    double deltalogprob = -dev2_bidimarray->GetLogProb(i,j);
                    deltalogprob -= GeneBranchDevLogLikelihood(i, j, get_ss);
                    double m = tuning * (Random::Uniform() - 0.5);
                    double e = exp(m);
                    (*dev2_bidimarray)(i,j) *= e;
                    deltalogprob += dev2_bidimarray->GetLogProb(i,j);
                    deltalogprob += GeneBranchDevLogLikelihood(i, j, get_ss);
                    deltalogprob += m;
                    int acc = (log(Random::Uniform()) < deltalogprob);
                    if (acc) {
                        nacc++;
                    } else {
                        (*dev2_bidimarray)(i,j) /= e;
                    }
                }
            }
        }
        return nacc / Ngene / Nbranch / nrep;
    }

    /*
    template<class SuffStat>
    double MoveDev12(double tuning, int nrep, SuffStat get_ss) {
        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {
            for (int i=0; i<Ngene; i++) {
                for (int j=0; j<Nbranch; j++)   {
                    double bk1 = (*dev1_bidimarray)(i,j);
                    double bk2 = (*dev2_bidimarray)(i,j);
                    double deltalogprob = -dev1_bidimarray->GetLogProb(i,j);
                    deltalogprob -= dev2_bidimarray->GetLogProb(i,j);
                    // deltalogprob -= GeneBranchDevLogLikelihood(i, j, get_ss);
                    double q1 = 1.0;
                    double q2 = 1.0 / timescale->GetVal(j);
                    // q1*x1 + q2*x2 = cste
                    double y1 = q1*(*dev1_bidimarray)(i,j);
                    double y2 = q2*(*dev2_bidimarray)(i,j);
                    double y = y1 + y2;
                    double x = y1 / y;
                    double m = tuning * (Random::Uniform() - 0.5);
                    x += m;
                    while ((x<0) || (x>y))  {
                        if (x<0)    {
                            x = -x;
                        }
                        if (x>y)    {
                            x = 2*y-x;
                        }
                    }
                    (*dev1_bidimarray)(i,j) = x*y/q1;
                    (*dev2_bidimarray)(i,j) = (1-x)*y/q2;
                    deltalogprob += dev1_bidimarray->GetLogProb(i,j);
                    deltalogprob += dev2_bidimarray->GetLogProb(i,j);
                    // deltalogprob += GeneBranchDevLogLikelihood(i, j, get_ss);
                    // deltalogprob += m;
                    int acc = (log(Random::Uniform()) < deltalogprob);
                    if (acc) {
                        nacc++;
                    } else {
                        (*dev1_bidimarray)(i,j) = bk1;
                        (*dev2_bidimarray)(i,j) = bk2;
                    }
                }
            }
        }
        return nacc / Ngene / Nbranch / nrep;
    }
    */

    template<class Array, class LogProbF, class UpdateF>
    double Mean2ArrayCompensatoryMove(Array& array, double tuning, int nrep,
            LogProbF logprob = [] () {return 0;},
            UpdateF update = [] () {}) {

        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {

            double deltalogprob = - array.GetLogProb();
            deltalogprob -= logprob();
            
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            for (int i=0; i<array.GetSize(); i++) {
                array[i] *= e;
            }
            double m2 = tuning * (Random::Uniform() - 0.5);
            double e2 = exp(m2);
            dev_mean2 /= e2;

            update();
            deltalogprob += array.GetLogProb();
            deltalogprob += array.GetSize() * m - m2;
            deltalogprob += logprob();

            int acc = (log(Random::Uniform()) < deltalogprob);
            if (acc) {
                nacc++;
            } else {
                for (int i=0; i<array.GetSize(); i++) {
                    array[i] /= e;
                }
                dev_mean2 *= e2;
                update();
            }
        }
        double ret = nacc / nrep;
        return ret;
    }

    template<class LogProbF, class UpdateF>
    double CompensatoryMove(double tuning, int nrep,
            LogProbF logprob = [] () {return 0;},
            UpdateF update = [] () {}) {

        double nacc = 0;
        for (int rep=0; rep<nrep; rep++)    {

            double deltalogprob = - gene_array->GetLogProb() - branch_array->GetLogProb();
            deltalogprob -= logprob();
            
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            for (int i=0; i<Ngene; i++) {
                (*gene_array)[i] *= e;
            }
            for (int j=0; j<Nbranch; j++)   {
                (*branch_array)[j] /= e;
            }

            update();
            deltalogprob += gene_array->GetLogProb() + branch_array->GetLogProb();
            deltalogprob += (Ngene - Nbranch) * m;
            deltalogprob += logprob();

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
                update();
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

    const Selector<double>* timescale;
    
    int relative;

    int fixgene_hypermean;
    int fixbranch_hypermean;

    double gene_hypermean;
    double gene_hyperinvshape;
    double branch_hypermean;
    double branch_hyperinvshape;
    double dev_invshape;
    double dev_mean2;
    double refmean;
    double dev_invshape2;

    IIDGamma *gene_array;
    IIDGamma *branch_array;
    BidimProduct* mean_bidimarray;
    BidimHomogeneousSelector<double>* mean2_bidimarray;
    GammaBidimArray<BidimProduct>* dev1_bidimarray;
    GammaBidimArray<BidimHomogeneousSelector<double>>* dev2_bidimarray;
    BidimWeightedSum* dev_bidimarray;

    GammaSuffStat gene_hypersuffstat;
    GammaSuffStat branch_hypersuffstat;
};

