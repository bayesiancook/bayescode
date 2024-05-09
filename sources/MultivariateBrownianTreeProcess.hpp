#pragma once

#include "NodeArray.hpp"
#include "BranchArray.hpp"
#include "Random.hpp"
#include "CovMatrix.hpp"
#include "ContinuousData.hpp"
#include "MultivariateNormalSuffStat.hpp"

class MultivariateBrownianTreeProcess : public SimpleNodeArray<vector<double> >  {

    public:

    MultivariateBrownianTreeProcess(const NodeSelector<double>& intimetree, const CovMatrix& insigma, const vector<double>& inrootmean, const vector<double>& inrootvar) :
        SimpleNodeArray<vector<double>>(intimetree.GetTree()),
        timetree(intimetree),
        sigma(insigma),
        rootmean(inrootmean),
        rootvar(inrootvar),
        clamp(intimetree.GetNnode(), vector<bool>(insigma.GetDim(),false))  {
            Assign(GetRoot());
            Sample();
    }

    void Assign(const Link* from)   {
        (*this)[from->GetNode()->GetIndex()].assign(GetDim(), 0);
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
            Assign(link->Out());
        }
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    int GetDim() const {
        return sigma.GetDim();
    }

    void SetAndClamp(const ContinuousData& data, int index, int fromindex)  {
        int k = 0;
        int n = 0;
        RecursiveSetAndClamp(GetRoot(), data, index, fromindex, k, n);
		cerr << data.GetCharacterName(fromindex) << " : " << n-k << " out of " << n << " missing\n";
    }

    void RecursiveSetAndClamp(const Link* from, const ContinuousData& data, int index, int fromindex, int& k, int& n)   {

		if(from->isLeaf()){
			n++;
			int tax = data.GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			if (tax != -1)	{
				double tmp = data.GetState(tax, fromindex);
				if (tmp != -1)	{
					k++;
                    (*this)[from->GetNode()->GetIndex()][index] = log(tmp);
                    clamp[from->GetNode()->GetIndex()][index] = true;
				}
			}
			else	{
				cerr << "set and clamp : " << from->GetNode()->GetName() << " not found\n";
			}
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetAndClamp(link->Out(), data, index, fromindex, k, n);
		}
	}

    void Shift(int index, double delta) {
        for (int i=0; i<GetNnode(); i++)   {
            if (! clamp[i][index])  {
                (*this)[i][index] += delta;
            }
        }
    }

    void GetContrast(const Link* from, vector<double>& contrast) const {
            double dt = timetree.GetVal(from->Out()->GetNode()->GetIndex()) - timetree.GetVal(from->GetNode()->GetIndex());
            double scaling = sqrt(dt);
            const vector<double>& up = GetVal(from->GetNode()->GetIndex());
            const vector<double>& down = GetVal(from->Out()->GetNode()->GetIndex());
            for (int i=0; i<GetDim(); i++)  {
                contrast[i] += (up[i] - down[i]) / scaling;
            }
    }

    void GetIndependentContrasts(map<const Link*, vector<double>>& contrasts,
            map<const Link*, double>& lengths) const {
        RecursiveGetIndependentContrasts(GetRoot(), contrasts, lengths);
    }

    vector<double> RecursiveGetIndependentContrasts(const Link* from, 
            map<const Link*, vector<double>>& contrasts,
            map<const Link*, double>& lengths) const    {
        if (from->isLeaf()) {
            return GetVal(from->GetNode()->GetIndex());
        }
        if (from->Next()->Next()->Next() != from)   {
            cerr << "error in get independent contrasts: tree is not binary\n";
            exit(1);
        }

        const Link* leftlink = from->Next();
        vector<double> leftval = RecursiveGetIndependentContrasts(leftlink->Out(), contrasts, lengths);
        double leftdt = timetree.GetVal(leftlink->GetNode()->GetIndex()) - timetree.GetVal(leftlink->Out()->GetNode()->GetIndex());

        const Link* rightlink = from->Next()->Next();
        vector<double> rightval = RecursiveGetIndependentContrasts(rightlink->Out(), contrasts, lengths);
        double rightdt = timetree.GetVal(rightlink->GetNode()->GetIndex()) - timetree.GetVal(rightlink->Out()->GetNode()->GetIndex());

        vector<double> nodeval(GetDim(),0);
        for (int i=0; i<GetDim(); i++)  {
            nodeval[i] = (leftval[i]/leftdt + rightval[i]/rightdt) / (1.0/leftdt + 1.0/rightdt);
        }
        double length = leftdt + rightdt;
        vector<double>& contrast = contrasts[from];
        for (int i=0; i<GetDim(); i++)  {
            contrast[i] += (leftval[i] - rightval[i]) / sqrt(leftdt + rightdt);
        }
        lengths[from] += length;
        return nodeval;
    }

    void Sample()   {
        RecursiveSample(GetRoot());
    }

    void RecursiveSample(const Link* from)  {
        LocalSample(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveSample(link->Out());
        }
    }

    void LocalSample(const Link* from) {
        if (from->isRoot()) {
            const vector<bool>& cl = clamp[from->GetNode()->GetIndex()];
            vector<double>& val = (*this)[from->GetNode()->GetIndex()];
            for (int i=0; i<GetDim(); i++)  {
                if (! cl[i])    {
                    val[i] = sqrt(rootvar[i]) * Random::sNormal() + rootmean[i];
                }
            }
        }
        else    {
            double dt = timetree.GetVal(from->Out()->GetNode()->GetIndex()) - timetree.GetVal(from->GetNode()->GetIndex());
            if (dt <= 0)    {
                cerr << "error: negative time interval\n";
                exit(1);
            }
            double scaling = sqrt(dt);

            const vector<double>& initval = (*this)[from->Out()->GetNode()->GetIndex()];
            vector<double>& finalval = (*this)[from->GetNode()->GetIndex()];
            const vector<bool>& cl = clamp[from->GetNode()->GetIndex()];

            // draw multivariate normal from sigma
            vector<double> contrast(GetDim(), 0);
            sigma.MultivariateNormalSample(contrast);

            // not conditional on clamped entries
            for (int i=0; i<GetDim(); i++)  {
                if (! cl[i])    {
                    finalval[i] = initval[i] + scaling*contrast[i];
                }
            }
        }
    }

    double GetLogProb() const {
        return RecursiveGetLogProb(GetRoot());
    }

    double RecursiveGetLogProb(const Link* from) const  {
        double total = GetLocalLogProb(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += RecursiveGetLogProb(link->Out());
        }
        return total;
    }

    double GetLocalLogProb(const Link* from) const  {

        // Normal(rootmean, rootvar)
        if (from->isRoot()) {
            const vector<double>& val = GetVal(from->GetNode()->GetIndex());
            double total = 0;
            for (int i=0; i<GetDim(); i++)  {
                double delta = val[i] - rootmean[i];
                total -= 0.5 * (log(2*Pi*rootvar[i]) + delta*delta/rootvar[i]);
            }
            return total;
        }

        // X_down ~ Normal(X_up, sigma*dt)
        // X = (X_down - X_up)
        // Y = (X_down - X_up)/sqrt(dt)
        // P(Y)dY = p(X)dX
        // p(X) = p(Y) dY/dX = p(Y) / sqrt(dt)^GetDim()
        // log P(X) = log P(Y) - 0.5 * GetDim() * log(dt)

        double dt = timetree.GetVal(from->Out()->GetNode()->GetIndex()) - timetree.GetVal(from->GetNode()->GetIndex());
        double scaling = sqrt(dt);

        const vector<double>& up = GetVal(from->GetNode()->GetIndex());
        const vector<double>& down = GetVal(from->Out()->GetNode()->GetIndex());

        vector<double> contrast(GetDim(), 0);
        for (int i=0; i<GetDim(); i++)  {
            contrast[i] = (up[i] - down[i])/scaling;
        }
        return sigma.logMultivariateNormalDensity(contrast) - 0.5*GetDim()*log(dt);
    }

    double GetNodeLogProb(const Link* from) const   {
        double total = GetLocalLogProb(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += GetLocalLogProb(link->Out());
        }
        return total;
    }

    void GetSampleCovarianceMatrix(CovMatrix& covmat, int& n) const    {
        RecursiveGetSampleCovarianceMatrix(GetRoot(), covmat, n);
    }

    void RecursiveGetSampleCovarianceMatrix(const Link* from, CovMatrix& covmat, int& n) const  {

        if (! from->isRoot())   {
            vector<double> contrast(GetDim(), 0);
            GetContrast(from, contrast);
            for (int i=0; i<GetDim(); i++)  {
                for (int j=0; j<GetDim(); j++)  {
                    covmat.add(i, j, contrast[i]*contrast[j]);
                }
            }
            n++;
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveGetSampleCovarianceMatrix(link->Out(), covmat, n);
        }
    }

    void AddSuffStat(MultivariateNormalSuffStat& to) const {
        GetSampleCovarianceMatrix(to.covmat, to.n);
    }

    void GetSumOfContrasts(vector<double>& sum) const    {
        return RecursiveSumOfContrasts(GetRoot(), sum);
    }

    void RecursiveSumOfContrasts(const Link* from, vector<double>& sum) const {
        if (!from->isRoot())    {
            vector<double> contrast(GetDim(), 0);
            GetContrast(from, contrast);
            for (int i=0; i<GetDim(); i++)  {
                sum[i] += contrast[i];
            }
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveSumOfContrasts(link->Out(), sum);
        }
    }

    template<class Update, class LogProb> void SingleNodeMove(int index, double tuning, Update update, LogProb logprob)   {
        RecursiveSingleNodeMove(index, tuning, GetRoot(), update, logprob);
    }

    template<class Update, class LogProb> void RecursiveSingleNodeMove(int index, double tuning, const Link* from, Update update, LogProb logprob)    {

        if (! clamp[from->GetNode()->GetIndex()][index])    {
            LocalSingleNodeMove(index, tuning, from, update, logprob);
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveSingleNodeMove(index, tuning, link->Out(), update, logprob);
        }
        if (! clamp[from->GetNode()->GetIndex()][index])    {
            LocalSingleNodeMove(index, tuning, from, update, logprob);
        }
    }

    template<class Update, class LogProb> double LocalSingleNodeMove(int index, double tuning, const Link* from, Update update, LogProb logprob) {
        double logprob1 = logprob(from);
        double delta = tuning * (Random::Uniform() - 0.5);
        (*this)[from->GetNode()->GetIndex()][index] += delta;
        update(from);
        double logprob2 = logprob(from);

        double deltalogprob = logprob2 - logprob1;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted)   {
            (*this)[from->GetNode()->GetIndex()][index] -= delta;
            update(from);
        }
        return ((double) accepted);
    }

    template<class Update, class LogProb>
    double FilterMove(int index, int nspan, double min_delta, double max_delta, Update update, LogProb logprob)   {
        std::vector<std::vector<double>> proposal(GetTree().GetNnode(), std::vector<double>(2*nspan+1, 0));
        std::vector<std::vector<double>> logcondl(GetTree().GetNnode(), std::vector<double>(2*nspan+1, 0));
        double delta = min_delta + (max_delta - min_delta)*Random::Uniform();
        for (int i=0; i<GetTree().GetNnode(); i++)  {
            double center = (*this)[i][index];
            for (int k=0; k<2*nspan+1; k++)   {
                proposal[i][k] = center + delta*(k - nspan)/nspan;
            }
        }
        BackwardFilterMove(GetRoot(), index, nspan, update, logprob, proposal, logcondl);
        double nacc = 0;
        double ntry = 0;
        FilterSampleRoot(index, nspan, update, logprob, proposal, logcondl, nacc, ntry);
        ForwardFilterMove(GetRoot(), index, nspan, update, logprob, proposal, logcondl, nacc, ntry);
        return nacc/ntry;
    }

    template<class Update, class LogProb>
    void BackwardFilterMove(const Link* from, int index, int nspan, Update update, LogProb logprob, std::vector<std::vector<double>>& proposal, std::vector<std::vector<double>>& logcondl)  {

        int from_index = from->GetNode()->GetIndex();
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            BackwardFilterMove(link->Out(), index, nspan, update, logprob, proposal, logcondl);
            std::vector<double> tmp = BackwardPropagate(link->Out(), index, nspan, update, logprob, proposal, logcondl);
            for (size_t k=0; k<logcondl[from_index].size(); k++)    {
                logcondl[from_index][k] += tmp[k];
            }
        }

        if (from->isRoot()) {
            for (size_t k=0; k<logcondl[from_index].size(); k++)    {
                (*this)[from_index][index] = proposal[from_index][k];
                logcondl[from_index][k] += GetLocalLogProb(from);
            }
            (*this)[from_index][index] = proposal[from_index][nspan];
        }
    }

    template<class Update, class LogProb>
    void ForwardFilterMove(const Link* from, int index, int nspan, Update update, LogProb logprob, std::vector<std::vector<double>>& proposal, std::vector<std::vector<double>>& logcondl, double& nacc, double& ntry)  {

        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            if (clamp[link->Out()->GetNode()->GetIndex()][index])    {
                if (! link->Out()->isLeaf())   {
                    std::cerr << "error in forward filter: only tips can be clamped\n";
                    exit(1);
                }
                update(link->Out());
            }
            else    {
                ForwardPropagate(link->Out(), index, nspan, update, logprob, proposal, logcondl, nacc, ntry);
                ForwardFilterMove(link->Out(), index, nspan, update, logprob, proposal, logcondl, nacc, ntry);
            }
        }
    }

    template<class Update, class LogProb>
    std::vector<double> ClampedBackwardPropagate(const Link* from, int index, int nspan, Update update, LogProb logprob, std::vector<std::vector<double>>& proposal, std::vector<std::vector<double>>& logcondl)  {

        int from_index = from->GetNode()->GetIndex();
        int parent_index = from->Out()->GetNode()->GetIndex();

        std::vector<double> tmp(logcondl[from_index].size(), 0);
        for (size_t k=0; k<logcondl[from_index].size(); k++)    {
            (*this)[parent_index][index] = proposal[from_index][k];
            update(from);
            tmp[k] = GetLocalLogProb(from) + logprob(from);
            if (std::isnan(tmp[k])) {
                cerr << "in ClampedBackward: return value is nan\n";
                exit(1);
            }
        }
        (*this)[parent_index][index] = proposal[parent_index][nspan];
        update(from);
        return tmp;
    }


    template<class Update, class LogProb>
    std::vector<double> UnclampedBackwardPropagate(const Link* from, int index, int nspan, Update update, LogProb logprob, std::vector<std::vector<double>>& proposal, std::vector<std::vector<double>>& logcondl)  {

        // up[k] = sum_l q_kl * down[l]

        int from_index = from->GetNode()->GetIndex();
        int parent_index = from->Out()->GetNode()->GetIndex();

        std::vector<std::vector<double>> logl(logcondl[from_index].size(), std::vector<double>(logcondl[from_index].size(), 0));

        for (size_t k=0; k<logcondl[from_index].size(); k++)    {
            (*this)[parent_index][index] = proposal[parent_index][k];
            for (size_t l=0; l<logcondl[from_index].size(); l++)    {
                (*this)[from_index][index] = proposal[from_index][l];
                update(from);
                logl[k][l] = GetLocalLogProb(from) + logprob(from) + logcondl[from_index][l];
                if (std::isnan(logl[k][l]))  {
                    cerr << "logl is nan\n";
                    cerr << GetLocalLogProb(from) << '\t' << logprob(from) << '\t' << logcondl[from_index][l] << '\n';
                    exit(1);
                }
            }
        }
        // double max = logl[0][0];
        std::vector<double> logtmp(logcondl[from_index].size(), 0);
        std::vector<double> tmp(logcondl[from_index].size(), 0);
        for (size_t k=0; k<logcondl[from_index].size(); k++)    {
            double max = logl[k][0];
            for (size_t l=0; l<logcondl[from_index].size(); l++)    {
                if (max < logl[k][l])   {
                    max = logl[k][l];
                }
            }
            for (size_t l=0; l<logcondl[from_index].size(); l++)    {
                tmp[k] += exp(logl[k][l] - max);
            }
            logtmp[k] = log(tmp[k]) + max;
        }

        double max = logtmp[0];
        for (size_t k=0; k<logcondl[from_index].size(); k++)    {
            if (max < logtmp[k])    {
                max = logtmp[k];
            }
        }

        if (std::isnan(max))    {
            cerr << "max is nan\n";
            exit(1);
        }
        if (std::isinf(max))    {
            cerr << "max is inf\n";
            exit(1);
        }

        (*this)[from_index][index] = proposal[from_index][nspan];
        (*this)[parent_index][index] = proposal[parent_index][nspan];
        update(from);

        return logtmp;
    }

    template<class Update, class LogProb>
    std::vector<double> BackwardPropagate(const Link* from, int index, int nspan, Update update, LogProb logprob, std::vector<std::vector<double>>& proposal, std::vector<std::vector<double>>& logcondl)  {
        if (clamp[from->GetNode()->GetIndex()][index]) {
            return ClampedBackwardPropagate(from, index, nspan, update, logprob, proposal, logcondl);
        }
        return UnclampedBackwardPropagate(from, index, nspan, update, logprob, proposal, logcondl);
    }

    template<class Update, class LogProb>
    void ForwardPropagate(const Link* from, int index, int nspan, Update update, LogProb logprob, std::vector<std::vector<double>>& proposal, std::vector<std::vector<double>>& logcondl, double& nacc, double& ntry)  {

        int from_index = from->GetNode()->GetIndex();
        std::vector<double> logtmp(logcondl[from_index].size(), 0);
        for (size_t l=0; l<logcondl[from_index].size(); l++) {
            (*this)[from_index][index] = proposal[from_index][l];
            update(from);
            logtmp[l] = GetLocalLogProb(from) + logprob(from) + logcondl[from_index][l];
        }
        double max = logtmp[0];
        for (size_t l=0; l<logcondl[from_index].size(); l++) {
            if (max < logtmp[l])    {
                max = logtmp[l];
            }
        }
        double tot = 0;
        std::vector<double> tmp(logcondl[from_index].size(), 0);
        for (size_t l=0; l<logcondl[from_index].size(); l++) {
            tmp[l] = exp(logtmp[l] - max);
            tot += tmp[l];
        }
        for (size_t l=0; l<logcondl[from_index].size(); l++) {
            tmp[l] /= tot;
        }
        int choose = Random::DrawFromDiscreteDistribution(tmp);
        if ((choose < 0) || (choose >= int(logcondl[from_index].size()))) {
            std::cerr << "error: choose out of range\n";
            exit(1);
        }
        (*this)[from_index][index] = proposal[from_index][choose];
        if (choose != nspan)    {
            nacc++;
        }
        ntry++;
        update(from);
    }

    template<class Update, class LogProb>
    void FilterSampleRoot(int index, int nspan, Update update, LogProb logprob, std::vector<std::vector<double>>& proposal, std::vector<std::vector<double>>& logcondl, double& nacc, double& ntry)  {

        int root_index = GetRoot()->GetIndex();

        double max = logcondl[root_index][0];
        for (size_t k=0; k<logcondl[root_index].size(); k++)    {
            if (max < logcondl[root_index][k])    {
                max = logcondl[root_index][k];
            }
        }

        std::vector<double> tmp(logcondl[root_index].size(), 0);
        double tot = 0;
        for (size_t k=0; k<logcondl[root_index].size(); k++)    {
            tmp[k] = exp(logcondl[root_index][k] - max);
            tot += tmp[k];
        }
        for (size_t k=0; k<logcondl[root_index].size(); k++)    {
            tmp[k] /= tot;
        }

        int choose = Random::DrawFromDiscreteDistribution(tmp);
        if ((choose < 0) || (choose >= int(logcondl[root_index].size()))) {
            std::cerr << "error: choose out of range\n";
            exit(1);
        }
        if (choose != nspan)    {
            nacc++;
        }
        ntry++;
        (*this)[root_index][index] = proposal[root_index][choose];
    }

    private:

    const NodeSelector<double>& timetree;
    const CovMatrix& sigma;
    const vector<double>& rootmean;
    const vector<double>& rootvar;
    vector<vector<bool>> clamp;
};

class MVBranchExpoLengthArray : public SimpleBranchArray<double>    {

    public:

    MVBranchExpoLengthArray(const NodeSelector<vector<double>>& innodetree, const NodeSelector<double>& inchrono, int inidx) :
        SimpleBranchArray<double>(innodetree.GetTree()),
        nodetree(innodetree),
        chrono(inchrono),
        idx(inidx)  {
            Update();
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    double GetTotalLength() const {
        return RecursiveGetTotalLength(GetRoot());
    }

    double RecursiveGetTotalLength(const Link* from) const {
        double tot = 0;
        if (! from->isRoot())   {
            tot += GetVal(from->GetBranch()->GetIndex());
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            tot += RecursiveGetTotalLength(link->Out());
        }
        return tot;
    }

    void Update()   {
        RecursiveUpdate(GetRoot());
    }

    void RecursiveUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveUpdate(link->Out());
        }
    }

    void LocalUpdate(const Link* from)  {
        if (!from->isRoot()) {
            double up = nodetree.GetVal(from->GetNode()->GetIndex())[idx];
            double down = nodetree.GetVal(from->Out()->GetNode()->GetIndex())[idx];
            double mean = 0.5 * (exp(up) + exp(down));
            // double mean = (exp(up) - exp(down)) / (up - down);
            double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            if (dt <= 0)    {
                cerr << "error: negative time on chronogram\n";
                exit(1);
            }
            (*this)[from->GetBranch()->GetIndex()] = mean * dt;
        }
    }

    void LocalNodeUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            LocalUpdate(link->Out());
        }
    }

    private:
    const NodeSelector<vector<double>>& nodetree;
    const NodeSelector<double>& chrono;
    int idx;
};

class MVBranchExpoMeanArray : public SimpleBranchArray<double>    {

    public:

    MVBranchExpoMeanArray(const NodeSelector<vector<double>>& innodetree, int inidx) :
        SimpleBranchArray<double>(innodetree.GetTree()),
        nodetree(innodetree),
        idx(inidx)  {
            Update();
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    double GetMean() const  {
        return GetTotal() / GetTree().GetNbranch();
    }

    double GetTotal() const {
        return RecursiveGetTotal(GetRoot());
    }

    double RecursiveGetTotal(const Link* from) const {
        double tot = 0;
        if (! from->isRoot())   {
            tot += GetVal(from->GetBranch()->GetIndex());
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            tot += RecursiveGetTotal(link->Out());
        }
        return tot;
    }

    void Update()   {
        RecursiveUpdate(GetRoot());
    }

    void RecursiveUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveUpdate(link->Out());
        }
    }

    void LocalUpdate(const Link* from)  {
        if (!from->isRoot()) {
            double up = nodetree.GetVal(from->GetNode()->GetIndex())[idx];
            double down = nodetree.GetVal(from->Out()->GetNode()->GetIndex())[idx];
            double mean = 0.5 * (exp(up) + exp(down));
            // double mean = (exp(up) - exp(down)) / (up - down);
            (*this)[from->GetBranch()->GetIndex()] = mean;
        }
    }

    void LocalNodeUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            LocalUpdate(link->Out());
        }
    }

    private:
    const NodeSelector<vector<double>>& nodetree;
    int idx;
};


