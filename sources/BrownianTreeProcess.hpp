#pragma once

#include "NodeArray.hpp"
#include "BranchArray.hpp"
#include "Random.hpp"
#include "SuffStat.hpp"

class BrownianTreeProcess : public SimpleNodeArray<double>  {

    public:

    BrownianTreeProcess(const NodeSelector<double>& intimetree, double intau):
        SimpleNodeArray<double>(intimetree.GetTree()),
        timetree(intimetree),
        tau(intau)  {
            Sample();
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    void SetTau(double intau)   {
        tau = intau;
    }

    void Sample()   {
        RecursiveSample(GetRoot(),0);
    }

    void RecursiveSample(const Link* from, double start)  {
        double start2 = LocalForwardSample(from, start);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveSample(link->Out(), start2);
        }
    }

    double LocalForwardSample(const Link* from, double start)   {
        if (from->isRoot()) {
            (*this)[from->GetNode()->GetIndex()] = start;
            return start;
        }
        double dt = timetree.GetVal(from->Out()->GetNode()->GetIndex()) - timetree.GetVal(from->GetNode()->GetIndex());
        if (dt <= 0)    {
            cerr << "error: negative time interval\n";
            exit(1);
        }
        double y = start + sqrt(1.0/tau*dt) * Random::sNormal();
        (*this)[from->GetNode()->GetIndex()] = y;
        return y;
    }

    void Shift(double d)    {
        for (int i=0; i<GetNnode(); i++)   {
            (*this)[i] += d;
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
        // X_down ~ Normal(X_up, sigma*dt)
        if (from->isRoot()) {
            return 0;
        }
        double dt = timetree.GetVal(from->Out()->GetNode()->GetIndex()) - timetree.GetVal(from->GetNode()->GetIndex());
        double localsigma = dt/tau;
        double delta = (GetVal(from->Out()->GetNode()->GetIndex()) - GetVal(from->GetNode()->GetIndex()));
        return -0.5 * (log(localsigma) + delta*delta/localsigma);
    }

    double GetNodeLogProb(const Link* from) const   {
        double total = GetLocalLogProb(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += GetLocalLogProb(link->Out());
        }
        return total;
    }

    double GetSumOfContrasts() const    {
        return RecursiveSumOfContrasts(GetRoot());
    }

    double RecursiveSumOfContrasts(const Link* from) const {
        double total = 0;
        if (!from->isRoot())    {
            double dt = timetree.GetVal(from->Out()->GetNode()->GetIndex()) - timetree.GetVal(from->GetNode()->GetIndex());
            total += (GetVal(from->Out()->GetNode()->GetIndex()) - GetVal(from->GetNode()->GetIndex())) / sqrt(dt);
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += RecursiveSumOfContrasts(link->Out());
        }
        return total;
    }

    double LocalProposeMove(int i, double tuning)  {
        (*this)[i] += tuning * (Random::Uniform() - 0.5);
        return 0;
    }

    private:

    const NodeSelector<double>& timetree;
    double tau;
};

/*
class BranchExpoMeanArray : public SimpleBranchArray<double>    {

    public:

    BranchExpoMeanArray(const NodeSelector<double>& innodetree, double inrootshift) : 
        SimpleBranchArray<double>(innodetree.GetTree()),
        nodetree(innodetree),
        rootshift(inrootshift)  {
            Update();
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

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
        if (! from->isRoot())   {
            double tmp = 0.5 * (nodetree.GetVal(from->GetNode()->GetIndex()) + nodetree.GetVal(from->Out()->GetNode()->GetIndex()));
            (*this)[from->GetBranch()->GetIndex()] = tmp * rootshift;
        }
    }

    private:
    const NodeSelector<double>& nodetree;
    double rootshift;

};
*/

class BranchExpoLengthArray : public SimpleBranchArray<double>    {

    public:

    BranchExpoLengthArray(const NodeSelector<double>& innodetree, const NodeSelector<double>& inchrono, double inrootshift) : 
        SimpleBranchArray<double>(innodetree.GetTree()),
        nodetree(innodetree),
        chrono(inchrono),
        rootshift(inrootshift)  {
            Update();
    }

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    void SetRootRate(double in) {
        rootshift = in;
    }

    double GetTotalLength() {
        return RecursiveGetTotalLength(GetRoot());
    }

    double RecursiveGetTotalLength(const Link* from)    {
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
            double a = nodetree.GetVal(from->GetNode()->GetIndex());
            double b = nodetree.GetVal(from->Out()->GetNode()->GetIndex());
            double mean = (exp(b) - exp(a)) / (b-a);
            // double mean = 0.5 * (nodetree.GetVal(from->GetNode()->GetIndex()) + nodetree.GetVal(from->Out()->GetNode()->GetIndex()));
            double dt = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            if (dt <= 0)    {
                cerr << "error: negative time on chronogram\n";
                exit(1);
            }
            (*this)[from->GetBranch()->GetIndex()] = mean * dt * rootshift;
        }
    }

    void LocalNodeUpdate(const Link* from)  {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            LocalUpdate(link->Out());
        }
    }

    private:
    const NodeSelector<double>& nodetree;
    const NodeSelector<double>& chrono;
    double rootshift;
};

/*
class NormalSuffStat : public SuffStat  {

    public:
    NormalSuffStat() : sum2(0), n(0) {}
    ~NormalSuffStat() {}

    void Clear() {
        sum2 = 0;
        n = 0;
    }

    void AddSuffStat(double x)  {
        sum2 += x * x;
        n += c;
    }

    void Add(const NormalSuffStat &from) {
        sum2 += from.sum2;
        n += from.n;
    }

    NormalSuffStat &operator+=(const NormalSuffStat &from) {
        Add(from);
        return *this;
    }

    void AddSuffStat(const BrownianTreeProcess& from)   {
        n += from.GetTree().GetNbranch();
        sum2 += from.GetSumOfContrasts();
    }

    //! return object size, when put into an MPI buffer
    unsigned int GetMPISize() const { return 2; }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << sum2 << n; }

    //! read object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> sum2 >> n; }

    //! read a GammaSuffStat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        double temp;
        buffer >> temp;
        sum2 += temp;
        int tmp;
        buffer >> tmp;
        n += tmp;
    }

    //! return log prob, as a function of the given precision parameter
    double GetLogProb(double tau) const   {
        return 0.5 * (n * log(2*Pi*sigma) - sum2 / sigma);
    }

    double GetSum2() const { return sum2; }
    int GetN() const { return n; }

    double sum2;
    int n;
};
*/


