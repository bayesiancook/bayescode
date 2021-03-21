#pragma once

#include "NodeArray.hpp"
#include "BranchArray.hpp"
#include "MPIBuffer.hpp"
#include "Random.hpp"

class Chronogram : public SimpleNodeArray<double>   {

    public:

    Chronogram(const Tree &intree) : SimpleNodeArray<double>(intree)    {
        Sample();
    }

    ~Chronogram() {}

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    //! sample all entries from prior
    void Sample() {
        double age = RecursiveSample(GetRoot());
        Rescale(1.0 / age);
        // (*this)[i] = Random::GammaSample(shape, scale);
    }

    double RecursiveSample(const Link* from)  {
        double max = 0;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            double tmp = RecursiveSample(link->Out());
            if (max < tmp)  {
                max = tmp;
            }
        }
        double age = max;
        if (!from->isLeaf())    {
           age += Random::GammaSample(1.0, 1.0);
        }
        (*this)[from->GetNode()->GetIndex()] = age;
        return age;
    }

    void Rescale(double f)  {
        for (int i=0; i<GetNnode(); i++)    {
            (*this)[i] *= f;
        }
    }

    double GetDeltaTime(const Link* from) const {
        if (from->isRoot()) {
            return 0;
        }
        return GetVal(from->Out()->GetNode()->GetIndex()) - GetVal(from->GetNode()->GetIndex());
    }

    /*
    //! get total log prob summed over all branches
    double GetLogProb() {
        return 0;
    }

    //! get log prob for a given branch
    double GetLogProb(int index) { return 0; }
    */

    double LocalProposeMove(const Link* from, double tuning)  {
        if (from->isLeaf()) {
            cerr << "error in chronogram: move proposed on leaf node\n";
            exit(1);
        }
        double t = GetVal(from->GetNode()->GetIndex());
        double max = GetVal(from->Out()->GetNode()->GetIndex());
        double min = 0;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            double tmp = GetVal(link->Out()->GetNode()->GetIndex()); 
            if (min < tmp)  {
                min = tmp;
            }
        }
        t += tuning * (max-min) * (Random::Uniform() - 0.5);
        while ((t < min) || (t > max))  {
            if (t < min)    {
                t = 2*min - t;
            }
            if (t > max)    {
                t = 2*max - t;
            }
        }
        (*this)[from->GetNode()->GetIndex()] = t;
        return 0;
    }

    template<class Update, class LogProb> void MoveTimes(Update update, LogProb logprob)    {
        RecursiveMoveTimes(1.0, GetRoot(), update, logprob);
    }

    template<class Update, class LogProb> void RecursiveMoveTimes(double tuning, const Link* from, Update update, LogProb logprob)    {
        if ((! from->isRoot()) && (! from->isLeaf()))  {
            LocalMoveTime(tuning, from, update, logprob);
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveMoveTimes(tuning, link->Out(), update, logprob);
        }
        if ((! from->isRoot()) && (! from->isLeaf()))  {
            LocalMoveTime(tuning, from, update, logprob);
        }
    }

    template<class Update, class LogProb> double LocalMoveTime(double tuning, const Link* from, Update update, LogProb logprob) {
        double logprob1 = logprob(from);
        double bk = GetVal(from->GetNode()->GetIndex());
        double loghastings = LocalProposeMove(from, tuning);
        update(from);
        double logprob2 = logprob(from);

        double deltalogprob = logprob2 - logprob1 + loghastings;
        int accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted)   {
            (*this)[from->GetNode()->GetIndex()] = bk;
            update(from);
        }
        return ((double) accepted);
    }
};


class NodeAgeToBranchLengthArray : public SimpleBranchArray<double>    {

    public:

    NodeAgeToBranchLengthArray(const NodeSelector<double>& innodetree) : 
        SimpleBranchArray<double>(innodetree.GetTree()),
        nodetree(innodetree)    {
            Update();
    }

    ~NodeAgeToBranchLengthArray() {}

    const Link *GetRoot() const { return GetTree().GetRoot(); }

    void Update()   {
        RecursiveUpdate(GetRoot());
    }

    void RecursiveUpdate(const Link* from)   {
        LocalUpdate(from);
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveUpdate(link->Out());
        }
    }

    void LocalUpdate(const Link* link)  {
        if (! link->isRoot())   {
            double tmp = nodetree.GetVal(link->Out()->GetNode()->GetIndex()) - nodetree.GetVal(link->GetNode()->GetIndex());
            if (tmp <= 0)   {
                cerr << "error: negative delta age : " << tmp << '\n';
                exit(1);
            }
            (*this)[link->GetBranch()->GetIndex()] = tmp;
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
};

