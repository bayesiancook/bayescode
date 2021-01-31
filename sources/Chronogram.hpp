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
        double age = max + Random::GammaSample(1.0, 1.0);
        (*this)[from->GetNode()->GetIndex()] = age;
        return age;
    }

    void Rescale(double f)  {
        for (int i=0; i<GetNnode(); i++)    {
            (*this)[i] *= f;
        }
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

