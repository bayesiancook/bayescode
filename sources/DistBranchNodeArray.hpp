#pragma once

#include "BranchArray.hpp"
#include "NodeArray.hpp"
#include <list>

class DistBranchNodeArray   {

    public:

    DistBranchNodeArray(const Tree& intree) : tree(intree), branchlength(intree.GetNbranch(), list<double>()), branchdist(intree.GetNbranch(), list<double>()), nodedist(intree.GetNnode(), list<double>()), counter(0)  {}

    virtual const Tree &GetTree() const { return tree; }
    const Link* GetRoot() const { return tree.GetRoot(); }
    int GetNbranch() const { return tree.GetNbranch(); }
    int GetNnode() const { return tree.GetNnode(); }

    void AddFromChrono(const NodeSelector<double>& chrono, const NodeSelector<vector<double>>& nodeval, int index)    {
        RecursiveAddFromChrono(GetRoot(), chrono, nodeval, index);
        counter++;
    }

    void RecursiveAddFromChrono(const Link* from, const NodeSelector<double>& chrono, const NodeSelector<vector<double>>& nodeval, int index)    {
        if (! from->isRoot())   {
            double bl = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            branchlength[from->GetBranch()->GetIndex()].push_front(bl);
        }
        nodedist[from->GetNode()->GetIndex()].push_front(nodeval.GetVal(from->GetNode()->GetIndex())[index]);
        if (! from->isRoot())   {
            double up = nodeval.GetVal(from->GetNode()->GetIndex())[index];
            double down = nodeval.GetVal(from->Out()->GetNode()->GetIndex())[index];
            double val = log(0.5 * (exp(up) + exp(down)));
            branchdist[from->GetBranch()->GetIndex()].push_front(val);
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveAddFromChrono(link->Out(), chrono, nodeval, index);
        }
    }

    void Sort() {
        for (int i=0; i<GetNbranch(); i++) {
            branchlength[i].sort();
            branchdist[i].sort();
        }
        for (int i=0; i<GetNnode(); i++)    {
            nodedist[i].sort();
        }
    }

    double GetBranchLengthMean(int index) const  {
		const list<double>& l = branchlength[index];
        double mean = 0;
		for (auto val : l) {
            mean += val; 
		}
        mean /= l.size();
        return mean;
    }

    double GetBranchLengthQuantile(int index, double q) const {
		const list<double>& l = branchlength[index];
		auto i = l.begin();
		int n = ((int) (((double) l.size())*q));
		for (int j=0; j<n; j++)	{
			i++;
		}
        return *i;
    }

    double GetBranchMean(int index) const  {
		const list<double>& l = branchdist[index];
        double mean = 0;
		for (auto val : l) {
            mean += val; 
		}
        mean /= l.size();
        return mean;
    }

    double GetBranchQuantile(int index, double q) const {
		const list<double>& l = branchdist[index];
		auto i = l.begin();
		int n = ((int) (((double) l.size())*q));
		for (int j=0; j<n; j++)	{
			i++;
		}
        return *i;
    }

    double GetNodeMean(int index) const  {
		const list<double>& l = nodedist[index];
        double mean = 0;
		for (auto val : l) {
            mean += val; 
		}
        mean /= l.size();
        return mean;
    }

    double GetNodeQuantile(int index, double q) const {
		const list<double>& l = nodedist[index];
		auto i = l.begin();
		int n = ((int) (((double) l.size())*q));
		for (int j=0; j<n; j++)	{
			i++;
		}
        return *i;
    }

    void MedianToStream(ostream& os, bool with_node_name = true) const {
        RecursiveMedianToStream(os, GetRoot(), with_node_name);
        os << ";\n";
    }

    void RecursiveMedianToStream(ostream& os, const Link* from, bool with_node_name) const {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName();
            if (with_node_name) {
                os << "_";
            }
        }
        else    {
            os << "(";
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                RecursiveMedianToStream(os,link->Out(), with_node_name);
                if (link->Next() != from)   {
                    os << ",";
                }
            }
            os << ")";
        }
        if (with_node_name) {
            os << exp(GetNodeQuantile(from->GetNode()->GetIndex(), 0.5));
        }
        if (! from->isRoot())    {
            os << ":";
            os << GetBranchLengthMean(from->GetBranch()->GetIndex());
        }
    }

    void TabulateNodeMedianToStream(ostream& os) const  {
        RecursiveTabulateNodeMedianToStream(os, GetRoot());
    }

    void RecursiveTabulateNodeMedianToStream(ostream& os, const Link* from) const {
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveTabulateNodeMedianToStream(os,link->Out());
        }
        os << tree.GetLeftMost(from) << '\t' << tree.GetRightMost(from);
        if (! from->isRoot())   {
            os << '\t' << GetBranchLengthMean(from->GetBranch()->GetIndex());
        }
        else    {
            os << '\t' << 0;
        }
        os << '\t' << exp(GetNodeQuantile(from->GetNode()->GetIndex(), 0.5));
        os << '\n';
    }

    void TabulateBranchMedianToStream(ostream& os) const    {
        RecursiveTabulateBranchMedianToStream(os, GetRoot());
    }

    void RecursiveTabulateBranchMedianToStream(ostream& os, const Link* from) const {
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveTabulateBranchMedianToStream(os,link->Out());
        }
        if (! from->isRoot())   {
            os << tree.GetLeftMost(from) << '\t' << tree.GetRightMost(from);
            os << '\t' << GetBranchLengthMean(from->GetBranch()->GetIndex());
            os << '\t' << exp(GetBranchQuantile(from->GetBranch()->GetIndex(), 0.5));
            os << '\n';
        }
    }

    private:

    const Tree& tree;
    vector<list<double>> branchlength;
    vector<list<double>> branchdist;
    vector<list<double>> nodedist;
    int counter;
};

