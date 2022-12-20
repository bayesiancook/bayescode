

void RecursiveToNewick(ostream& os, const Link* from, const BranchSelector<double>& ds, const BranchSelector<double>& dnds) {
    if (from->isLeaf()) {
        os << from->GetNode()->GetName();
        os << "_";
    }
    else    {
        os << "(";
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveToNewick(os, link->Out(), ds, dnds);
            if (link->Next() != from)   {
                os << ",";
            }
        }
        os << ")";
    }
    if (! from->isRoot())	{
        os << dnds.GetVal(from->GetBranch()->GetIndex());
        os << ":";
        os << ds.GetVal(from->GetBranch()->GetIndex());
    }
    else    {
        os << dnds.GetVal(from->Next()->GetBranch()->GetIndex());
    }
}

void ToNewick(ostream& os, const BranchSelector<double>& ds, const BranchSelector<double>& dnds)    {
    RecursiveToNewick(os, ds.GetTree().GetRoot(), ds, dnds);
    os << ";\n";
}

void RecursiveBranchToNewick(ostream& os, const Link* from, const BranchSelector<double>& v)    {
    if (from->isLeaf()) {
        os << from->GetNode()->GetName();
    }
    else    {
        os << "(";
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveBranchToNewick(os, link->Out(), v);
            if (link->Next() != from)   {
                os << ",";
            }
        }
        os << ")";
    }
    if (! from->isRoot())	{
        os << ":";
        os << v.GetVal(from->GetBranch()->GetIndex());
    }
}

void BranchToNewick(ostream& os, const BranchSelector<double>& v)   {
    RecursiveBranchToNewick(os, v.GetTree().GetRoot(), v);
    os << ";\n";
}

void RecursiveTabulate(ostream& os, const Link* from, const BranchSelector<double>& v, bool leaf)  {
    if (leaf)   {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName() << '\t' << v.GetVal(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    else    {
        if (! from->isRoot())   {
            os << v.GetTree().GetLeftMost(from) << '\t' << v.GetTree().GetRightMost(from) << '\t' << v.GetVal(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    for (const Link* link=from->Next(); link!=from; link=link->Next())  {
        RecursiveTabulate(os, link->Out(), v, leaf);
    }
}

void Tabulate(ostream& os, const BranchSelector<double>& v, bool leaf) {
    RecursiveTabulate(os, v.GetTree().GetRoot(), v, leaf);
}

void RecursiveTabulate(ostream& os, const Tree* tree, const Link* from, const vector<double>& v, bool leaf)  {
    if (leaf)   {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName() << '\t' << v.at(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    else    {
        if (! from->isRoot())   {
            os << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\t' << v.at(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    for (const Link* link=from->Next(); link!=from; link=link->Next())  {
        RecursiveTabulate(os, tree, link->Out(), v, leaf);
    }
}

void Tabulate(ostream& os, const Tree* tree, const vector<double>& v, bool leaf) {
    RecursiveTabulate(os, tree, tree->GetRoot(), v, leaf);
}

