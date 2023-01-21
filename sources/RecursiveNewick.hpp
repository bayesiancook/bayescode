

template<class Branch2Branch, class Branch2Node>
void RecursiveToNewick(ostream& os, const Link* from, Branch2Branch branch2branch, Branch2Node branch2node)	{
    if (from->isLeaf()) {
        os << from->GetNode()->GetName();
        os << "_";
    }
    else    {
        os << "(";
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveToNewick(os, link->Out(), branch2branch, branch2node);
            if (link->Next() != from)   {
                os << ",";
            }
        }
        os << ")";
    }
    if (! from->isRoot())	{
        os << branch2node(from->GetBranch()->GetIndex());
        os << ":";
        os << branch2branch(from->GetBranch()->GetIndex());
    }
    else    {
        os << branch2node(from->Next()->GetBranch()->GetIndex());
    }
}

template<class Branch2Branch, class Branch2Node>
void ToNewick(ostream& os, const Tree& tree, Branch2Branch branch2branch, Branch2Node branch2node)	{
    RecursiveToNewick(os, tree.GetRoot(), branch2branch, branch2node);
    os << ";\n";
}

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

void RecursiveToNewick(ostream& os, const Link* from, const vector<double>& ds, const vector<double>& dnds) {
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
        os << dnds.at(from->GetBranch()->GetIndex());
        os << ":";
        os << ds.at(from->GetBranch()->GetIndex());
    }
    else    {
        os << dnds.at(from->Next()->GetBranch()->GetIndex());
    }
}

void ToNewick(ostream& os, const Tree& tree, const vector<double>& ds, const vector<double>& dnds)    {
    RecursiveToNewick(os, tree.GetRoot(), ds, dnds);
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

void RecursiveTabulate(ostream& os, const Link* from, const BranchSelector<double>& vs, const BranchSelector<double>& vn, bool leaf)  {
    if (leaf)   {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName() << '\t' << vs.GetVal(from->GetBranch()->GetIndex()) << '\t' << vn.GetVal(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    else    {
        if (! from->isRoot())   {
            os << vs.GetTree().GetLeftMost(from) << '\t' << vs.GetTree().GetRightMost(from) << '\t' << vs.GetVal(from->GetBranch()->GetIndex()) << '\t' << vn.GetVal(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    for (const Link* link=from->Next(); link!=from; link=link->Next())  {
        RecursiveTabulate(os, link->Out(), vs, vn, leaf);
    }
}

void Tabulate(ostream& os, const BranchSelector<double>& vs, const BranchSelector<double>& vn, bool leaf) {
    RecursiveTabulate(os, vs.GetTree().GetRoot(), vs, vn, leaf);
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

void RecursiveTabulate(ostream& os, const Tree* tree, const Link* from, const vector<double>& vs, const vector<double>& vn, bool leaf)  {
    if (leaf)   {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName() << '\t' << vs.at(from->GetBranch()->GetIndex()) << '\t' << vn.at(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    else    {
        if (! from->isRoot())   {
            os << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\t' << vs.at(from->GetBranch()->GetIndex()) << '\t' << vn.at(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    for (const Link* link=from->Next(); link!=from; link=link->Next())  {
        RecursiveTabulate(os, tree, link->Out(), vs, vn, leaf);
    }
}

void Tabulate(ostream& os, const Tree* tree, const vector<double>& vs, const vector<double>& vn, bool leaf) {
    RecursiveTabulate(os, tree, tree->GetRoot(), vs, vn, leaf);
}

void RecursiveTabulate(ostream& os, const Tree* tree, const Link* from, const vector<double>& vs, const vector<double>& vn, const vector<double>& vg, bool leaf)  {
    if (leaf)   {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName() << '\t' << vs.at(from->GetBranch()->GetIndex()) << '\t' << vn.at(from->GetBranch()->GetIndex()) << '\t' << vg.at(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    else    {
        if (! from->isRoot())   {
            os << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\t' << vs.at(from->GetBranch()->GetIndex()) << '\t' << vn.at(from->GetBranch()->GetIndex()) << '\t' << vg.at(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    for (const Link* link=from->Next(); link!=from; link=link->Next())  {
        RecursiveTabulate(os, tree, link->Out(), vs, vn, vg, leaf);
    }
}

void Tabulate(ostream& os, const Tree* tree, const vector<double>& vs, const vector<double>& vn, const vector<double>& vg, bool leaf) {
    RecursiveTabulate(os, tree, tree->GetRoot(), vs, vn, vg, leaf);
}

