#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <map>
#include <string>
#include "Random.hpp"
#include "StringStreamUtils.hpp"
using namespace std;

class TaxonSet;  // forward decl

/**
 * \brief a Node of a phylogenetic tree
 *
 * In itself, and like the Branch class,
 * the Node class does not give access to the local structure of the tree
 * (does not point to anything -- only pointed to by surrounding instances of
 * the Link class -- see Tree and Link). Branch and Node are merely storing the
 * name (string) and index (int) of the local branch or node.
 */

class Node {
  private:
    int index;
    std::string name;

  public:
    //! default constructor
    Node() : index(0), name("") {}
    //! constructor with string (giving the name of the node)
    Node(std::string s) : index(0), name(std::move(s)) {}
    //! copy constructor (pointer-style)
    Node(const Node *from) : index(from->index), name(from->name) {}

    virtual ~Node() = default;

    virtual std::string GetName() const { return name; }
    virtual void SetName(std::string inname) { name = inname; }
    int GetIndex() const { return index; }
    void SetIndex(int i) { index = i; }
};

/**
 * \brief a Branch of a phylogenetic tree
 *
 * In itself, and like the Node class,
 * the Branch class does not give access to the local structure of the tree.
 * (does not point to anything -- only pointed to by surrounding instances of
 * the Link class -- see Tree and Link). Branch and Node are merely storing the
 * name (string) and index (int) of the local branch or node
 */

class Branch {
  private:
    int index;
    std::string name;

  public:
    //! default constructor
    Branch() : index(0), name("") {}
    //! constructor with string (giving the name of the branch)
    Branch(std::string s) : index(0), name(std::move(s)) {}
    //! copy constructor (pointer style)
    Branch(const Branch *from) : index(from->index), name(from->name) {}

    virtual ~Branch() = default;

    virtual std::string GetName() const { return name; }
    virtual void SetName(std::string inname) { name = inname; }
    int GetIndex() const { return index; }
    void SetIndex(int i) { index = i; }
};

/**
 * \brief Link provides the building block for the chained structure of pointers
 * defining a tree topology invariant by re-rooting (such as stored by the Tree
 * class).
 *
 * The global structure defined by the Link pointers has the following property.
 * (1) Each branch has two links, one at each end; which access to each other
 * through the Link* out pointer. (2) The links of all branches stemming from a
 * given node form a circular chained structure, through the Link* next pointer.
 * (3) Links corresponding to tip nodes are their own next.
 * (4) There is a root link, which does not point to any branch, and is just
 * inserted in the circular chained structure at a given node of the tree. The
 * root link is its own out. Rerooting of the tree can be done by just inserting
 * this root node anywhere in the structure. From there, recursive functions
 * will make a complete traversal of the tree, from the current position of the
 * root down to all tips.
 */

class Link {
  private:
    Link *next;
    Link *out;
    Branch *branch;
    Node *node;
    int index;

  public:
    Link() {
        next = out = this;
        branch = nullptr;
        node = nullptr;
    }

    Link(const Link * /*unused*/) {
        next = out = this;
        node = nullptr;
        branch = nullptr;
    }

    Link *Next() const { return next; }
    Link *Out() const { return out; }
    Branch *GetBranch() const { return branch; }
    Node *GetNode() const { return node; }

    void SetBranch(Branch *inbranch) { branch = inbranch; }
    void SetNode(Node *innode) { node = innode; }

    void SetIndex(int i) { index = i; }

    int GetIndex() const { return index; }

    void SetOut(Link *inout) { out = inout; }

    void SetNext(Link *innext) { next = innext; }

    void AppendTo(Link *link) {
        if (link != nullptr) {
            link->next = this;
        }
    }

    void Insert(Link *link) {  // insert link after this
        link->next = next;
        next = link;
    }

    void InsertOut(Link *link) {  // insert link as out
        link->out = this;
        out = link;
    }

    bool isLeaf() const { return (next == this); }

    bool isUnary() const { return (next->Next() == this && !isLeaf()); }

    bool isRoot() const { return (out == this); }

    // degree : number of branches connecting to the node associated to this link
    int GetDegree() const {
        int d = 1;
        const Link *link = next;
        while (link != this) {
            d++;
            link = link->next;
        }
        return d;
    }

    const Link *GetUp(int &d) const {
        std::cerr << "in getup\n";
        exit(1);
        d = 1;
        const Link *link = Out();
        // const Link* link = link->Out();
        while (link->GetDegree() == 2) {
            link = link->Next()->Out();
            d++;
        }
        return link;
    }
};

/** \brief A phylogenetic tree
 *
 * The tree structure is encoded in the chained structure defined by the links
 * (see Link class), which themselves point to branches and nodes (Branch and
 * Node class).
 *
 */

class Tree {
  public:
    //! default constructor: set member pointers to 0
    Tree();

    //! copy constructor: recursively clone the entire chained structure of Link,
    //! Branch and Node objects
    Tree(const Tree *from);

    //! create a tree by reading from file (newick format expected)
    Tree(std::string filename);
    Tree(std::istream& is);

    //! recursively delete the whole chained structure of Link, Branch and Node
    //! objects
    ~Tree() /*override*/;

    //! \brief register all leaves of the tree with an external TaxonSet
    //
    //! The taxon set defines an indexing system for taxon names, with indices
    //! ranging over 0..P-1. The tree is recursively traversed, and each leaf's
    //! name is looked for in the map of the taxon set. If not found : an error is
    //! produced; otherwise, the leaf's index is set equal to the index of the
    //! corresponding taxon in TaxonSet.
    void RegisterWith(const TaxonSet *taxset);

    //! defines a global indexing system over all nodes, branches and links (such
    //! that tip node indices are in correspondance with the indexing provided by
    //! the TaxonSet).
    void SetIndices() {
        Nlink = 0;
        Nnode = GetSize();
        Nbranch = 0;
        linkmap.clear();
        nodemap.clear();
        branchmap.clear();
        SetIndices(GetRoot(), Nlink, Nnode, Nbranch);
    }

    //! const access to the root link
    Link *GetRoot() const /*override*/ { return root; }

    //! reroot: after rerooting, root->next == from
    void RootAt(Link *from);

    //! output to stream (newick format)
    void ToStream(std::ostream &os) const;
    void ToStreamWithBranchIndex(std::ostream &os) const;

    //! return total number of links
    int GetNlink() const { return Nlink; }

    //! return total number of branches
    int GetNbranch() const { return Nbranch; }

    //! return total number of nodes
    int GetNnode() const { return Nnode; }

    //! return total number of tips
    unsigned int GetSize() const { return GetSize(GetRoot()); }

    // interprets the (string) name field of the branch pointed to by given link
    // as a branch length (as a float, or double) -- not safe
    double GetBranchLength(const Link *link) const { return atof(GetBranchName(link).c_str()); }
    double GetBranchLength(int index) const { return atof(GetBranchName(index).c_str()); }

  private:
    // return const pointer to node with given index
    const Node *GetNode(int index) const {
        map<int, const Node *>::const_iterator i = nodemap.find(index);
        if (i == nodemap.end()) {
            cerr << "error in Tree::GetNode(int): not found\n";
            exit(1);
        }
        return i->second;
    }

    // return const pointer to branch with given index
    const Branch *GetBranch(int index) const {
        map<int, const Branch *>::const_iterator i = branchmap.find(index);
        if (i == branchmap.end()) {
            cerr << "error in Tree::GetBranch(int): not found\n";
            exit(1);
        }
        return i->second;
    }

    // return const pointer to link with given index
    Link *GetLink(int index) const {
        map<int, Link *>::const_iterator i = linkmap.find(index);
        if (i == linkmap.end()) {
            cerr << "error in Tree::GetLink(int): not found\n";
        }
        return i->second;
    }

    // recursively calculates the maximum height (or depth) from the given link
    // down to all of its downstream tips
    double GetMaxHeight(const Link *from) const {
        double max = 0;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            double tmp = GetMaxHeight(link->Out());
            if (max < tmp) {
                max = tmp;
            }
        }
        if (!from->isRoot()) {
            max += GetBranchLength(from);
        }
        return max;
    }

    // recursively calculates the minimum height (or depth) from the given link
    // down to all of its downstream tips
    double GetMinHeight(const Link *from) const {
        double min = -1;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            double tmp = GetMinHeight(link->Out());
            if ((min == -1) || (min > tmp)) {
                min = tmp;
            }
        }
        if (!from->isRoot()) {
            min += GetBranchLength(from);
        }
        return min;
    }

    // output to stream, renormalizing all branch lengths by given normalization
    // factor (and restoring branch lengths afterwards)
    void ToStreamRenorm(const Link *from, std::ostream &os, double normfactor) const {
        if (from->isLeaf()) {
            os << GetNodeName(from);
        } else {
            os << "(";
            for (const Link *link = from->Next(); link != from; link = link->Next()) {
                ToStreamRenorm(link->Out(), os, normfactor);
                if (link->Next() != from) {
                    os << ",";
                }
            }
            os << ")";
            os << GetNodeName(from);
        }
        if (from->isRoot()) {
            os << ";\n";
        } else {
            os << ":";
            os << GetBranchLength(from) * normfactor;
        }
    }

    std::string GetBranchName(const Link *link) const /*override*/ {
        return link->GetBranch()->GetName();
    }

    std::string GetNodeName(const Link *link) const /*override*/ {
        return link->GetNode()->GetName();
    }

    std::string GetBranchName(int index) const  {
        return GetBranch(index)->GetName();
    }

    std::string GetNodeName(int index) const    {
        return GetNode(index)->GetName();
    }


    void EraseInternalNodeName();
    void EraseInternalNodeName(Link *from);

    int GetSize(const Link *from) const {
        if (from->isLeaf()) {
            return 1;
        }
        int total = 0;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += GetSize(link->Out());
        }
        return total;

        return 0;
    }

    int GetFullSize(const Link *from) const {
        if (from->isLeaf()) {
            return 1;
        }
        int total = 1;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += GetFullSize(link->Out());
        }
        return total;

        return 0;
    }

    virtual const Link *GetLCA(std::string tax1, std::string tax2) const {
        bool found1 = false;
        bool found2 = false;
        const Link *link = RecursiveGetLCA(GetRoot(), tax1, tax2, found1, found2);
        return link;
    }

    virtual const Link *GetLCA(const Link *from1, const Link *from2) const {
        bool found1 = false;
        bool found2 = false;
        const Link *link = RecursiveGetLCA(GetRoot(), from1, from2, found1, found2);
        return link;
    }

    void Subdivide(Link *from, int Ninterpol);

    std::string Reduce(const Link *from = nullptr) {
        if (from == nullptr) {
            from = GetRoot();
        }
        if (from->isLeaf()) {
            std::cerr << from->GetNode()->GetName() << '\n';
            ;
            return from->GetNode()->GetName();
        }
        std::string name = "None";
        for (Link *link = from->Next(); link != from; link = link->Next()) {
            std::string tmp = Reduce(link->Out());
            if (tmp == "diff") {
                name = "diff";
            } else if (name == "None") {
                name = tmp;
            } else if (name != tmp) {
                name = "diff";
            }
        }
        std::cerr << '\t' << name << '\n';
        from->GetNode()->SetName(name);
        return name;

        return "";
    }

    void PrintReduced(std::ostream &os, const Link *from = nullptr) {
        if (from == nullptr) {
            from = GetRoot();
        }
        if (from->GetNode()->GetName() != "diff") {
            os << from->GetNode()->GetName();
        } else {
            os << '(';
            for (const Link *link = from->Next(); link != from; link = link->Next()) {
                PrintReduced(os, link->Out());
                if (link->Next() != from) {
                    os << ',';
                }
            }
            os << ')';
        }
        if (from->isRoot()) {
            os << ";\n";
        }
    }

    const Link *ChooseInternalNode() const {
        int n = CountInternalNodes(GetRoot());
        int m = (int)(n * Random::Uniform());
        const Link *tmp;
        const Link *chosen = ChooseInternalNode(GetRoot(), tmp, m);
        if (chosen == nullptr) {
            std::cerr << "error in choose internal node: null pointer\n";
            exit(1);
        }
        return chosen;
    }

    int CountInternalNodes(const Link *from) const;
    const Link *ChooseInternalNode(const Link *from, const Link *&fromup, int &n) const;
    int CountNodes(const Link *from) const;
    const Link *ChooseNode(const Link *from, const Link *&fromup, int &n) const;

    // delete the leaf pointed to by the next link and set everything right.
    void DeleteNextLeaf(Link *previous);

    // delete the (assumed) non-branching Node pointed to by from set everything
    // right.
    void DeleteUnaryNode(Link *from);

    // recursive function called by RegisterWith
    bool RegisterWith(const TaxonSet *taxset, Link *from, int &tot);

    map<int, const Node *> nodemap;
    map<int, const Branch *> branchmap;
    map<int, Link *> linkmap;

    // recursively set the system of node, branch and link indices
    void SetIndices(Link *from, int &linkindex, int &nodeindex, int &branchindex) {
        if (!from->isRoot()) {
            from->GetBranch()->SetIndex(branchindex);
            branchmap[branchindex] = from->GetBranch();
            branchindex++;
        }

        if (!from->isLeaf()) {
            from->GetNode()->SetIndex(nodeindex);
            nodemap[nodeindex] = from->GetNode();
            nodeindex++;
        } else {
            nodemap[from->GetNode()->GetIndex()] = from->GetNode();
        }

        if (!from->isRoot()) {
            from->Out()->SetIndex(linkindex);
            linkmap[linkindex] = from->Out();
            linkindex++;
        }
        from->SetIndex(linkindex);
        linkmap[linkindex] = from;
        linkindex++;

        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            SetIndices(link->Out(), linkindex, nodeindex, branchindex);
        }
    }

    // returns 0 if not found
    // returns link if found (then found1 and found2 must
    const Link *RecursiveGetLCA(const Link *from, std::string tax1, std::string tax2, bool &found1,
                                bool &found2) const {
        const Link *ret = nullptr;
        if (from->isLeaf()) {
            // found1 |= (from->GetNode()->GetName() == tax1);
            // found2 |= (from->GetNode()->GetName() == tax2);
            std::string name1 = GetLeafNodeName(from).substr(0, tax1.size());
            std::string name2 = GetLeafNodeName(from).substr(0, tax2.size());
            found1 |= static_cast<int>(name1 == tax1);
            found2 |= static_cast<int>(name2 == tax2);
            /*
              found1 |= (GetLeafNodeName(from) == tax1);
              found2 |= (GetLeafNodeName(from) == tax2);
            */
            if (ret == nullptr) {
                if (found1 && found2) {
                    ret = from;
                }
            }
        } else {
            for (const Link *link = from->Next(); link != from; link = link->Next()) {
                bool tmp1 = false;
                bool tmp2 = false;
                const Link *ret2 = RecursiveGetLCA(link->Out(), tax1, tax2, tmp1, tmp2);
                found1 |= static_cast<int>(tmp1);
                found2 |= static_cast<int>(tmp2);
                if (ret2 != nullptr) {
                    if (ret != nullptr) {
                        std::cerr << "error : found node twice\n";
                        std::cerr << tax1 << '\t' << tax2 << '\n';
                        ToStream(std::cerr, ret2->Out());
                        std::cerr << '\n';
                        ToStream(std::cerr, ret->Out());
                        std::cerr << '\n';
                        exit(1);
                    }
                    ret = ret2;
                }
            }
            if (ret == nullptr) {
                if (found1 && found2) {
                    ret = from;
                }
            }
        }
        return ret;
    }

    const Link *RecursiveGetLCA(const Link *from, const Link *from1, const Link *from2,
                                bool &found1, bool &found2) const {
        const Link *ret = nullptr;
        found1 |= static_cast<int>(from == from1);
        found2 |= static_cast<int>(from == from2);
        if (ret == nullptr) {
            if (found1 && found2) {
                ret = from;
            }
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            bool tmp1 = false;
            bool tmp2 = false;
            const Link *ret2 = RecursiveGetLCA(link->Out(), from1, from2, tmp1, tmp2);
            found1 |= static_cast<int>(tmp1);
            found2 |= static_cast<int>(tmp2);
            if (ret2 != nullptr) {
                if (ret != nullptr) {
                    std::cerr << "error : found node twice\n";
                    std::cerr << from1 << '\t' << from2 << '\n';
                    ToStream(std::cerr, ret2->Out());
                    std::cerr << '\n';
                    ToStream(std::cerr, ret->Out());
                    std::cerr << '\n';
                    exit(1);
                }
                ret = ret2;
            }
        }
        if (ret == nullptr) {
            if (found1 && found2) {
                ret = from;
            }
        }
        return ret;
    }

    void ReadFromStream(std::istream &is);
    // reading a tree from a stream:
    // recursively invokes the two following functions

    Link *ParseGroup(std::string input, Link *from);
    // a group is an expression of one of the two following forms:
    //
    //  (Body)Node_name
    //  (Body)Node_name:Branch_name
    //
    // where Body is a list, and Node_name and Branch_name are 2     std::strings
    // Node_name may be an empty     std::string
    //
    // Node_name cannot contain the ':' character, but Branch_name can
    // thus, if the group reads "(BODY)A:B:C"
    // then Node_name = "A" and Branch_name = "B:C"

    Link *ParseList(std::string input, Node *node);
    // a list is an expression of the form X1,X2,...Xn
    // where Xi is a group

    void RecursiveClone(const Link *from, Link *to);
    // Used by Tree(const Tree* from)
    // only clone the Links, and their mutual relations
    // does not copy the Node or Branch objects

    // deletes the link structure
    // does not delete the Node or Branch objects
    void RecursiveDelete(Link *from);

    void SetRoot(Link *link) { root = link; }

    int GetNinternalNode() const { return RecursiveGetNinternalNode(GetRoot()); }

    void ToStream(std::ostream &os, const Link *from) const;
    void ToStreamWithBranchIndex(std::ostream &os, const Link *from) const;
    double ToStreamSimplified(std::ostream &os, const Link *from) const;

  public:

    const Link *GetLeftMostLink(const Link *from) const {
        if (from->isLeaf()) {
            return from;
        }
        return GetLeftMostLink(from->Next()->Out());
    }

    const Link *GetRightMostLink(const Link *from) const {
        if (from->isLeaf()) {
            return from;
        }
        const Link *link = from->Next();
        while (link->Next() != from) {
            link = link->Next();
        }
        return GetRightMostLink(link->Out());
    }

    std::string GetLeftMost(const Link *from) const {
        if (from->isLeaf()) {
            return GetNodeName(from);
        }
        return GetLeftMost(from->Next()->Out());
    }

    std::string GetRightMost(const Link *from) const {
        if (from->isLeaf()) {
            return GetNodeName(from);
        }
        const Link *link = from->Next();
        while (link->Next() != from) {
            link = link->Next();
        }
        return GetRightMost(link->Out());
    }

  private: 
    static void Simplify() { simplify = true; }

    void PrintTab(std::ostream &os) const { RecursivePrintTab(os, GetRoot()); }

    void RecursivePrintTab(std::ostream &os, const Link *from) const {
        os << GetLeftMost(from) << '\t' << GetRightMost(from) << '\t' << GetNodeName(from) << '\n';
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursivePrintTab(os, link->Out());
        }
    }

    int RecursiveGetNinternalNode(const Link *from) const {
        int n = 0;
        if (!from->isLeaf()) {
            n++;
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            n += RecursiveGetNinternalNode(link->Out());
        }
        return n;
    }

    int RecursiveGetNnode(const Link *from) const {
        int n = 1;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            n += RecursiveGetNnode(link->Out());
        }
        return n;
    }

    string GetLeafNodeName(const Link *link) const { return GetNodeName(link); }

    static bool simplify;

    // data fields
    // just 2 pointers, to the root and to a list of taxa
    Link *root;
    int Nlink;
    int Nbranch;
    int Nnode;
};

#endif  // TREE_H
