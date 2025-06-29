#pragma once

// includes root component as a separate item
class NucPathSuffStatBranchArray : public SimpleBranchArray<NucPathSuffStat>    {

  public:
    //! constructor (param: tree)
    NucPathSuffStatBranchArray(const Tree &intree)
        : SimpleBranchArray<NucPathSuffStat>(intree), tree(intree) {}
    ~NucPathSuffStatBranchArray() {}

    const Tree &GetTree() const /*override*/ { return tree; }

    int GetNbranch() const { return GetTree().GetNbranch(); }

    const NucPathSuffStat& GetRootVal() const {return rootsuffstat;}
    NucPathSuffStat& GetRootValRef() { return rootsuffstat; }

    //! set all suff stats to 0
    void Clear() {
        rootsuffstat.Clear();
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Clear();
        }
    }

    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    void AddSuffStat(const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                     const MGOmegaCodonSubMatrix &rootcodonsubmatrix,
                     const NodeSelector<PathSuffStat> &pathsuffstatarray) {
        RecursiveAddSuffStat(GetTree().GetRoot(), codonsubmatrixarray, rootcodonsubmatrix,
                             pathsuffstatarray);
    }

    void RecursiveAddSuffStat(const Link *from,
                              const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                              const MGOmegaCodonSubMatrix &rootcodonsubmatrix,
                              const NodeSelector<PathSuffStat> &pathsuffstatarray) {
        if (from->isRoot()) {
            rootsuffstat.AddSuffStat(rootcodonsubmatrix, 
                pathsuffstatarray.GetVal(from->GetNode()->GetIndex()));
        }
        else    {
            (*this)[from->GetBranch()->GetIndex()].AddSuffStat(
                codonsubmatrixarray.GetVal(from->GetBranch()->GetIndex()),
                pathsuffstatarray.GetVal(from->GetNode()->GetIndex()));
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveAddSuffStat(link->Out(), codonsubmatrixarray, rootcodonsubmatrix,
                                 pathsuffstatarray);
        }
    }

    void Add(const NucPathSuffStatBranchArray& from)	{
        rootsuffstat.Add(from.rootsuffstat);
        for (int i=0; i<GetNbranch(); i++) {
            (*this)[i].Add(from.GetVal(i));
        }
    }

    double GetLogProb(const BranchArray<SubMatrix> &matrixarray, const SubMatrix& rootmatrix, CodonStateSpace& cod) const {
        double total = rootsuffstat.GetLogProb(rootmatrix, cod);
        for (int i=0; i<GetNbranch(); i++) {
            total += GetVal(i).GetLogProb(matrixarray.GetVal(i), cod);
        }
        return total;
    }

    //! return array size when put into an MPI buffer
    unsigned int GetMPISize() const { return (1 + GetNbranch()) * rootsuffstat.GetMPISize(); }

    //! put array into MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        buffer << rootsuffstat;
        for (int i = 0; i < GetNbranch(); i++) {
            buffer << GetVal(i);
        }
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        buffer >> rootsuffstat;
        for (int i = 0; i < GetNbranch(); i++) {
            buffer >> (*this)[i];
        }
    }

    //! get an array from MPI buffer and then add it to this array
    void Add(const MPIBuffer &buffer) {
        rootsuffstat += buffer;
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i] += buffer;
        }
    }

  private:
    const Tree &tree;
    NucPathSuffStat rootsuffstat;
};

