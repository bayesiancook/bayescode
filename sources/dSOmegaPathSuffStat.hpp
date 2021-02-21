#pragma once

#include "PathSuffStat.hpp"

// in relative time
class dSOmegaPathSuffStat : public SuffStat {

    public:

    dSOmegaPathSuffStat() {}
    ~dSOmegaPathSuffStat() {}

    void Clear()    {
        nsyn = nnonsyn = 0;
        bsyn = bnonsyn = 0;
    }

    void AddSuffStat(const OmegaCodonSubMatrix &codonsubmatrix, const PathSuffStat &pathsuffstat, double branchlength, double omega) {
        int ncodon = codonsubmatrix.GetNstate();
        const CodonStateSpace *statespace = codonsubmatrix.GetCodonStateSpace();

        const std::map<pair<int, int>, double> &paircount = pathsuffstat.GetPairCountMap();
        const std::map<int, double> &waitingtime = pathsuffstat.GetWaitingTimeMap();

        double tmpbsyn = 0;
        double tmpbnonsyn = 0;
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            double totsynrate = 0;
            double totnonsynrate = 0;
            int a = i->first;
            for (int b = 0; b < ncodon; b++) {
                if (b != a) {
                    if (codonsubmatrix(a, b) != 0) {
                        if (!statespace->Synonymous(a, b)) {
                            totnonsynrate += codonsubmatrix(a, b);
                        }
                        else    {
                            totsynrate += codonsubmatrix(a, b);
                        }
                    }
                }
            }
            tmpbsyn += i->second * totsynrate;
            tmpbnonsyn += i->second * totnonsynrate;
        }
        tmpbsyn /= branchlength;
        tmpbnonsyn /= branchlength*omega;
        bsyn += tmpbsyn;
        bnonsyn += tmpbnonsyn;

        for (std::map<pair<int, int>, double>::const_iterator i = paircount.begin(); i != paircount.end(); i++) {
            if (!statespace->Synonymous(i->first.first, i->first.second)) {
                nnonsyn += i->second;
            }
            else    {
                nsyn += i->second;
            }
        }
    }

    double GetLogProb(double l, double omega) const { 
        return (nsyn + nnonsyn)*log(l) + nnonsyn*log(omega) - l*(bsyn + bnonsyn*omega);
    }

    void Add(const dSOmegaPathSuffStat &from) {
        nsyn += from.nsyn;
        nnonsyn += from.nnonsyn;
        bsyn += from.bsyn;
        bnonsyn += from.bnonsyn;
    }

    int GetCount() const {
        return nsyn + nnonsyn;
    }

    int GetSynCount() const {
        return nsyn;
    }

    int GetNonSynCount() const  {
        return nnonsyn;
    }

    double GetSynBeta() const  {
        return bsyn;
    }

    double GetNonSynBeta() const    {
        return bnonsyn;
    }

    double GetBeta(double omega) const  {
        return bsyn + omega*bnonsyn;
    }

    dSOmegaPathSuffStat &operator+=(const dSOmegaPathSuffStat &from) {
        Add(from);
        return *this;
    }

    //! return size when put into an MPI buffer
    unsigned int GetMPISize() const { return 4; }

    //! put current value of count and beta into an MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << nsyn << nnonsyn << bsyn << bnonsyn; }

    //! get value from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> nsyn >> nnonsyn >> bsyn >> bnonsyn; }

    //! get a PoissonSuffStat from MPI buffer and then add it to this object
    void Add(const MPIBuffer &buffer) {
        double tmp;
        // int tmp;
        buffer >> tmp;
        nsyn += tmp;
        buffer >> tmp;
        nnonsyn += tmp;

        double temp;
        buffer >> temp;
        bsyn += temp;
        buffer >> temp;
        bnonsyn += temp;
    }

    void ToStream(ostream& os) const { os << nsyn << nnonsyn << bsyn << bnonsyn; }
    void FromStream(istream& is) { is >> nsyn >> nnonsyn >> bsyn >> bnonsyn; }

    private:

    // int nsyn;
    // int nnonsyn;
    double nsyn;
    double nnonsyn;
    double bsyn;
    double bnonsyn;
};

ostream& operator<<(ostream& os, const dSOmegaPathSuffStat& suffstat)    {
    suffstat.ToStream(os);
    return os;
}

istream& operator>>(istream& is, dSOmegaPathSuffStat& suffstat)    {
    suffstat.FromStream(is);
    return is;
}

class dSOmegaPathSuffStatBranchArray : public SimpleBranchArray<dSOmegaPathSuffStat>    {

  public:
    //! constructor (param: tree)
    dSOmegaPathSuffStatBranchArray(const Tree &intree)
        : SimpleBranchArray<dSOmegaPathSuffStat>(intree) {}
    ~dSOmegaPathSuffStatBranchArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Clear();
        }
    }

    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    void AddSuffStat(const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                     const NodeSelector<PathSuffStat> &pathsuffstatarray,
                     const BranchSelector<double>& branchlength, const BranchSelector<double>& branchomega) {
        RecursiveAddSuffStat(GetTree().GetRoot(), codonsubmatrixarray, pathsuffstatarray, branchlength, branchomega);
    }

    void RecursiveAddSuffStat(const Link *from,
                              const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                              const NodeSelector<PathSuffStat> &pathsuffstatarray,
                              const BranchSelector<double>& branchlength, const BranchSelector<double>& branchomega) {
        if (!from->isRoot()) {
            (*this)[from->GetBranch()->GetIndex()].AddSuffStat(
                codonsubmatrixarray.GetVal(from->GetBranch()->GetIndex()),
                pathsuffstatarray.GetVal(from->GetNode()->GetIndex()),
                branchlength.GetVal(from->GetBranch()->GetIndex()),
                branchomega.GetVal(from->GetBranch()->GetIndex()));
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveAddSuffStat(link->Out(), codonsubmatrixarray, pathsuffstatarray, branchlength, branchomega);
        }
    }

    void Add(const dSOmegaPathSuffStatBranchArray& from)    {
        for (int i=0; i<GetNbranch(); i++) {
            (*this)[i].Add(from.GetVal(i));
        }
    }

    //! return total log prob over array, given an array of omega_i's of same size
    double GetLogProb(const BranchSelector<double>& branchlength, const BranchSelector<double>& branchomega) const  {
        double total = 0;
        for (int i=0; i<GetNbranch(); i++) {
            total += GetVal(i).GetLogProb(branchlength.GetVal(i), branchomega.GetVal(i));
        }
        return total;
    }

    //! return array size when put into an MPI buffer
    unsigned int GetMPISize() const { return GetVal(0).GetMPISize() * GetNbranch(); }

    //! put array into MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        for (int i=0; i<GetNbranch(); i++) {
            buffer << GetVal(i);
        }
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        for (int i=0; i<GetNbranch(); i++) {
            buffer >> (*this)[i];
        }
    }

    //! get an array from MPI buffer and then add it to this array
    void Add(const MPIBuffer &buffer) {
        for (int i=0; i<GetNbranch(); i++) {
            (*this)[i] += buffer;
        }
    }

    void ToStream(ostream& os) const {
        for (int i=0; i<GetNbranch(); i++) {
            os << GetVal(i) << '\t';
        }
        os << '\n';
    }

    void FromStream(istream& is)    {
        for (int i=0; i<GetNbranch(); i++) {
            is >> (*this)[i];
        }
    }



};
