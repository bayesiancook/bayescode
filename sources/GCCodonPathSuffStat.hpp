#pragma once

#include "PathSuffStat.hpp"
#include "RelativePathSuffStat.hpp"
#include "CodonSubMatrix.hpp"
#include "PoissonSuffStat.hpp"
// #include "GeneBranchArray,hpp"

// in relative time
class GCCodonPathSuffStat : public SuffStat {

    public:

    GCCodonPathSuffStat() {
        Clear();
    }
    ~GCCodonPathSuffStat() {}

    void Clear()    {
        ngc = nat = 0;
        bgc = bat = 0;
    }

    void AddSuffStat(const OmegaCodonSubMatrix &codonsubmatrix, const PathSuffStat &pathsuffstat, double branchlength)  {
        // int ncodon = codonsubmatrix.GetNstate();
        const CodonStateSpace *statespace = codonsubmatrix.GetCodonStateSpace();

        const std::map<pair<int, int>, double> &paircount = pathsuffstat.GetPairCountMap();
        // const std::map<int, double> &waitingtime = pathsuffstat.GetWaitingTimeMap();

        /*
        double tmpbgc = 0;
        double tmpbat = 0;
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            double totgcrate = 0;
            double totatrate = 0;
            int a = i->first;
            for (int b = 0; b < ncodon; b++) {
                if (b != a) {
                    if (codonsubmatrix(a, b) != 0) {
                        int pos = statespace->GetDifferingPosition(a,b);
                        int d = statespace->GetCodonPosition(pos,b);
                        if ((d == 0) || (d == 3))   {
                            totatrate += codonsubmatrix(a,b);
                        }
                        else    {
                            totgcrate += codonsubmatrix(a,b);
                        }
                    }
                }
            }
            tmpbgc += i->second * totgcrate;
            tmpbat += i->second * totatrate;
        }
        tmpbgc /= branchlength;
        tmpbat /= branchlength;
        bgc += tmpbgc;
        bat += tmpbat;
        */

        for (std::map<pair<int, int>, double>::const_iterator i = paircount.begin(); i != paircount.end(); i++) {
            int a = i->first.first;
            int b = i->first.second;
            int pos = statespace->GetDifferingPosition(a,b);
            int d = statespace->GetCodonPosition(pos,b);
            if ((d == 0) || (d == 3))   {
                nat += i->second;
            }
            else    {
                ngc += i->second;
            }
        }
    }

    void Add(const GCCodonPathSuffStat &from)   {
        ngc += from.ngc;
        nat += from.nat;
        bgc += from.bgc;
        bat += from.bat;
    }

    void Add(double gccount, double atcount, double gcbeta, double atbeta)  {
        ngc += gccount;
        nat += atcount;
        bgc += gcbeta;
        bat += atbeta;
    }

    double GetGCCount() const {
        return ngc;
    }

    double GetATCount() const {
        return nat;
    }

    double GetGCBeta() const {
        return bgc;
    }

    double GetATBeta() const {
        return bat;
    }

    GCCodonPathSuffStat &operator+=(const GCCodonPathSuffStat &from) {
        Add(from);
        return *this;
    }

    double GetGC() const    {
        return ngc / (ngc + nat);
    }

    void Normalize(double factor)   {
        ngc *= factor;
        nat *= factor;
        bgc *= factor;
        bat *= factor;
    }
        
    //! return size when put into an MPI buffer
    unsigned int GetMPISize() const { return 4; }

    //! put current value of count and beta into an MPI buffer
    void MPIPut(MPIBuffer &buffer) const { buffer << ngc << nat << bgc << bat; }

    //! get value from MPI buffer
    void MPIGet(const MPIBuffer &buffer) { buffer >> ngc >> nat >> bgc >> bat; }

    //! get a PoissonSuffStat from MPI buffer and then add it to this object
    void Add(const MPIBuffer &buffer) {
        double tmp;
        buffer >> tmp;
        ngc += tmp;
        buffer >> tmp;
        nat += tmp;

        double temp;
        buffer >> temp;
        bgc += temp;
        buffer >> temp;
        bat += temp;
    }

    void ToStream(ostream& os) const { os << ngc << '\t' << nat << '\t' << bgc << '\t' << bat; }
    void FromStream(istream& is) { is >> ngc >> nat >> bgc >> bat; }

    private:

    double ngc;
    double nat;
    double bgc;
    double bat;
};

ostream& operator<<(ostream& os, const GCCodonPathSuffStat& suffstat)    {
    suffstat.ToStream(os);
    return os;
}

istream& operator>>(istream& is, GCCodonPathSuffStat& suffstat)    {
    suffstat.FromStream(is);
    return is;
}

class GCCodonPathSuffStatBranchArray : public SimpleBranchArray<GCCodonPathSuffStat>    {

  public:

    //! constructor (param: tree)
    GCCodonPathSuffStatBranchArray(const Tree &intree)
        : SimpleBranchArray<GCCodonPathSuffStat>(intree) {
            Clear();
    }

    //! copy constructor
    GCCodonPathSuffStatBranchArray(const GCCodonPathSuffStatBranchArray& from) :
        SimpleBranchArray<GCCodonPathSuffStat>(from.GetTree()) {
            Clear();
    }

    ~GCCodonPathSuffStatBranchArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Clear();
        }
    }

    void Normalize(double factor)   {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Normalize(factor);
        }
    }

    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    void AddSuffStat(const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                     const NodeSelector<PathSuffStat> &pathsuffstatarray,
                     const BranchSelector<double>& branchlength)    {
        RecursiveAddSuffStat(GetTree().GetRoot(), codonsubmatrixarray, pathsuffstatarray, branchlength);
    }

    void RecursiveAddSuffStat(const Link *from,
                              const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                              const NodeSelector<PathSuffStat> &pathsuffstatarray,
                              const BranchSelector<double>& branchlength)   {
        if (!from->isRoot()) {
            (*this)[from->GetBranch()->GetIndex()].AddSuffStat(
                codonsubmatrixarray.GetVal(from->GetBranch()->GetIndex()),
                pathsuffstatarray.GetVal(from->GetNode()->GetIndex()),
                branchlength.GetVal(from->GetBranch()->GetIndex()));
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveAddSuffStat(link->Out(), codonsubmatrixarray, pathsuffstatarray, branchlength);
        }
    }

    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    //! site-homogeneous version
    void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix,
                     const NodeSelector<PathSuffStat> &pathsuffstatarray,
                     const BranchSelector<double>& branchlength)  {
        RecursiveAddSuffStat(GetTree().GetRoot(), codonsubmatrix, pathsuffstatarray, branchlength);
    }

    void RecursiveAddSuffStat(const Link *from,
                              const MGOmegaCodonSubMatrix& codonsubmatrix,
                              const NodeSelector<PathSuffStat> &pathsuffstatarray,
                              const BranchSelector<double>& branchlength) {
        if (!from->isRoot()) {
            (*this)[from->GetBranch()->GetIndex()].AddSuffStat(
                codonsubmatrix,
                pathsuffstatarray.GetVal(from->GetNode()->GetIndex()),
                branchlength.GetVal(from->GetBranch()->GetIndex()));
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveAddSuffStat(link->Out(), codonsubmatrix, pathsuffstatarray, branchlength);
        }
    }

    void Add(const GCCodonPathSuffStatBranchArray& from)    {
        for (int i=0; i<GetNbranch(); i++) {
            (*this)[i].Add(from.GetVal(i));
        }
    }

    void GetGC(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)   {
            into[i] = GetVal(i).GetGC();
        }
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

