#pragma once

#include "PathSuffStat.hpp"
#include "RelativePathSuffStat.hpp"
#include "CodonSubMatrix.hpp"
#include "PoissonSuffStat.hpp"
// #include "GeneBranchArray,hpp"

// in relative time
class dSOmegaPathSuffStat : public SuffStat {

    public:

    dSOmegaPathSuffStat() {
        Clear();
    }
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

    void AddSuffStat(const OmegaCodonSubMatrix &codonsubmatrix, const RelativePathSuffStat &pathsuffstat, double omega) {
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
        tmpbnonsyn /= omega;
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

    double GetLogProbdSIntegrated(double l, double omega, double dt, double nu) const   {
        //double alpha = dt / nu;
        double alpha = 1.0 / nu;
        double alphapost = alpha + nsyn + nnonsyn;
        double betapost = alpha + l*(bsyn + bnonsyn*omega);
        return alpha*log(alpha) - Random::logGamma(alpha) - alphapost*log(betapost) + Random::logGamma(alphapost) + (nsyn + nnonsyn)*log(l) + nnonsyn*log(omega);
    }

    double GetLogProbOmIntegrated(double l, double omega, double dt, double nu) const   {
        // double alpha = dt / nu;
        double alpha = 1.0 / nu;
        double alphapost = alpha + nnonsyn;
        double betapost = alpha + l*bnonsyn*omega;
        return alpha*log(alpha) - Random::logGamma(alpha) - alphapost*log(betapost) + Random::logGamma(alphapost) + (nsyn + nnonsyn)*log(l) + nnonsyn*log(omega) - l*bsyn;
    }

    void Add(const dSOmegaPathSuffStat &from, double omega = 1.0) {
        nsyn += from.nsyn;
        nnonsyn += from.nnonsyn;
        bsyn += from.bsyn;
        bnonsyn += omega*from.bnonsyn;
    }

    void Add(double syncount, double synbeta, double nonsyncount, double nonsynbeta)    {
        nsyn += syncount;
        nnonsyn += nonsyncount;
        bsyn += synbeta;
        bnonsyn += nonsynbeta;
    }

    void TodSSuffStat(PoissonSuffStat& suffstat, double omega) const  {
        suffstat.AddSuffStat(nsyn + nnonsyn, bsyn + omega*bnonsyn);
    }

    void ToOmSuffStat(PoissonSuffStat& suffstat, double l) const  {
        suffstat.AddSuffStat(nnonsyn, l*bnonsyn);
    }

    void AddWNdSSuffStat(PoissonSuffStat& suffstat, double l, double omega) const   {
        suffstat.AddSuffStat(nsyn + nnonsyn, l*(bsyn + bnonsyn*omega));
    }

    void AddWNOmSuffStat(PoissonSuffStat& suffstat, double l, double omega) const   {
        suffstat.AddSuffStat(nnonsyn, l*omega*bnonsyn);
    }

    double GetCount() const {
        return nsyn + nnonsyn;
    }

    double GetSynCount() const {
        return nsyn;
    }

    double GetNonSynCount() const  {
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

    double GetdNdS() const    {
        if ((!bsyn) || (!bnonsyn) || (!nsyn))    {
            return 0;
        }
        return (nnonsyn / bnonsyn) / (nsyn / bsyn);
    }

    double GetdS() const    {
        if (! bsyn) {
            return 0;
        }
        return nsyn / bsyn;
    }

    double GetdN() const    {
        if (! bnonsyn)  {
            return 0;
        }
        return nnonsyn / bnonsyn;
    }

    void Normalize(double factor)   {
        nsyn *= factor;
        nnonsyn *= factor;
        bsyn *= factor;
        bnonsyn *= factor;
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

    void ToStream(ostream& os) const { os << nsyn << '\t' << nnonsyn << '\t' << bsyn << '\t' << bnonsyn; }
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
        : SimpleBranchArray<dSOmegaPathSuffStat>(intree) {
            Clear();
    }

    //! copy constructor
    dSOmegaPathSuffStatBranchArray(const dSOmegaPathSuffStatBranchArray& from) :
        SimpleBranchArray<dSOmegaPathSuffStat>(from.GetTree()) {
            Clear();
    }

    ~dSOmegaPathSuffStatBranchArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetNbranch(); i++) {
            (*this)[i].Clear();
        }
    }

    void TodSSuffStat(BranchArray<PoissonSuffStat>& suffstat, const BranchSelector<double>& omega) const {
        for (int i = 0; i < GetNbranch(); i++) {
            GetVal(i).TodSSuffStat(suffstat[i], omega.GetVal(i));
        }
    }

    void ToOmSuffStat(BranchArray<PoissonSuffStat>& suffstat, const BranchSelector<double>& l) const {
        for (int i = 0; i < GetNbranch(); i++) {
            GetVal(i).ToOmSuffStat(suffstat[i], l.GetVal(i));
        }
    }

    void AddWNdSSuffStat(BranchArray<PoissonSuffStat>& suffstat, const BranchSelector<double>& length, const BranchSelector<double>& omega, const BranchSelector<double>& wnom) const {
        for (int i = 0; i < GetNbranch(); i++) {
            GetVal(i).AddWNdSSuffStat(suffstat[i], length.GetVal(i), omega.GetVal(i)*wnom.GetVal(i));
        }
    }

    void AddWNOmSuffStat(BranchArray<PoissonSuffStat>& suffstat, const BranchSelector<double>& length, const BranchSelector<double>& omega, const BranchSelector<double>& wnds) const {
        for (int i = 0; i < GetNbranch(); i++) {
            GetVal(i).AddWNOmSuffStat(suffstat[i], length.GetVal(i)*wnds.GetVal(i), omega.GetVal(i));
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

    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    //! site-homogeneous version
    void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix,
                     const NodeSelector<PathSuffStat> &pathsuffstatarray,
                     const BranchSelector<double>& branchlength, double omega)  {
        RecursiveAddSuffStat(GetTree().GetRoot(), codonsubmatrix, pathsuffstatarray, branchlength, omega);
    }

    void RecursiveAddSuffStat(const Link *from,
                              const MGOmegaCodonSubMatrix& codonsubmatrix,
                              const NodeSelector<PathSuffStat> &pathsuffstatarray,
                              const BranchSelector<double>& branchlength, double omega) {
        if (!from->isRoot()) {
            (*this)[from->GetBranch()->GetIndex()].AddSuffStat(
                codonsubmatrix,
                pathsuffstatarray.GetVal(from->GetNode()->GetIndex()),
                branchlength.GetVal(from->GetBranch()->GetIndex()),
                omega);
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveAddSuffStat(link->Out(), codonsubmatrix, pathsuffstatarray, branchlength, omega);
        }
    }

    // RELATIVE version
    //! compute omega suff stats and do a member-wise addition -- for Muse and
    //! Gaut codon matrices
    void AddSuffStat(const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                     const NodeSelector<RelativePathSuffStat> &pathsuffstatarray,
                     const BranchSelector<double>& branchomega) {
        RecursiveAddSuffStat(GetTree().GetRoot(), codonsubmatrixarray, pathsuffstatarray, branchomega);
    }

    void RecursiveAddSuffStat(const Link *from,
                              const BranchSelector<MGOmegaCodonSubMatrix> &codonsubmatrixarray,
                              const NodeSelector<RelativePathSuffStat> &pathsuffstatarray,
                              const BranchSelector<double>& branchomega) {
        if (!from->isRoot()) {
            (*this)[from->GetBranch()->GetIndex()].AddSuffStat(
                codonsubmatrixarray.GetVal(from->GetBranch()->GetIndex()),
                pathsuffstatarray.GetVal(from->GetNode()->GetIndex()),
                branchomega.GetVal(from->GetBranch()->GetIndex()));
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveAddSuffStat(link->Out(), codonsubmatrixarray, pathsuffstatarray, branchomega);
        }
    }

    void Add(const dSOmegaPathSuffStatBranchArray& from, double omega = 1)    {
        for (int i=0; i<GetNbranch(); i++) {
            (*this)[i].Add(from.GetVal(i), omega);
        }
    }

    void Add(const Tree& syncount, const Tree& synbeta, const Tree& nonsyncount, const Tree& nonsynbeta)    {
        for (int i=0; i<GetNbranch(); i++) {
            (*this)[i].Add(syncount.GetBranchLength(i), synbeta.GetBranchLength(i), nonsyncount.GetBranchLength(i), nonsynbeta.GetBranchLength(i));
        }
    }

    void GetdNdS(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)   {
            into[i] = GetVal(i).GetdNdS();
        }
    }

    double GetTotaldS() const {
        double M = 0;
        double L = 0;
        for (int i=0; i<GetNbranch(); i++)   {
            M += GetVal(i).GetSynCount();
            L += GetVal(i).GetSynBeta();
        }
        return M/L;
    }

    double GetTotaldN() const {
        double M = 0;
        double L = 0;
        for (int i=0; i<GetNbranch(); i++)   {
            M += GetVal(i).GetNonSynCount();
            L += GetVal(i).GetNonSynBeta();
        }
        return M/L;
    }

    double GetMeandNdS() const  {
        return GetTotaldN() / GetTotaldS();
    }

    void GetdS(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)   {
            into[i] = GetVal(i).GetdS();
        }
    }

    void GetdN(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)   {
            into[i] = GetVal(i).GetdN();
        }
    }

    void GetSynCount(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)  {
            into[i] = GetVal(i).GetSynCount();
        }
    }

    void GetSynBeta(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)  {
            into[i] = GetVal(i).GetSynBeta();
        }
    }

    void GetNonSynCount(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)  {
            into[i] = GetVal(i).GetNonSynCount();
        }
    }

    void GetNonSynBeta(BranchArray<double>& into) const    {
        for (int i=0; i<GetNbranch(); i++)  {
            into[i] = GetVal(i).GetNonSynBeta();
        }
    }

    double GetMaxZS(const BranchSelector<double>& dS) const {
        double max = 0;
        for (int i=0; i<GetNbranch(); i++)  {
            double obs = GetVal(i).GetSynCount();
            double exp = GetVal(i).GetSynBeta()*dS.GetVal(i);
            double y = (obs - exp) / (sqrt(exp) + 1);
            double z = fabs(y);
            if (max < z)    {
                max = z;
            }
        }
        return max;
    }

    double GetMaxZN(const BranchSelector<double>& dS, const BranchSelector<double>& beta, double omega) const {
        double max = 0;
        for (int i=0; i<GetNbranch(); i++)  {
            double obs = GetVal(i).GetNonSynCount();
            double exp = GetVal(i).GetNonSynBeta()*dS.GetVal(i)*beta.GetVal(i)*omega;
            double y = (obs - exp) / (sqrt(exp) + 1);
            double z = fabs(y);
            if (max < z)    {
                max = z;
            }
        }
        return max;
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

