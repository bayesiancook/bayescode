
#include "PathSuffStat.hpp"
#include "CodonSuffStat.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include "CodonSubMatrix.hpp"
#include "StrandSymmetricIrreversibleSubMatrix.hpp"
#include "GeneBranchGammaEffects.hpp"

class MeanNucPathSuffStat : public SuffStat {
  public:
    MeanNucPathSuffStat()
        : rootcount(4, 0), paircount(4, vector<double>(4, 0)), pairbeta(4, vector<double>(4, 0)) {}

    ~MeanNucPathSuffStat() {}

    //! set suff stat to 0
    void Clear() {
        for (int i = 0; i < Nnuc; i++) {
            rootcount[i] = 0;
            for (int j = 0; j < Nnuc; j++) {
                paircount[i][j] = 0;
                pairbeta[i][j] = 0;
            }
        }
    }

    double GetRootCount(int i) const    {
        return rootcount[i];
    }

    double GetPairCount(int i, int j) const {
        return paircount[i][j];
    }

    double GetPairBeta(int i, int j) const  {
        return pairbeta[i][j];
    }

    //! Note that the resulting 4x4 nuc path suff stat depends on other aspects of
    //! the codon matrix (e.g. the value of omega)
    void AddSuffStat(const NucCodonSubMatrix &codonmatrix, const PathSuffStat &codonpathsuffstat) {
        const CodonStateSpace *cod = codonmatrix.GetCodonStateSpace();
        const SubMatrix *nucmatrix = codonmatrix.GetNucMatrix();

        // root part
        const std::map<int, double> &codonrootcount = codonpathsuffstat.GetRootCountMap();
        for (std::map<int, double>::const_iterator i = codonrootcount.begin();
             i != codonrootcount.end(); i++) {
            int codon = i->first;
            rootcount[cod->GetCodonPosition(0, codon)] += i->second;
            rootcount[cod->GetCodonPosition(1, codon)] += i->second;
            rootcount[cod->GetCodonPosition(2, codon)] += i->second;
        }

        const std::map<int, double> &waitingtime = codonpathsuffstat.GetWaitingTimeMap();
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            int codon = i->first;
            for (int c2 = 0; c2 < cod->GetNstate(); c2++) {
                if (c2 != codon) {
                    int pos = cod->GetDifferingPosition(codon, c2);
                    if (pos < 3) {
                        int n1 = cod->GetCodonPosition(pos, codon);
                        int n2 = cod->GetCodonPosition(pos, c2);
                        pairbeta[n1][n2] +=
                            i->second * codonmatrix(codon, c2) / (*nucmatrix)(n1, n2);
                    }
                }
            }
        }

        const std::map<pair<int, int>, double> &codonpaircount = codonpathsuffstat.GetPairCountMap();
        for (std::map<pair<int, int>, double>::const_iterator i = codonpaircount.begin();
             i != codonpaircount.end(); i++) {
            int cod1 = i->first.first;
            int cod2 = i->first.second;
            int pos = cod->GetDifferingPosition(cod1, cod2);
            if (pos == 3) {
                cerr << "error in codon conj path suffstat\n";
                exit(1);
            }
            int n1 = cod->GetCodonPosition(pos, cod1);
            int n2 = cod->GetCodonPosition(pos, cod2);
            paircount[n1][n2] += i->second;
        }
    }

    void AddSuffStat(const NucCodonSubMatrix &codonmatrix, const RelativePathSuffStat &codonpathsuffstat, double length) {
        const CodonStateSpace *cod = codonmatrix.GetCodonStateSpace();
        const SubMatrix *nucmatrix = codonmatrix.GetNucMatrix();

        // root part
        const std::map<int, double> &codonrootcount = codonpathsuffstat.GetRootCountMap();
        for (std::map<int, double>::const_iterator i = codonrootcount.begin();
             i != codonrootcount.end(); i++) {
            int codon = i->first;
            rootcount[cod->GetCodonPosition(0, codon)] += i->second;
            rootcount[cod->GetCodonPosition(1, codon)] += i->second;
            rootcount[cod->GetCodonPosition(2, codon)] += i->second;
        }

        const std::map<int, double> &waitingtime = codonpathsuffstat.GetWaitingTimeMap();
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            int codon = i->first;
            for (int c2 = 0; c2 < cod->GetNstate(); c2++) {
                if (c2 != codon) {
                    int pos = cod->GetDifferingPosition(codon, c2);
                    if (pos < 3) {
                        int n1 = cod->GetCodonPosition(pos, codon);
                        int n2 = cod->GetCodonPosition(pos, c2);
                        pairbeta[n1][n2] +=
                            length * i->second * codonmatrix(codon, c2) / (*nucmatrix)(n1, n2);
                    }
                }
            }
        }

        const std::map<pair<int, int>, double> &codonpaircount = codonpathsuffstat.GetPairCountMap();
        for (std::map<pair<int, int>, double>::const_iterator i = codonpaircount.begin();
             i != codonpaircount.end(); i++) {
            int cod1 = i->first.first;
            int cod2 = i->first.second;
            int pos = cod->GetDifferingPosition(cod1, cod2);
            if (pos == 3) {
                cerr << "error in codon conj path suffstat\n";
                exit(1);
            }
            int n1 = cod->GetCodonPosition(pos, cod1);
            int n2 = cod->GetCodonPosition(pos, cod2);
            paircount[n1][n2] += i->second;
        }
    }

    //! \brief return the log probability as a function of a nucleotide matrix
    //!
    //! The codon state space is given as an argument (the nucleotide matrix or
    //! \brief return the log probability as a function of a nucleotide matrix
    //!
    //! The codon state space is given as an argument (the nucleotide matrix or
    //! the suff stat themselves do not know the genetic code)
    double GetLogProb(const SubMatrix &mat, const CodonStateSpace &cod) const {
        double total = 0;
        // root part
        int nroot = 0;
        auto rootstat = mat.GetStationary();
        for (int i = 0; i < Nnuc; i++) {
            total += rootcount[i] * log(rootstat[i]);
            nroot += rootcount[i];
        }
        total -= nroot / 3 * log(cod.GetNormStat(rootstat));

        // non root part
        for (int i = 0; i < Nnuc; i++) {
            for (int j = 0; j < Nnuc; j++) {
                if (i != j) {
                    total += paircount[i][j] * log(mat(i, j));
                    total -= pairbeta[i][j] * mat(i, j);
                }
            }
        }

        return total;
    }

    //! add another nucpath suff stat to this
    void Add(const MeanNucPathSuffStat &from) {
        for (int i = 0; i < Nnuc; i++) {
            rootcount[i] += from.rootcount[i];
        }
        for (int i = 0; i < Nnuc; i++) {
            for (int j = 0; j < Nnuc; j++) {
                paircount[i][j] += from.paircount[i][j];
                pairbeta[i][j] += from.pairbeta[i][j];
            }
        }
    }

    //! add another nuc path suffstat to this, operator version
    MeanNucPathSuffStat &operator+=(const MeanNucPathSuffStat &from) {
        Add(from);
        return *this;
    }

  private:
    std::vector<double> rootcount;
    std::vector<vector<double>> paircount;
    std::vector<vector<double>> pairbeta;
};

class NucRatesGeneBranchArray : public SimpleBidimArray<vector<double>> {

    public:

    NucRatesGeneBranchArray(int Ngene, int Nbranch, 
            const GeneBranchGammaEffects& inrhoAC,
            const GeneBranchGammaEffects& inrhoAG,
            // const GeneBranchGammaEffects& inrhoAT,
            const GeneBranchGammaEffects& inrhoCA,
            const GeneBranchGammaEffects& inrhoCG,
            const GeneBranchGammaEffects& inrhoCT) :
        SimpleBidimArray<vector<double>>(Ngene, Nbranch, vector<double>(5,0)),
        rhoAC(inrhoAC), rhoAG(inrhoAG), rhoCA(inrhoCA), rhoCG(inrhoCG), rhoCT(inrhoCT)  {
        // rhoAC(inrhoAC), rhoAG(inrhoAG), rhoAT(inrhoAT), rhoCA(inrhoCA), rhoCG(inrhoCG), rhoCT(inrhoCT)  {
        Update();
    }

    void Update()   {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                Update(i,j);
            }
        }
    }

    void UpdateRow(int i)   {
        for (int j=0; j<GetNcol(); j++) {
            Update(i,j);
        }
    }

    void UpdateColumn(int j)    {
        for (int i=0; i<GetNrow(); i++) {
            Update(i,j);
        }
    }

    void Update(int i, int j)   {
        (*this)(i,j)[0] = rhoAC.GetVal(i,j);
        (*this)(i,j)[1] = rhoAG.GetVal(i,j);
        (*this)(i,j)[2] = 1.0;
        // (*this)(i,j)[2] = rhoAT.GetVal(i,j);
        (*this)(i,j)[3] = rhoCA.GetVal(i,j);
        (*this)(i,j)[4] = rhoCG.GetVal(i,j);
        (*this)(i,j)[5] = rhoCT.GetVal(i,j);
    }

    private:

    const GeneBranchGammaEffects& rhoAC;
    const GeneBranchGammaEffects& rhoAG;
    // const GeneBranchGammaEffects& rhoAT;
    const GeneBranchGammaEffects& rhoCA;
    const GeneBranchGammaEffects& rhoCG;
    const GeneBranchGammaEffects& rhoCT;
};

class NucMatrixGeneBranchArray : public BidimArray<StrandSymmetricIrreversibleSubMatrix>  {

    public:
    NucMatrixGeneBranchArray(int inNgene, int inNbranch, const BidimSelector<vector<double>>& inrates, bool innormalise) :
        Ngene(inNgene), Nbranch(inNbranch), rates(inrates), normalise(innormalise),
        matrixarray(Ngene, vector<StrandSymmetricIrreversibleSubMatrix*>(Nbranch, (StrandSymmetricIrreversibleSubMatrix*) 0)) {
        Create();
    }

    ~NucMatrixGeneBranchArray() {
        Delete();
    }

    int GetNrow() const {
        return Ngene;
    }

    int GetNcol() const {
        return Nbranch;
    }

    const StrandSymmetricIrreversibleSubMatrix &GetVal(int i, int j) const { return *matrixarray[i][j]; }
    StrandSymmetricIrreversibleSubMatrix &operator()(int i, int j) { return *matrixarray[i][j]; }

    void UpdateMatrices()   {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    void UpdateMatrixRow(int i)   {
        for (int j=0; j<GetNcol(); j++) {
            matrixarray[i][j]->CorruptMatrix();
        }
    }

    void UpdateMatricColumn(int j)  {
        for (int i=0; i<GetNrow(); i++) {
            matrixarray[i][j]->CorruptMatrix();
        }
    }
    private:

    void Create()   {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                matrixarray[i][j] = new StrandSymmetricIrreversibleSubMatrix(rates.GetVal(i,j), normalise);
            }
        }
    }
    
    void Delete()	{
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                delete matrixarray[i][j];
            }
	    }
    }

    int Ngene;
    int Nbranch;
    const BidimSelector<vector<double>>& rates;
    bool normalise;
    vector<vector<StrandSymmetricIrreversibleSubMatrix *>> matrixarray;
};

class CodonMatrixGeneBranchArray : public BidimArray<MGOmegaCodonSubMatrix>  {

    public:
    CodonMatrixGeneBranchArray(int inNgene, int inNbranch,
            const CodonStateSpace* incodonstatespace,
            const BidimSelector<StrandSymmetricIrreversibleSubMatrix>& innucmat,
            const GeneBranchGammaEffects& inomega) :
        Ngene(inNgene), Nbranch(inNbranch), 
        codonstatespace(incodonstatespace),
        nucmat(innucmat), omega(inomega),
        matrixarray(Ngene, vector<MGOmegaCodonSubMatrix*>(Nbranch, (MGOmegaCodonSubMatrix*) 0)) {
        Create();
    }

    ~CodonMatrixGeneBranchArray() {
        Delete();
    }

    int GetNrow() const {
        return Ngene;
    }

    int GetNcol() const {
        return Nbranch;
    }

    const MGOmegaCodonSubMatrix &GetVal(int i, int j) const { return *matrixarray[i][j]; }
    MGOmegaCodonSubMatrix &operator()(int i, int j) { return *matrixarray[i][j]; }

    void UpdateMatrices()   {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                matrixarray[i][j]->SetOmega(omega.GetVal(i,j));
                matrixarray[i][j]->CorruptMatrix();
            }
        }
    }

    void UpdateMatrixRow(int i)   {
        for (int j=0; j<GetNcol(); j++) {
            matrixarray[i][j]->SetOmega(omega.GetVal(i,j));
            matrixarray[i][j]->CorruptMatrix();
        }
    }

    void UpdateMatricColumn(int j)  {
        for (int i=0; i<GetNrow(); i++) {
            matrixarray[i][j]->SetOmega(omega.GetVal(i,j));
            matrixarray[i][j]->CorruptMatrix();
        }
    }
    private:

    void Create()   {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                matrixarray[i][j] = new MGOmegaCodonSubMatrix(codonstatespace, 
                        &nucmat.GetVal(i,j), omega.GetVal(i,j));
            }
        }
    }
    
    void Delete()   {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                delete matrixarray[i][j];
            }
        }
    }
    
    int Ngene;
    int Nbranch;
    const CodonStateSpace* codonstatespace;
    const BidimSelector<StrandSymmetricIrreversibleSubMatrix>& nucmat;
    const GeneBranchGammaEffects& omega;
    vector<vector<MGOmegaCodonSubMatrix *>> matrixarray;
};

class PathSuffStatGeneBranchArray : public SimpleBidimArray<PathSuffStat>   {

    public:

    PathSuffStatGeneBranchArray(int Ngene, int Nbranch, int Nstate) : 
        SimpleBidimArray<PathSuffStat>(Ngene, Nbranch, PathSuffStat(Nstate)) {
        Clear();
    }

    void Clear()    {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                (*this)(i,j).Clear();
            }
        }
    }

    void Get(int gene, const NodeSelector<PathSuffStat>& nodearray)   {
        for (int j=0; j<GetNcol(); j++) {
            (*this)(gene,j).Clear();
        }
        RecursiveAdd(gene, nodearray.GetTree().GetRoot(), nodearray);
    }

    void RecursiveAdd(int gene, const Link *from, const NodeSelector<PathSuffStat>& nodearray)    {
        if (!from->isRoot()) {
            (*this)(gene,from->GetBranch()->GetIndex()).Add(nodearray.GetVal(from->GetNode()->GetIndex()));
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursiveAdd(gene, link->Out(), nodearray);
        }
    }

    //! return total log prob over array, given an array of omega_i's of same size
    double GetLogProb(const CodonMatrixGeneBranchArray& mat) const  {
        double total = 0;
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                total += GetVal(i,j).GetLogProb(mat.GetVal(i,j));
            }
        }
        return total;
    }
};

class dSOmegaPathSuffStatGeneBranchArray : public SimpleBidimArray<dSOmegaPathSuffStat> {

    public:

    dSOmegaPathSuffStatGeneBranchArray(int Ngene, int Nbranch) : 
        SimpleBidimArray<dSOmegaPathSuffStat>(Ngene, Nbranch, dSOmegaPathSuffStat()) {
        Clear();
    }

    ~dSOmegaPathSuffStatGeneBranchArray() {}

    void Clear()    {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                (*this)(i,j).Clear();
            }
        }
    }

    void AddSuffStat(const PathSuffStatGeneBranchArray& pathss, const CodonMatrixGeneBranchArray& mat, 
            const GeneBranchGammaEffects& ds, const GeneBranchGammaEffects& om)  {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                (*this)(i,j).AddSuffStat(mat.GetVal(i,j), pathss.GetVal(i,j), 
                        ds.GetVal(i,j), om.GetVal(i,j));
            }
        }
    }

    void AddTo(dSOmegaPathSuffStatBranchArray& to)    {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                to[j].Add((*this)(i,j));
            }
        }
    }
};

class NucPathSuffStatGeneBranchArray : public SimpleBidimArray<MeanNucPathSuffStat> {

    public:

    NucPathSuffStatGeneBranchArray(int Ngene, int Nbranch) : 
        SimpleBidimArray<MeanNucPathSuffStat>(Ngene, Nbranch, MeanNucPathSuffStat()) {
        Clear();
    }

    ~NucPathSuffStatGeneBranchArray() {}

    void Clear()    {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                (*this)(i,j).Clear();
            }
        }
    }

    void AddSuffStat(const PathSuffStatGeneBranchArray& pathss, const CodonMatrixGeneBranchArray& mat)  {
        for (int i=0; i<GetNrow(); i++) {
            for (int j=0; j<GetNcol(); j++) {
                (*this)(i,j).AddSuffStat(mat.GetVal(i,j), pathss.GetVal(i,j));
            }
        }
    }
};

