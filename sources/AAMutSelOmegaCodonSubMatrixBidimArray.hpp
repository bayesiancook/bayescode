
#pragma once


/**
 * \brief An array of mutation-selection codon matrices (with omega)
 *
 * The array takes a single nucleotide matrix, an array of amino-acid fitness
 * profiles and either a single omega or an array of omega values and construct
 * mutation-selection matrices accordingly.
 */

class AAMutSelOmegaCodonSubMatrixBidimArray : public BidimArray<SubMatrix>,
                                         public BidimArray<AAMutSelOmegaCodonSubMatrix> {
  public:
    //! constructor with a nucleotide matrix, an array of amino-acid fitness
    //! profiles and a single omega value (for all matrices)
    AAMutSelOmegaCodonSubMatrixBidimArray(const CodonStateSpace *incodonstatespace,
                                     const SubMatrix *innucmatrix,
                                     const Selector<vector<double>> *inaafitnessarray,
                                     const Selector<double> *inomegaarray)
        : codonstatespace(incodonstatespace),
          nucmatrix(innucmatrix),
          aafitnessarray(inaafitnessarray),
          omegaarray(inomegaarray),
          matrixbidimarray(inaafitnessarray->GetSize(), std::vector<AAMutSelOmegaCodonSubMatrix*>(inomegaarray->GetSize())) {
        Create();
    }

    ~AAMutSelOmegaCodonSubMatrixBidimArray() { Delete(); }

    //! return array size
    int GetNrow() const { return aafitnessarray->GetSize(); }
    int GetNcol() const { return omegaarray->GetSize(); }

    //! const access to matrix i
    const AAMutSelOmegaCodonSubMatrix &GetVal(int i, int j) const { return *matrixbidimarray[i][j]; }
    //! non-const access to matrix i
    AAMutSelOmegaCodonSubMatrix &operator()(int i, int j) { return *matrixbidimarray[i][j]; }

    //! const acess to nucleotide matrix
    const SubMatrix &GetNucMatrix() const { return *nucmatrix; }

    //! update all matrices
    void UpdateCodonMatrices() {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) {
                    (*this)(i,j).SetOmega(omegaarray->GetVal(j));
                    (*this)(i,j).CorruptMatrix();
            }
        }
    }

  private:
    void Create() {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) {
                matrixbidimarray[i][j] = new AAMutSelOmegaCodonSubMatrix(codonstatespace, nucmatrix,
                                                                 aafitnessarray->GetVal(i),
                                                                 omegaarray->GetVal(j), 1.0);
	    }
        }
    }

    void Delete() {
        for (int i = 0; i < this->GetNrow(); i++) {
            for (int j = 0; j < this->GetNcol(); j++) {
                delete matrixbidimarray[i][j];
            }
        }
    }

    const CodonStateSpace *codonstatespace;
    const SubMatrix *nucmatrix;
    const Selector<vector<double>> *aafitnessarray;
    const Selector<double> *omegaarray;
    vector<vector<AAMutSelOmegaCodonSubMatrix *>> matrixbidimarray;
};

