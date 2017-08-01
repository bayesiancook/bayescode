#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include "Eigen/Dense"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "Random.hpp"

using Vector = double *;
using ConstVect = const double *;
using Matrix = double **;
using EMatrix = Eigen::MatrixXd;
using EVector = Eigen::VectorXd;

class SubMatrix {

  protected:
    // these 2 pure virtual functions are the most essential component of the
    // SubMatrix class
    // see GTRSubMatrix.cpp and CodonSubMatrix.cpp for examples

    // ComputeArray(int state) is in charge of computing the row of the rate
    // matrix
    // corresponding to all possible rates of substitution AWAY from state
    //
    virtual void ComputeArray(int state) const = 0;

    // ComputeStationary() is in charge of computing the vector of stationary
    // probabilities
    // (equilibirum frequencies)
    // of the substitution process
    virtual void ComputeStationary() const = 0;

  public:
    static const int UniSubNmax = 500;
	static int		nunisubcount;
	static int		GetUniSubCount() {return nunisubcount;}

    static int nuni;
    static int nunimax;
    static int diagcount;

    static double GetMeanUni() { return ((double)nunimax) / nuni; }

    SubMatrix(int inNstate, bool innormalise = false);
    virtual ~SubMatrix();

    void Create();

    double operator()(int /*i*/, int /*j*/) const;
    // const double *GetRow(int i) const;
    EVector GetRow(int i) const;

    double Stationary(int i) const;

    // const double *GetStationary() const;
    const EVector &GetEigenStationary() const;
    double *GetStationary() const {
        if (!statflag) {
            UpdateStationary();
        }
        return oldStationary;
    }

    int GetNstate() const { return Nstate; }

    double GetRate() const;
    void ScalarMul(double e);

    bool isNormalised() const { return normalise; }
    void Normalise() const;

    virtual void CorruptMatrix();
    void UpdateMatrix() const;

    void ActivatePowers() const;
    void InactivatePowers() const;
    double Power(int n, int i, int j) const;
    double GetUniformizationMu() const;

    virtual void ToStream(std::ostream &os) const;
    void CheckReversibility() const;

    int GetDiagStat() const { return ndiagfailed; }

    void BackwardPropagate(const double *up, double *down, double length) const;
    void ForwardPropagate(const double *down, double *up, double length) const;
    // virtual void     FiniteTime(int i0, double* down, double length);

    // double **GetQ() const { return Q; }

    static double meanz;
    static double maxz;
    static double nz;

	// uniformization resampling methods
	// CPU level 1
	void GetFiniteTimeTransitionProb(int state, double* down, double efflength) const;
	double 			GetFiniteTimeTransitionProb(int stateup, int statedown, double efflength) const;
	int 			DrawUniformizedTransition(int state, int statedown, int n) const;
	int 			DrawUniformizedSubstitutionNumber(int stateup, int statedown, double efflength) const;
	//

	// used by accept-reject resampling method
	// CPU level 1
	int 			DrawOneStep(int state) const;
	double			DrawWaitingTime(int state) const;
	int 			DrawFromStationary() const;

	// double SuffStatLogProb(SuffStat* suffstat);

  protected:
    void UpdateRow(int state) const;
    void UpdateStationary() const;

    void ComputePowers(int N) const;
    void CreatePowers(int n) const;

    bool ArrayUpdated() const;

    int Diagonalise() const;

    // data members

    mutable bool powflag;
    mutable bool diagflag;
    mutable bool statflag;
    mutable bool *flagarray;

    int Nstate;
    mutable int npow;
    mutable double UniMu;

    double ***mPow;

    // Q : the infinitesimal generator matrix
    // mutable double **Q;
    mutable EMatrix Q;            // Q : the infinitesimal generator matrix

    // the stationary probabilities of the matrix
    // mutable double *mStationary;
    mutable EVector mStationary;  // the stationary probabilities of the matrix
    double *oldStationary;

    mutable Eigen::EigenSolver<EMatrix> solver;

    bool normalise;

    // an auxiliary matrix
    mutable double **aux;

  protected:

    // v : eigenvalues
    // vi : imaginary part
    // u : the matrix of eigen vectors
    // invu : the inverse of u

    /*
    mutable double **u;
    mutable double **invu;
    mutable double *v;
    mutable double *vi;
    */

    mutable EMatrix u;     // u : the matrix of eigen vectors
    mutable EMatrix invu;  // invu : the inverse of u
    mutable EVector v;     // v : eigenvalues
    mutable EVector vi;    // vi : imaginary part

    mutable int ndiagfailed;
};

//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------

inline double SubMatrix::operator()(int i, int j) const {
    if (!flagarray[i]) {
        UpdateRow(i);
    }
    return Q(i, j);
    // return Q[i][j];
}

inline EVector SubMatrix::GetRow(int i) const {
    if (!flagarray[i]) {
        UpdateRow(i);
    }
    return Q.row(i);
}

/*
inline const double *SubMatrix::GetRow(int i) const {
    if (!flagarray[i]) {
        UpdateRow(i);
    }
    return Q[i];
}
*/

inline const EVector &SubMatrix::GetEigenStationary() const {
// inline const double *SubMatrix::GetStationary() const {
    if (!statflag) {
        UpdateStationary();
    }
    return mStationary;
}

inline double SubMatrix::Stationary(int i) const {
    if (!statflag) {
        UpdateStationary();
    }
    return mStationary[i];
}

inline void SubMatrix::CorruptMatrix() {
    diagflag = false;
    statflag = false;
    for (int k = 0; k < Nstate; k++) {
        flagarray[k] = false;
    }
    InactivatePowers();
}

inline bool SubMatrix::ArrayUpdated() const {
    bool qflag = true;
    for (int k = 0; k < Nstate; k++) {
        qflag &= static_cast<int>(flagarray[k]);
    }
    return qflag;
}

inline void SubMatrix::UpdateStationary() const {
    ComputeStationary();
    for (int i = 0; i < Nstate; i++) {
        oldStationary[i] = mStationary[i];
    }
    statflag = true;
}

inline void SubMatrix::UpdateRow(int state) const {
    if (isNormalised()) {
        UpdateMatrix();
    } else {
        if (!statflag) {
            UpdateStationary();
        }
        ComputeArray(state);
        flagarray[state] = true;
    }
}

inline void SubMatrix::BackwardPropagate(const double *up, double *down, double length) const {
    if (! diagflag)  {
        Diagonalise();
    }

    int matSize = GetNstate();

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] = 0;
    }
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += invu(i,j) * up[j];
            // aux[i] += invu[i][j] * up[j];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] *= exp(length * v[i]);
    }

    for (int i = 0; i < GetNstate(); i++) {
        down[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            down[i] += u(i,j) * aux[j];
            // down[i] += u[i][j] * aux[j];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        if (std::isnan(down[i])) {
            std::cerr << "error in back prop\n";
            for (int j = 0; j < GetNstate(); j++) {
                std::cerr << up[j] << '\t' << down[j] << '\t' << Stationary(j) << '\n';
            }
            exit(1);
        }
    }
    double maxup = 0;
    for (int k = 0; k < matSize; k++) {
        if (up[k] < 0) {
            std::cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
            // down[k] = 0;
        }
        if (maxup < up[k]) {
            maxup = up[k];
        }
    }
    double max = 0;
    for (int k = 0; k < matSize; k++) {
        if (down[k] < 0) {
            down[k] = 0;
        }
        if (max < down[k]) {
            max = down[k];
        }
    }
    if (maxup == 0) {
        std::cerr << "error in backward propagate: null up array\n";
        exit(1);
    }
    if (max == 0) {
        std::cerr << "error in backward propagate: null array\n";
        for (int k = 0; k < matSize; k++) {
            std::cerr << up[k] << '\t' << down[k] << '\n';
        }
        std::cerr << '\n';
        exit(1);
    }
    down[matSize] = up[matSize];

    delete[] aux;
}

inline void SubMatrix::ForwardPropagate(const double *down, double *up, double length) const {

    if (! diagflag)  {
        Diagonalise();
    }

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += down[j] * u(j,i);
            // aux[i] += down[j] * u[j][i];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] *= exp(length * v[i]);
    }

    for (int i = 0; i < GetNstate(); i++) {
        up[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            up[i] += aux[j] * invu(j,i);
            // up[i] += aux[j] * invu[j][i];
        }
    }

    delete[] aux;
}

inline double SubMatrix::GetFiniteTimeTransitionProb(int stateup, int statedown, double efflength)	const {

    if (! diagflag) {
        Diagonalise();
    }

	double tot = 0;
	for (int i=0; i<GetNstate(); i++)	{
		tot += invu(stateup,i) * exp(efflength * v[i]) * invu(i,statedown);
		// tot += invu[stateup][i] * exp(efflength * v[i]) * invu[i][statedown];
	}
	return tot;
}

inline void SubMatrix::GetFiniteTimeTransitionProb(int state, double* p, double efflength)	const {

	double* p1 = new double[GetNstate()];
	for (int k=0; k<GetNstate(); k++)	{
		p1[k] = 0;
	}
	p1[state] = 1;
	ForwardPropagate(p1,p,efflength);
	double tot = 0;
	for (int k=0; k<GetNstate(); k++)	{
		tot += p[k];
	}
	delete[] p1;
	if (fabs(1 - tot) > 1e-5)	{
		std::cerr << "error in forward propagate: normalization : " << tot << '\t' << fabs(1 - tot) << '\n';
		std::cerr << "eff length : " << efflength << '\n';
		ToStream(std::cerr);
		exit(1);
	}
}

inline int SubMatrix::DrawUniformizedTransition(int state, int statedown, int n)	const {

	double* p = new double[GetNstate()];
	double tot = 0;
	for (int l=0; l<GetNstate(); l++)	{
		tot += Power(1,state,l) * Power(n,l,statedown);
		p[l] = tot;
	}

	double s = tot * Random::Uniform();
	int k = 0;
	while ((k<GetNstate()) && (s > p[k]))	{
		k++;
	}
	delete[] p;
	if (k == GetNstate())	{
		std::cerr << "error in DrawUniformizedTransition: overflow\n";
		throw;
	}
	return k;
}

inline int SubMatrix::DrawUniformizedSubstitutionNumber(int stateup, int statedown, double efflength)	const {

	double mu = GetUniformizationMu();
	double fact = exp(- efflength * mu);
	int m = 0;
	double total = (stateup==statedown) * fact;
	double Z = GetFiniteTimeTransitionProb(stateup,statedown,efflength);
	double q = Random::Uniform() * Z;
	
	while ((m<UniSubNmax) && (total < q)) 	{
		m++;
		fact *= mu * efflength / m;
		total += Power(m,stateup,statedown) * fact;
		if ((total-Z)>1e-12)	{
			std::cerr << "error in DrawUniformizedSubstitutionNumber: normalising constant\n";
			std::cerr << total << '\t' << Z << '\n';
			std::cerr << mu << '\n';
			std::cerr << m << '\n';
			std::cerr << stateup << '\t' << statedown << '\n';

			ToStream(std::cerr);
			CheckReversibility();
			throw;
		}
	}
	if (m == UniSubNmax)	{
		nunisubcount++;
	}
	return m;
}

inline double SubMatrix::DrawWaitingTime(int state)	const {

    EVector row = GetRow(state);
	// const double* row = GetRow(state);
	double t = Random::sExpo() / (-row[state]);
	return t;
}

inline int SubMatrix::DrawOneStep(int state)	const {

    EVector row = GetRow(state);
	double p = -row[state] * Random::Uniform();
	int k = -1;
	double tot = 0;
	do	{
		k++;
		if (k != state)	{
			tot += row[k];
		}
	}
	while ((k<GetNstate()) && (tot < p));
	if (k == GetNstate())	{
		std::cerr << "error in DrawOneStep\n";
		std::cerr << GetNstate() << '\n';
		for (int k=0; k<GetNstate(); k++)	{
			std::cerr << row[k] << '\n';
		}
		exit(1);
	}
	return k;
}

inline int SubMatrix::DrawFromStationary()	const {

	double p = Random::Uniform();
	int k = 0;
	double tot = mStationary[k];
	while ((k<GetNstate()) && (tot < p))	{
		k++;
		if (k == GetNstate())	{
			std::cerr << "error in DrawFromStationary\n";
			double tot = 0;
			for (int l=0; l<GetNstate(); l++)	{
				std::cerr << mStationary[l] << '\t';
				tot += mStationary[l];
			}
			std::cerr << '\n';
			std::cerr << "total : " << tot << '\t' << 1-tot << '\n';
			exit(1);
		}
		tot += mStationary[k];
	}
	return k;
}

#endif  // SUBMATRIX_H
