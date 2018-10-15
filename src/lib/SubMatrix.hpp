#pragma once

// #include "Eigen/Dense"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "Random.hpp"

// using EMatrix = Eigen::MatrixXd;
// using EVector = Eigen::VectorXd;
//

/** \brief A generic substitution matrix class.
 *
 * SubMatrix represents the generator a Nstate*Nstate Markovian substitution
 * process. It implements most of the standard methods that are necessary for
 * likelihood calculation and stochastic character mapping along phylogenies
 * (diagonalisation, exponentiation, drawing a waiting time until next event,
 * choosing next event conditional on current event, etc). Specific substitution
 * processes can be implemented by deriving them from this abstact class (see
 * for instance GTRSubMatrix or MGOmegaCodonSubMatrix). This essentially
 * requires to implement the two following pure virtual functions of SubMatrix:
 * ComputeArray(int s), which is in charge of computing row s of the rate
 * matrix, and ComputeStationary(), which should calculate the equilbrium
 * frequencies of the process.
 */

class SubMatrix {
  public:
    //! \brief base constructor: should always give the dimension (Nstate) of the
    //! matrix and whether we want it to be normalized
    SubMatrix(int inNstate, bool innormalise = false);
    virtual ~SubMatrix();

    //! const access to entry at ith row and jth column (checked for current
    //! update status, see UpdateMatrix)
    double operator()(int /*i*/, int /*j*/) const;
    //! const access to a row of the matrix (checked for current update status)
    EVector GetRow(int i) const;

    //! const access to equilbrium frequency of state i (checked for current
    //! update status)
    double Stationary(int i) const;

    //! const access to equilbrium frequency vector (checked for current update
    //! status)
    const EVector &GetStationary() const;

    //! dimension of the statespace
    int GetNstate() const { return Nstate; }

    //! compute the row of the generator corresponding to rates away from given
    //! state (should be defined in derived classes)
    virtual void ComputeArray(int state) const = 0;

    //! compute the vector of stationary probabilities (equilibrium frequencies --
    //! should be defined in derived classes)
    virtual void ComputeStationary() const = 0;

    //! whether or not the matrix is normalized (i.e. the expected rate of
    //! substitution per unit of time is 1)
    bool isNormalised() const { return normalise; }

    //! normalize the matrix
    void Normalise() const;

    //! get the normalization factor
    double GetRate() const;

    //! multiply all entries by a scalar e
    void ScalarMul(double e);

    //! set the flags telling the matrix that it should recalculate its rates and
    //! all dependent variables
    virtual void CorruptMatrix();

    //! \brief recalculate all rates and dependent variables
    //!
    //! access to rates, equilibrium frequencies or exponentiation/diagonalisation
    //! variables (see methods GetRow, GetStationary, etc) is first checked for
    //! correct update of all rates and dependent variables, as follows: if
    //! CorruptMatrix has been called before the last access, then, corrupt flags
    //! are active, und Update functions are called automatically before giving
    //! access. Note that update can in practice be partial (only the rows that
    //! are required are recalculated). This is implemented by way of an array of
    //! update flags.
    void UpdateMatrix() const;

    //! a simple output stream function (mostly useful for tracing and debugging)
    virtual void ToStream(std::ostream &os) const;

    //! propagate a conditional likelihood vector in the tip-to-root direction
    //! (pruning algorithm)
    void BackwardPropagate(const double *up, double *down, double length) const;
    //! propagate a conditional likelihood vector in the root-to-tip direction
    //! (pruning algorithm)
    void ForwardPropagate(const double *down, double *up, double length) const;

    //! get vector of finite time transition probabilities from given state to all
    //! possible states down, along branch of efflength=length*rate
    void GetFiniteTimeTransitionProb(int state, double *down, double efflength) const;

    //! return the transition probability between stateup and statedown along
    //! branch of efflength = length*rate
    double GetFiniteTimeTransitionProb(int stateup, int statedown, double efflength) const;
    //! draw the uniformized number of transitions along the branch, conditional
    //! on begin and end states
    int DrawUniformizedSubstitutionNumber(int stateup, int statedown, double efflength) const;
    //! draw the state of the next event, given current state and given that total
    //! number of events until reaching statedown is n
    int DrawUniformizedTransition(int state, int statedown, int n) const;

    //! draw state of next event given current state
    int DrawOneStep(int state) const;
    //! draw state after total time, given current state
    int DrawFiniteTime(int state, double time) const;
    //! draw waiting time until next event, given current state
    double DrawWaitingTime(int state) const;
    //! draw state from equilibrium frequencies
    int DrawFromStationary() const;

  protected:
    static const int UniSubNmax = 500;
    static int nunisubcount;
    static int GetUniSubCount() { return nunisubcount; }

    static int nuni;
    static int nunimax;
    static int diagcount;
    static double diagerr;

    static double GetMeanUni() { return ((double)nunimax) / nuni; }

    void Create();

    void ActivatePowers() const;
    void InactivatePowers() const;
    double Power(int n, int i, int j) const;
    double GetUniformizationMu() const;

    void CheckReversibility() const;

    int GetDiagStat() const { return ndiagfailed; }

    static double meanz;
    static double maxz;
    static double nz;

    void UpdateRow(int state) const;
    void UpdateStationary() const;

    void ComputePowers(int N) const;
    void CreatePowers(int n) const;

    bool ArrayUpdated() const;

    int Diagonalise() const;
    int EigenDiagonalise() const;
    int OldDiagonalise() const;
    double CheckDiag() const;

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
    mutable double **ptrQ;
    mutable EMatrix Q;  // Q : the infinitesimal generator matrix

    // the stationary probabilities of the matrix
    mutable double *ptrStationary;
    mutable EVector mStationary;  // the stationary probabilities of the matrix

    mutable Eigen::EigenSolver<EMatrix> solver;

    bool normalise;

    // an auxiliary matrix
    mutable double **aux;

  protected:
    mutable double **ptru;
    mutable double **ptrinvu;
    mutable double *ptrv;

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
    if (!flagarray[i]) { UpdateRow(i); }
    return Q(i, j);
    // return Q[i][j];
}

inline EVector SubMatrix::GetRow(int i) const {
    if (!flagarray[i]) { UpdateRow(i); }
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

inline const EVector &SubMatrix::GetStationary() const {
    // inline const double *SubMatrix::GetStationary() const {
    if (!statflag) { UpdateStationary(); }
    return mStationary;
}

inline double SubMatrix::Stationary(int i) const {
    if (!statflag) { UpdateStationary(); }
    return mStationary[i];
}

inline void SubMatrix::CorruptMatrix() {
    diagflag = false;
    statflag = false;
    for (int k = 0; k < Nstate; k++) { flagarray[k] = false; }
    InactivatePowers();
}

inline bool SubMatrix::ArrayUpdated() const {
    bool qflag = true;
    for (int k = 0; k < Nstate; k++) { qflag &= static_cast<int>(flagarray[k]); }
    return qflag;
}

inline void SubMatrix::UpdateStationary() const {
    ComputeStationary();
    statflag = true;
}

inline void SubMatrix::UpdateRow(int state) const {
    if (isNormalised()) {
        UpdateMatrix();
    } else {
        if (!statflag) { UpdateStationary(); }
        ComputeArray(state);
        flagarray[state] = true;
    }
}

inline void SubMatrix::BackwardPropagate(const double *up, double *down, double length) const {
    if (!diagflag) { Diagonalise(); }

    int matSize = GetNstate();

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) { aux[i] = 0; }
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += invu(i, j) * up[j];
            // aux[i] += invu[i][j] * up[j];
        }
    }

    for (int i = 0; i < GetNstate(); i++) { aux[i] *= exp(length * v[i]); }

    for (int i = 0; i < GetNstate(); i++) { down[i] = 0; }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            down[i] += u(i, j) * aux[j];
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
        if (maxup < up[k]) { maxup = up[k]; }
    }
    double max = 0;
    for (int k = 0; k < matSize; k++) {
        if (down[k] < 0) { down[k] = 0; }
        if (max < down[k]) { max = down[k]; }
    }
    if (maxup == 0) {
        std::cerr << "error in backward propagate: null up array\n";
        exit(1);
    }
    if (max == 0) {
        std::cerr << "error in backward propagate: null array\n";
        for (int k = 0; k < matSize; k++) { std::cerr << up[k] << '\t' << down[k] << '\n'; }
        std::cerr << '\n';
        exit(1);
    }
    down[matSize] = up[matSize];

    delete[] aux;
}

inline void SubMatrix::ForwardPropagate(const double *down, double *up, double length) const {
    if (!diagflag) { Diagonalise(); }

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) { aux[i] = 0; }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += down[j] * u(j, i);
            // aux[i] += down[j] * u[j][i];
        }
    }

    for (int i = 0; i < GetNstate(); i++) { aux[i] *= exp(length * v[i]); }

    for (int i = 0; i < GetNstate(); i++) { up[i] = 0; }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            up[i] += aux[j] * invu(j, i);
            // up[i] += aux[j] * invu[j][i];
        }
    }

    delete[] aux;
}

inline double SubMatrix::GetFiniteTimeTransitionProb(
    int stateup, int statedown, double efflength) const {
    if (!diagflag) { Diagonalise(); }

    double tot = 0;
    for (int i = 0; i < GetNstate(); i++) {
        tot += u(stateup, i) * exp(efflength * v[i]) * invu(i, statedown);
        // tot += invu[stateup][i] * exp(efflength * v[i]) * invu[i][statedown];
    }
    return tot;
}

inline void SubMatrix::GetFiniteTimeTransitionProb(int state, double *p, double efflength) const {
    double *p1 = new double[GetNstate()];
    for (int k = 0; k < GetNstate(); k++) { p1[k] = 0; }
    p1[state] = 1;
    ForwardPropagate(p1, p, efflength);
    double tot = 0;
    for (int k = 0; k < GetNstate(); k++) { tot += p[k]; }
    delete[] p1;
    if (fabs(1 - tot) > 1e-4) {
        std::cerr << "error in forward propagate: normalization : " << tot << '\t' << fabs(1 - tot)
                  << '\n';
        std::cerr << "eff length : " << efflength << '\n';
        std::cerr << "diag error: " << CheckDiag() << '\n';
        // ToStream(std::cerr);
        exit(1);
    }
}

inline int SubMatrix::DrawUniformizedTransition(int state, int statedown, int n) const {
    double *p = new double[GetNstate()];
    double tot = 0;
    for (int l = 0; l < GetNstate(); l++) {
        tot += Power(1, state, l) * Power(n, l, statedown);
        p[l] = tot;
    }

    double s = tot * Random::Uniform();
    int k = 0;
    while ((k < GetNstate()) && (s > p[k])) { k++; }
    delete[] p;
    if (k == GetNstate()) {
        std::cerr << "error in DrawUniformizedTransition: overflow\n";
        throw;
    }
    return k;
}

inline int SubMatrix::DrawUniformizedSubstitutionNumber(
    int stateup, int statedown, double efflength) const {
    double mu = GetUniformizationMu();
    double fact = exp(-efflength * mu);
    int m = 0;
    double total = (stateup == statedown) * fact;
    double Z = GetFiniteTimeTransitionProb(stateup, statedown, efflength);
    double q = Random::Uniform() * Z;

    while ((m < UniSubNmax) && (total < q)) {
        m++;
        fact *= mu * efflength / m;
        total += Power(m, stateup, statedown) * fact;
        if ((total - Z) > 1e-5) {
            std::cerr << "error in DrawUniformizedSubstitutionNumber: normalising "
                         "constant\n";
            std::cerr << total << '\t' << Z << '\n';
            std::cerr << mu << '\n';
            std::cerr << m << '\n';
            std::cerr << stateup << '\t' << statedown << '\n';

            CheckReversibility();
            std::cerr << "diag error: " << CheckDiag() << '\n';
            // ToStream(std::cerr);
            throw;
        }
    }
    if (m == UniSubNmax) { nunisubcount++; }
    return m;
}

inline int SubMatrix::DrawFiniteTime(int initstate, double time) const {
    double t = 0;
    int s = initstate;
    while (t < time) {
        double dt = DrawWaitingTime(s);
        t += dt;
        if (t < time) { s = DrawOneStep(s); }
    }
    return s;
}

inline double SubMatrix::DrawWaitingTime(int state) const {
    EVector row = GetRow(state);
    // const double* row = GetRow(state);
    double t = Random::sExpo() / (-row[state]);
    return t;
}

inline int SubMatrix::DrawOneStep(int state) const {
    EVector row = GetRow(state);
    double p = -row[state] * Random::Uniform();
    int k = -1;
    double tot = 0;
    do {
        k++;
        if (k != state) { tot += row[k]; }
    } while ((k < GetNstate()) && (tot < p));
    if (k == GetNstate()) {
        std::cerr << "error in DrawOneStep\n";
        std::cerr << GetNstate() << '\n';
        for (int k = 0; k < GetNstate(); k++) { std::cerr << row[k] << '\n'; }
        exit(1);
    }
    return k;
}

inline int SubMatrix::DrawFromStationary() const {
    double p = Random::Uniform();
    int k = 0;
    double tot = mStationary[k];
    while ((k < GetNstate()) && (tot < p)) {
        k++;
        if (k == GetNstate()) {
            std::cerr << "error in DrawFromStationary\n";
            double tot = 0;
            for (int l = 0; l < GetNstate(); l++) {
                std::cerr << mStationary[l] << '\t';
                tot += mStationary[l];
            }
            std::cerr << '\n';
            std::cerr << "total : " << tot << '\t' << 1 - tot << '\n';
            exit(1);
        }
        tot += mStationary[k];
    }
    return k;
}
