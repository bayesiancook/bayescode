#pragma once

#include <iostream>
#include "linalg.hpp"
#include "Random.hpp"

class CovMatrix {

	public:

	static double GetMeanDiagError()	{
		if (ndiag)	{
			return toterror / ndiag;
		}
		return 0;
	}

	static double GetMaxDiagError()	{
		return maxerror;
	}

	static double GetMaxDiagN()	{
		return maxn;
	}

	CovMatrix(int indim) :
        value(indim, vector<double>(indim, 0)), 
        invvalue(indim, vector<double>(indim, 0)), 
        u(indim, vector<double>(indim, 0)), 
        invu(indim, vector<double>(indim, 0)), 
        v(indim,0),
        logv(indim,0),
        diagflag(false) {
            for (int i=0; i<GetDim(); i++)  {
                value[i][i] = 1.0;
            }
	}

    int GetDim() const {
        return int(value.size());
    }

	double	operator()(int i, int j) const {
		return value[i][j];
	}

    void setval(int i, int j, double val)   {
        value[i][j] = val;
    }

    void add(int i, int j, double val)  {
        value[i][j] += val;
    }

	CovMatrix& operator*=(double scal) {
		for (int i=0; i<GetDim(); i++){
            for (int j=0; j<GetDim(); j++)  {
                value[i][j] *= scal;
            }
        }
		diagflag = false;
		return *this;
	}

	virtual void MultivariateNormalSample(vector<double>& vec, bool inverse=false) const {
        vector<double> principalcomp(GetDim(), 0);
		for (int i=0; i<GetDim(); i++) {
            if (inverse)    {
                principalcomp[i] =  Random::sNormal() * sqrt(GetEigenVal()[i]);
            }
            else    {
                principalcomp[i] =  Random::sNormal() / sqrt(GetEigenVal()[i]);
            }
		}
		for (int i=0; i<GetDim(); i++) {
			vec[i]=0;
        }
		for (int i=0; i<GetDim(); i++) {
			for (int j=0; j<GetDim(); j++) {
				vec[j] +=  principalcomp[i] * GetEigenVect()[j][i];
			}
		}
	}

	double logMultivariateNormalDensity(const vector<double>& dval) const {

		double tXSX = 0;
		for (int i=0; i<GetDim() ; i++) {
			tXSX += GetInvMatrix()[i][i] * dval[i] * dval[i];
			for (int j=0; j<i ; j++) {
				tXSX += 2 * GetInvMatrix()[i][j] * dval[j] * dval[i];
			}
		}
		return -0.5 * (GetLogDeterminant() + tXSX);
	}

	const vector<double>& GetEigenVal() const {
		if (! diagflag) Diagonalise();
		return v;
	}

	const vector<double>& GetLogEigenVal() const {
		if (! diagflag) Diagonalise();
		return logv;
	}

	const vector<vector<double>>& GetEigenVect() const {
		if (! diagflag)	{
			Diagonalise();
		}

		return u;
	}

	double GetLogDeterminant() const {
		double ret = 0;
		for (int i=0; i<GetDim(); i++) {
			ret += GetLogEigenVal()[i];
		}
		return ret;
	}

	const vector<vector<double>>& GetInvEigenVect() const {
		if (! diagflag) Diagonalise();
		return invu;
	}

	void SetToZero()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}
        diagflag = false;
	}

	void SetToIdentity()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				value[i][j] = 0;
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			value[i][i] = 1;
		}
        diagflag = false;
	}

	void Project(int index, double** m)	{

		int k = 0;
		for (int i=0; i<GetDim(); i++)	{
			if (i != index)	{
				int l = 0;
				for (int j=0; j<GetDim(); j++)	{
					if (j != index)	{
						m[k][l] = value[i][j] - value[i][index] * value[j][index] / value[index][index];
						l++;
					}
				}
				k++;
			}
		}
	}

	const vector<vector<double>>& GetMatrix() const {
		return value;
	}

	const vector<vector<double>>& GetInvMatrix() const {
		if (! diagflag) Diagonalise();
		return invvalue;
	}

	bool isPositive() const {
		if (! diagflag) Diagonalise();
		bool r = true;
		for (int i=0; i<GetDim(); i++){
			if(GetEigenVal()[i] <= 1e-6){
				r = false;
			}
		}
		return r;
	}

	void CorruptDiag() const {
		diagflag = false;
	}

	double GetMax() const {
		double max = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tmp = fabs(value[i][j]);
				if (max < tmp)	{
					max = tmp;
				}
			}
		}
		return max;
	}

	//Set the matrix to it s inverse //loook si diagflag
	int Invert() {
        vector<vector<double>> a(GetDim(), vector<double>(GetDim(), 0));

		// copy value into a :
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				a[i][j] = value[i][j];
			}
		}

		// invert a into value
		// InvertMatrix(a, GetDim(), w, iw, value);
		double logdet = LinAlg::Gauss(a,GetDim(),value);

		// cerr << "check inverse : " << CheckInverse() << '\n';
		diagflag = false;
		if (std::isinf(logdet))	{
			cerr << "error in cov matrix: non invertible\n";
			return 1;
			exit(1);
		}
		return 0;
	}

	double CheckInverse() const {
		double max = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += value[i][k] * GetInvMatrix()[k][j];
				}
				if (i == j)	{
					tot --;
				}
				if (max < fabs(tot))	{
					max = fabs(tot);
				}
			}
		}
		return max;
	}

	int Diagonalise() const {

		int nmax = 1000;
		double epsilon = 1e-10;

		int n = LinAlg::DiagonalizeSymmetricMatrix(value,GetDim(),nmax,epsilon,v,u);
		if (maxn < n)	{
			maxn = n;
		}
		bool failed = (n == nmax);
		if (failed)	{
			cerr << "diag failed\n";
			cerr << n << '\n';
			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					cerr << value[i][j] << '\t';
				}
				cerr << '\n';
			}
			exit(1);
		}

		// normalise u
		for (int i=0; i<GetDim(); i++)	{
			double total = 0;
			for (int j=0; j<GetDim(); j++)	{
				total += u[j][i] * u[j][i];
			}
		}
		// u-1 = tu
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				invu[j][i] = u[i][j];
			}
		}

		for (int i=0; i<GetDim(); i++)	{
			logv[i] = log(v[i]);
		}

		LinAlg::Gauss(value,GetDim(),invvalue);

		diagflag = true;
		double tmp = CheckDiag();
		if (maxerror < tmp)	{
			maxerror = tmp;
		}
		toterror += tmp;
		ndiag ++;
		return failed;
	}

	double CheckDiag() const {
        vector<vector<double>> a(GetDim(), vector<double>(GetDim(), 0));
        vector<vector<double>> b(GetDim(), vector<double>(GetDim(), 0));

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += invu[i][k] * GetMatrix()[k][j];
				}
				a[i][j] = tot;
			}
		}

		double max = 0;

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tot = 0;
				for (int k=0; k<GetDim(); k++)	{
					tot += a[i][k] * u[k][j];
				}
				b[i][j] = tot;
				if (i != j)	{
					if (max < fabs(tot))	{
						max = fabs(tot);
					}
				}
			}
		}
		return max;
	}

	void SetToScatterMatrix(const vector<vector<double>>& invals)    {
        int df = invals.size();
		for (int i=0; i<GetDim(); i++) {
			for (int j=0; j<GetDim(); j++) {
				value[i][j] = 0;
				for (int l=0; l<df; l++) {
					value[i][j] += invals[l][i] * invals[l][j];
				}
			}
		}
		diagflag = false;
	}

    /*
	void PrintCorrelationCoefficients(ostream& os)	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i==j) 	{
					os << "\t-";
				}
				else	{
					os << '\t' << GetMatrix()[i][j] / sqrt(GetMatrix()[i][i] * GetMatrix()[j][j]);
				}
			}
			os << '\n';
		}
	}


	void PrintEigenVectors(ostream& os)	{
		os << "val";
		for (int j=0; j<GetDim(); j++)	{
			os << '\t' << j;
		}
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			os << v[i] << '\t';
			for (int j=0; j<GetDim(); j++)	{
				os << '\t' << u[i][j];
			}
			os << '\n';
		}
		os << '\n';
		os << "inverse eigenvector matrix\n";
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << '\t' << invu[i][j];
			}
			os << '\n';
		}
		os << '\n';

		// proportion of variance explained
		double** prop = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			prop[i] = new double[GetDim()];
		}
		for (int i=0; i<GetDim(); i++)	{
			double total = 0;
			for (int j=0; j<GetDim(); j++)	{
				prop[i][j] = invu[j][i] * invu[j][i] * v[j];
				// prop[i][j] = u[i][j] * u[i][j] * v[j];
				// prop[j][i] = invu[j][i] * invu[j][i] * v[j];
				total += prop[i][j];
			}
			for (int j=0; j<GetDim(); j++)	{
				prop[i][j] /= total;
			}
		}
		os << "proportion of variance explained\n";
		for (int i=0; i<GetDim(); i++)	{
			os << i;
			for (int j=0; j<GetDim(); j++)	{
				os << '\t' << prop[i][j];
			}
			os << '\n';
		}
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			delete[] prop[i];
		}
		delete[] prop;

	}
    */


	friend ostream& operator<<(ostream& os, const CovMatrix& r)  {
		for (int i=0; i<r.GetDim(); i++)	{
			for (int j=0; j<r.GetDim(); j++)	{
				os << r.GetMatrix()[i][j] << '\t';
			}
		}
		return os;
	}

	friend istream& operator>>(istream& is, CovMatrix& r)  {
		for (int i=0; i<r.GetDim(); i++)	{
			for (int j=0; j<r.GetDim(); j++)	{
				is >> r.value[i][j];
			}
		}
		return is;
	}

	private:

	vector<vector<double>> value;
    // inverse
	mutable vector<vector<double>> invvalue;
    // eigenvect
	mutable vector<vector<double>> u;
    // inv eigenvect
	mutable vector<vector<double>> invu;
    // eigenval
    mutable vector<double> v;
    // log eigenval
    mutable vector<double> logv;
	mutable bool diagflag;

	static int maxn;
	static int ndiag;
	static double maxerror;
	static double toterror;
};

