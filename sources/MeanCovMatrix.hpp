
#pragma once

#include "CovMatrix.hpp"

class MeanCovMatrix {

	public:

	MeanCovMatrix(int indim) : tex(false), dim(indim), size(0) {

        val.assign(dim, vector<vector<double>>(dim, vector<double>()));
        invval.assign(dim, vector<vector<double>>(dim, vector<double>()));
        slope.assign(dim, vector<vector<double>>(dim, vector<double>()));
        slope2.assign(dim, vector<vector<double>>(dim, vector<double>()));
        mean.assign(dim, vector<double>(dim,0));
        meaninv.assign(dim, vector<double>(dim,0));
        meanslope.assign(dim, vector<double>(dim,0));
        meanslope2.assign(dim, vector<double>(dim,0));
        varslope.assign(dim, vector<double>(dim,0));
        varslope2.assign(dim, vector<double>(dim,0));
        var.assign(dim, vector<double>(dim,0));
        varinv.assign(dim, vector<double>(dim,0));
        pp.assign(dim, vector<double>(dim,0));
        ppinv.assign(dim, vector<double>(dim,0));
        correl.assign(dim, vector<double>(dim,0));
        correl2.assign(dim, vector<double>(dim,0));
        partialcorrel.assign(dim, vector<double>(dim,0));
        partialcorrel2.assign(dim, vector<double>(dim,0));
		meanpropvar.assign(dim,0);
	}

	virtual ~MeanCovMatrix()	{
	}

	void SetLatex(bool b)	{
		tex = b;
	}

	int GetDim() const {
		return dim;
	}

	unsigned int GetSize() const	{
		return size;
	}

	double GetVal(int i, int j, int k) const	{
		if (val[i][j].size() != GetSize())	{
			cerr << "error : " << i << '\t' << j << '\t' << k << '\n';
			exit(1);
		}
		return val[i][j][k];
	}

	double GetInvVal(int i, int j, int k) const	{
		if (invval[i][j].size() != GetSize())	{
			cerr << "error : " << i << '\t' << j << '\t' << k << '\n';
			exit(1);
		}
		return invval[i][j][k];
	}

	virtual void Add(const CovMatrix& sample)	{
		sample.Diagonalise();
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				val[i][j].push_back(sample(i,j));
				invval[i][j].push_back(sample.GetInvMatrix()[i][j]);
			}
		}
		size++;
	}

	void ComputeSlopes()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					meanslope[i][j] = 0;
					varslope[i][j] = 0;
					for (unsigned int k=0; k<size; k++)	{
						double sxy = val[i][j][k];
						double sx = val[i][i][k];
						double sy = val[j][j][k];
						// double tmp = (sy + sqrt((sy - sx)*(sy - sx) - 4 *sxy *sxy - 2 * sx * sy)) / 2 / sxy;
						double tmp = (sy - sx + sqrt((sy - sx)*(sy - sx) + 4 *sxy *sxy)) / 2 / sxy;
						slope[i][j].push_back(tmp);
						meanslope[i][j] += tmp;
						varslope[i][j] += tmp * tmp;
					}
					meanslope[i][j] /= size;
					varslope[i][j] /= size;
					varslope[i][j] -= meanslope[i][j] * meanslope[i][j];
					sort(slope[i][j].begin(), slope[i][j].end());
				}
			}
		}
	}
	
	
	void ComputeSlopes2()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					meanslope2[i][j] = 0;
					varslope2[i][j] = 0;
					for (unsigned int k=0; k<size; k++)	{
						double sxy = val[i][j][k];
						double sx = val[i][i][k];
						// double sy = val[j][j][k];
						// double tmp = (sy + sqrt((sy - sx)*(sy - sx) - 4 *sxy *sxy - 2 * sx * sy)) / 2 / sxy;
						double tmp = sxy/sx;
						slope2[i][j].push_back(tmp);
						meanslope2[i][j] += tmp;
						varslope2[i][j] += tmp * tmp;
					}
					meanslope2[i][j] /= size;
					varslope2[i][j] /= size;
					varslope2[i][j] -= meanslope2[i][j] * meanslope2[i][j];
					sort(slope2[i][j].begin(), slope2[i][j].end());
				}
			}
		}
	}

	void Normalize()	{
		for (int i=0; i<GetDim(); i++)	{
			meanpropvar[i] = 0;
			for (int j=0; j<GetDim(); j++)	{
				mean[i][j] = 0;
				meaninv[i][j] = 0;
				var[i][j] = 0;
				varinv[i][j] = 0;
				pp[i][j] = 0;
				ppinv[i][j] = 0;
				correl[i][j] = 0;
				correl2[i][j] = 0;
				partialcorrel[i][j] = 0;
				partialcorrel2[i][j] = 0;
				if (val[i][j].size() != size)	{
					cerr << "error in meancovmatrix: non matching length\n";
					exit(1);
				}
				for (unsigned int k=0; k<size; k++)	{
					double& tmp = val[i][j][k];
					mean[i][j] += tmp;
					var[i][j] += tmp * tmp;
					if (tmp > 0)	{
						pp[i][j]++;
					}
					double r = tmp / sqrt(val[i][i][k] * val[j][j][k]);
					correl[i][j] += r;
					correl2[i][j] += r*r;
					double& inv = invval[i][j][k];
					meaninv[i][j] += inv;
					varinv[i][j] += inv * inv;
					if (inv < 0)	{
						ppinv[i][j] ++;
					}
					double pr = - inv / sqrt(invval[i][i][k] * invval[j][j][k]);
					partialcorrel[i][j] += pr;
					partialcorrel2[i][j] += pr * pr;

					meanpropvar[i] += 1.0 - 1.0 /(val[i][i][k] * invval[i][i][k]);
				}
				mean[i][j] /= size;
				var[i][j] /= size;
				var[i][j] -= mean[i][j] * mean[i][j];
				pp[i][j] /= size;
				meaninv[i][j] /= size;
				varinv[i][j] /= size;
				varinv[i][j] -= meaninv[i][j] * meaninv[i][j];
				ppinv[i][j] /= size;
				correl[i][j] /= size;
				correl2[i][j] /= size;
				partialcorrel[i][j] /= size;
				partialcorrel2[i][j] /= size;
				// sort(val[i][j].begin(), val[i][j].end());
				meanpropvar[i] /= size;
			}
		}
		ComputeSlopes();
		ComputeSlopes2();
	}

	double GetPropVariance(int i)	const {
		return meanpropvar[i];
		// return 1.0 - 1.0 / (mean[i][i] * meaninv[i][i]);
	}

	// MeanCovMatrix Project(bool* array);

	void PrintPropVariances(ostream& os) const {

		os << "proportion of variance of each trait explained by all other traits:\n";
		for (int i=0; i<GetDim(); i++)	{
			os << "trait " << i << " : " << GetPropVariance(i) << '\n';
		}
		os << '\n';
	}

	void PrintCovariances(ostream& os) const	{
		if (tex)	{
			os << "{\\sc covariance}";
			for (int i=0; i<GetDim(); i++)	{
				os << "&" << i;
			}
			os << "\\\\\n";
			os << "\\hline\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				os << i;
				for (int j=0; j<i; j++)	{
					os << "&-";
				}
				for (int j=i; j<GetDim(); j++)	{
					os << " & ";
					os << "$";
					os << ((double) ((int) (100 * mean[i][j]))) / 100;
					if ((pp[i][j] > 0.975) || (pp[i][j] < 0.025))	{
						os << "^*";
					}
					os << "$";
				}
				os << "\\\\\n";
			}
			os << "\\\\\n";
		}
		else	{
			os.precision(3);
			// output matrix of covariance
			os << "covariances\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					os << setw(7) << mean[i][j] << '\t';
				}
				os << '\n';
			}
			os << '\n';
		}
	}


	void PrintPrecisions(ostream& os) const	{
		os.precision(3);
		// output matrix of covariance
		os << "precisions\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << meaninv[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintR2(ostream& os) const	{
		if (tex)	{
			os << "{\\sc correlation ($R$)}";
			for (int i=0; i<GetDim(); i++)	{
				os << "&" << i;
			}
			os << "\\\\\n";
			os << "\\hline\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				os << i;
				for (int j=0; j<=i; j++)	{
					os << "&-";
				}
				for (int j=i+1; j<GetDim(); j++)	{
					os << " & ";
					os << "$";
					os << ((double) ((int) (100 * correl[i][j]))) / 100;
					// os << ((double) ((int) (100 * correl[i][j] * correl[i][j]))) / 100;
					if ((pp[i][j] > 0.975) || (pp[i][j] < 0.025))	{
						os << "^{**}";
					}
					else if ((pp[i][j] > 0.95) || (pp[i][j] < 0.05))	{
						os << "^*";
					}
					os << "$";
				}
				os << "\\\\\n";
			}
			os << "\\\\\n";
		}
		else	{

			os.precision(3);

			// correlation coefficients
			os << "r2\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					os << setw(7) << correl2[i][j] << '\t';
					// os << setw(7) << correl[i][j] * correl[i][j] << '\t';
				}
				os << '\n';
			}
			os << '\n';
		}

	}

	void PrintCorrel(ostream& os) const	{

		os.precision(3);

		// correlation coefficients
		os << "correlation coefficients\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << correl[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintPartialR2(ostream& os) const	{

		os.precision(3);

		// correlation coefficients
		os << "partial correlation coefficients\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << partialcorrel2[i][j] << '\t';
				// os << setw(7) << partialcorrel[i][j] * partialcorrel[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintPartialCorrel(ostream& os) const	{

		os.precision(3);

		// correlation coefficients
		os << "partial correlation coefficients\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << partialcorrel[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintPosteriorProbs(ostream& os) const	{

		if (tex)	{
			os << "{\\sc posterior prob.}";
			for (int i=0; i<GetDim(); i++)	{
				os << "&" << i;
			}
			os << "\\\\\n";
			os << "\\hline\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				os << i;
				for (int j=0; j<i; j++)	{
					os << "&-";
				}
				for (int j=i; j<GetDim(); j++)	{
					os << " & ";
					os << "$";
					os << ((double) ((int) (100 * pp[i][j]))) / 100;
					if ((pp[i][j] > 0.975) || (pp[i][j] < 0.025))	{
						os << "^*";
					}
					os << "$";
				}
				os << "\\\\\n";
			}
			os << "\\\\\n";
		}
		else	{
			os.precision(2);
			// pp
			os << "posterior probs\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					if (i != j)	{
						os << setw(3) << pp[i][j] << '\t';
					}
					else	{
						os << setw(3) << '-'  << '\t';
					}
				}
				os << '\n';
			}
			os << '\n';
			}
	}

	void PrintPrecisionsPosteriorProbs(ostream& os) const	{

		os.precision(2);

		// pp
		os << "posterior probs\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					os << setw(3) << ppinv[i][j] << '\t';
				}
				else	{
					os << setw(3) << '-'  << '\t';
				}
			}
			os << '\n';
		}
		os << '\n';
	}

    void PrintSlopesNe(ostream& os, int idxdNdS, int idxpiNpiS, int idxpiS) const   {

        int l = (int) (0.025 * size);
        int j = 0; 
        os << "relation\tmedian\tCI95min\tCI95max\n";
        os << "logdN/dS~logNe" << '\t' << meanslope[j][idxdNdS] << '\t' << slope[j][idxdNdS][l] << '\t' << slope[j][idxdNdS][size-1-l] << '\n';
        os << "logpiN/piS~logNe" << '\t' << meanslope[j][idxpiNpiS] << '\t' << slope[j][idxpiNpiS][l] << '\t' << slope[j][idxpiNpiS][size-1-l] << '\n';
        os << '\n';
        j = idxpiS;
        os << "logdN/dS~logpiS" << '\t' << meanslope[j][idxdNdS] << '\t' << slope[j][idxdNdS][l] << '\t' << slope[j][idxdNdS][size-1-l] << '\n';
        os << "logpiN/piS~logpiS" << '\t' << meanslope[j][idxpiNpiS] << '\t' << slope[j][idxpiNpiS][l] << '\t' << slope[j][idxpiNpiS][size-1-l] << '\n';
    }

	void PrintSlopes(ostream& os) const	{

		os.precision(3);

		// slopes
		os << "slopes \n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					int l = (int) (0.025 * size);
					// l = 0;
					// os << setw(7) << meanslope[j][i] << " ( " << setw(7) << slope[j][i][l] << " , " << setw(7) << slope[j][i][size-1-l] << " ) " << '\t';
					os << meanslope[j][i] << " : " << sqrt(varslope[j][i]) << " ( " << slope[j][i][l] << " , " << slope[j][i][size-1-l] << " ) " << '\t';
				}
				else	{
					os << '-' << '(' << '-' << " , " << '-' << " ) " << '\t';
				}

			}
			os << '\n';
			os << '\n';
		}
		os << '\n';
	}
	
	void PrintSlopes2(ostream& os) const	{

		os.precision(3);

		// slopes
		os << "slopes2 \n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					int l = (int) (0.025 * size);
					// l = 0;
					// os << setw(7) << meanslope[j][i] << " ( " << setw(7) << slope[j][i][l] << " , " << setw(7) << slope[j][i][size-1-l] << " ) " << '\t';
					os << meanslope2[j][i] << " : " << sqrt(varslope2[j][i]) << " ( " << slope2[j][i][l] << " , " << slope2[j][i][size-1-l] << " ) " << '\t';
				}
				else	{
					os << '-' << '(' << '-' << " , " << '-' << " ) " << '\t';
				}

			}
			os << '\n';
			os << '\n';
		}
		os << '\n';
	}

	void ToStream(ostream& os) const {
		if (tex)	{
			os << "\\begin{tabular}{";
			for (int i=0; i<GetDim() + 1; i++)	{
				os << "r";
			}
			os << "}\n";
			os << '\n';

			// PrintCovariances(os);
			PrintR2(os);
			// PrintPosteriorProbs(os);

			os << "\\end{tabular}\n";
		}
		else	{
			PrintCovariances(os);
			PrintCorrel(os);
			PrintPosteriorProbs(os);
			PrintPrecisions(os);
			PrintPartialCorrel(os);
			PrintPrecisionsPosteriorProbs(os);
			// PrintSlopes(os);
		}
	}

	friend ostream& operator<<(ostream& os, const MeanCovMatrix& m)	{
		m.ToStream(os);
		return os;
	}

	bool tex;
	int dim;
	unsigned int size;

	// arrays of values
	vector<vector<vector<double>>> val;
	vector<vector<vector<double>>> invval;
	vector<vector<vector<double>>> slope;
	vector<vector<vector<double>>> slope2;
	
	vector<vector<double>> mean;
	vector<vector<double>> meanslope;
	vector<vector<double>> meanslope2;
	vector<vector<double>> varslope;
	vector<vector<double>> varslope2;
	vector<vector<double>> var;
	vector<vector<double>> pp;
	vector<vector<double>> correl;
	vector<vector<double>> partialcorrel;
	vector<vector<double>> correl2;
	vector<vector<double>> partialcorrel2;
	vector<vector<double>> meaninv;
	vector<vector<double>> varinv;
	vector<vector<double>> ppinv;

	vector<double> meanpropvar;
};


