#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "DiffSelSparseModel.hpp"
using namespace std;

/**
 * \brief An MCMC sample for DiffSelSparseModel
 *
 * implements a read function, returning the MCMC estimate of the posterior probabilities of differential effects across sites
 */


class DiffSelSparseSample : public Sample {

	private:
	string modeltype;
	string datafile;
	string treefile;
    int ncond, nlevel, codonmodel, fixhyper;

	public:

	string GetModelType() override {return modeltype;}

	DiffSelSparseModel* GetModel() {return (DiffSelSparseModel*) model;}

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see Sample)
	DiffSelSparseSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	override {

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
        is >> datafile >> treefile >> ncond >> nlevel;
        is >> codonmodel;
        is >> fixhyper;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
		is >> chainevery >> chainuntil >> chainsaveall >> chainsize;

		// make a new model depending on the type obtained from the file
		if (modeltype == "DIFFSELSPARSE")	{
            model = new DiffSelSparseModel(datafile, treefile, ncond, nlevel, codonmodel);
            GetModel()->SetFitnessHyperMode(fixhyper);
		}
		else	{
			cerr << "error when opening file "  << name << '\n';
			exit(1);
		}

        GetModel()->Allocate();

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	//! \brief computes the posterior mean estimate (and the posterior standard deviation) of omega
	void ReadPP(double cutoff, int siteoffset)	{

        int Nsite = GetModel()->GetNsite();
        int Ncond = ncond;
        vector<vector<vector<double> > > pp(Ncond-1, vector<vector<double> > (Nsite, vector<double>(Naa,0)));
        vector<vector<double> > sitepp(Ncond-1, vector<double>(Nsite,0));
        for (int i=0; i<size; i++)  {
            cerr << '.';
            GetNextPoint();
            for (int k=1; k<Ncond; k++) {
                const vector<vector<int> > & toggle = GetModel()->GetCondToggleArray(k);
                for (int j=0; j<Nsite; j++) {
                    int n = 0;
                    for (int a=0; a<Naa; a++)   {
                        pp[k-1][j][a] += toggle[j][a];
                        n += toggle[j][a];
                    }
                    if (n)  {
                        sitepp[k-1][j] ++;
                    }
                }
            }
        }
        cerr << '\n';

        // normalization
        for (int k=1; k<Ncond; k++) {
            for (int j=0; j<Nsite; j++) {
                for (int a=0; a<Naa; a++)   {
                    pp[k-1][j][a] /= size;
                }
                sitepp[k-1][j] /= size;
            }
        }

        // write output
        for (int k=1; k<Ncond; k++) {
            ostringstream s;
            s << name << "_" << k << ".sitepp";
            ofstream os(s.str().c_str());
            os << "site\tsitepp";
            for (int a=0; a<Naa; a++)   {
                os << '\t' << AminoAcids[a];
            }
            os << '\n';
            for (int j=0; j<Nsite; j++) {
                if (sitepp[k-1][j] > cutoff) {
                    os << j + siteoffset << '\t' << int(100*sitepp[k-1][j]);
                    for (int a=0; a<Naa; a++)   {
                        os << '\t' << int(100*pp[k-1][j][a]);
                    }
                    os << '\n';
                }
            }
        }

        // write output
        for (int k=1; k<Ncond; k++) {
            ostringstream s;
            s << name << "_" << k << ".pp";
            ofstream os(s.str().c_str());
            os << "site\tAA\tpp\n";
            for (int j=0; j<Nsite; j++) {
                for (int a=0; a<Naa; a++)   {
                    if (pp[k-1][j][a] > cutoff) {
                        os << j+siteoffset << '\t' << AminoAcids[a] << '\t' << int(100*pp[k-1][j][a]) << '\n';
                    }
                }
            }
        }

        cerr << "post probs written in .pp files\n";
        cerr << '\n';
    }
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

    int siteoffset = 0;
    double cutoff = 0.90;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				until = atoi(argv[i]);
			}
            else if (s == "-c") {
                i++;
                cutoff = atof(argv[i]);
            }
            else if (s == "-s") {
                i++;
                siteoffset = atoi(argv[i]);
            }
			else	{
					if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readglobom [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	DiffSelSparseSample* sample = new DiffSelSparseSample(name,burnin,every,until);
    sample->ReadPP(cutoff,siteoffset);
}


