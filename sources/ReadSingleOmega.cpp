#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "SingleOmegaModel.hpp"
using namespace std;


class SingleOmegaSample : public Sample {

	private:
	string modeltype;
	string datafile;
	string treefile;

	public:

	string GetModelType() {return modeltype;}

	SingleOmegaModel* GetModel() {return (SingleOmegaModel*) model;}

	SingleOmegaSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> datafile >> treefile;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
		is >> chainevery >> chainuntil >> chainsize;

		// make a new model depending on the type obtained from the file
		if (modeltype == "SINGLEOMEGA")	{
			model = new SingleOmegaModel(datafile,treefile);
		}
		else	{
			cerr << "error when opening file "  << name << '\n';
			exit(1);
		}

        GetModel()->Allocate();

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);

		/*
		cerr << "UPDATE\n";
		model->Update();
		*/

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read()	{

        cerr << size << " points to read\n";

        double meanomega = 0;
        double varomega = 0;

        for (int i=0; i<size; i++)  {
            cerr << '.';
            GetNextPoint();
            double om = GetModel()->GetOmega();
            meanomega += om;
            varomega += om*om;
        }
        cerr << '\n';
        meanomega /= size;
        varomega /= size;
        varomega -= meanomega*meanomega;

        cout << "posterior mean omega : " << meanomega << '\t' << sqrt(varomega) << '\n';
    }
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

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

	SingleOmegaSample* sample = new SingleOmegaSample(name,burnin,every,until);
    sample->Read();
}

