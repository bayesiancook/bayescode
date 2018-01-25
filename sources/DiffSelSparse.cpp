#include <cmath>
#include <fstream>
#include "DiffSelSparseModel.hpp"
#include "Chain.hpp"
using namespace std;


class DiffSelSparseChain: public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int codonmodel, category, level, fixglob, fixvar;

  public:
    DiffSelSparseModel* GetModel() {
        return static_cast<DiffSelSparseModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    DiffSelSparseChain(string indata, string intree, int incategory,
                                              int inlevel, int inevery, int inuntil, 
					      int incodonmodel,
                                              string inname, int force)
        : modeltype("DIFFSELSPARSE"),
          datafile(indata),
          treefile(intree),
          codonmodel(incodonmodel),
          category(incategory),
          level(inlevel)	{
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    DiffSelSparseChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new DiffSelSparseModel(datafile, treefile, category, level, codonmodel);
        GetModel()->Allocate();
        GetModel()->Unfold(true);
        cerr << "-- Reset" << endl;
        Reset(force);
        cerr << "-- initial ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile >> category >> level;
        is >> codonmodel;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "DIFFSELSPARSE") {
            model = new DiffSelSparseModel(
                datafile, treefile, category, level, codonmodel);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        model->Update();
        GetModel()->Unfold(true);
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\t' << category << '\t' << level << '\n';
        param_os << codonmodel << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';

        model->ToStream(param_os);
    }
};


int main(int argc, char* argv[]) {

    cerr << "-- Parsing command line arguments\n";

    // this is an already existing chain on the disk; reopen and restart
    if (argc == 2 && argv[1][0] != '-') {
        string name = argv[1];
        cerr << "-- Trying to reopen existing chain named " << name << " on disk\n";
    }

    // this is a new chain
    else {

	string datafile = "";
	string treefile = "";
	int ncond = 2;
	int nlevel = 1;
	int codonmodel = 1;

	string name = "";
	int every = 1;
	int until = -1;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-d")	{
				i++;
				datafile = argv[i];
			}
			else if ((s == "-t") || (s == "-T"))	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-ncond")	{
				i++;
				ncond = atoi(argv[i]);
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
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
	}
	catch(...)	{
		cerr << "error in command\n";
		exit(1);
	}

	DiffSelSparseChain* chain = new DiffSelSparseChain(datafile,treefile,ncond,nlevel,every,until,codonmodel,name,true);
	chain->Start();
    }
}
