#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneSparseConditionOmegaModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneSparseConditionOmegaChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int ncond, nlevel;
    double pipos, pineg;
    int blmode, nucmode;

  public:
    MultiGeneSparseConditionOmegaModel *GetModel() {
        return static_cast<MultiGeneSparseConditionOmegaModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    MultiGeneSparseConditionOmegaChain(string indatafile, string intreefile, int inncond, int innlevel, double inpipos, double inpineg,
		    		 int inblmode, int innucmode,
                                 int inevery, int inuntil, string inname, int force, int inmyid,
                                 int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENESPARSECONDOMEGA"),
          datafile(indatafile),
          treefile(intreefile),
          ncond(inncond),
          nlevel(innlevel), pipos(inpipos), pineg(inpineg) {
	blmode = inblmode;
	nucmode = innucmode;
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    MultiGeneSparseConditionOmegaChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneSparseConditionOmegaModel(datafile, treefile, ncond, nlevel, pipos, pineg, myid, nprocs);
        GetModel()->SetAcrossGenesModes(blmode,nucmode);
        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        GetModel()->Update();
        Reset(force);

        if (!myid) {
            model->Trace(cerr);
        }
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        is >> ncond >> nlevel;
        is >> pipos >> pineg;
        is >> blmode >> nucmode;
        int tmp;
        is >> tmp;
        if (tmp) {
	    cerr << "error when reading model\n";
	    exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENESPARSECONDOMEGA") {
            model =
                new MultiGeneSparseConditionOmegaModel(datafile, treefile, ncond, nlevel, pipos, pineg, myid, nprocs);
                GetModel()->SetAcrossGenesModes(blmode,nucmode);
        } else {
            cerr << "error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        GetModel()->FromStream(is);
        GetModel()->Update();
        if (!myid) {
            cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
            model->Trace(cerr);
        }
    }

    void Save() override {
        if (!myid) {
            ofstream param_os((name + ".param").c_str());
            param_os << GetModelType() << '\n';
            param_os << datafile << '\t' << treefile << '\n';
            param_os << ncond << '\t' << nlevel << '\n';
            param_os << pipos << '\t' << pineg << '\n';
	    param_os << blmode << '\t' << nucmode << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) override {
        MultiGeneChain::MakeFiles(force);
        ofstream gos((name + ".gene").c_str());
        ofstream bos((name + ".cond").c_str());
        ofstream bgos((name + ".condgene").c_str());
        ofstream nameos((name + ".genelist").c_str());
        GetModel()->PrintGeneList(nameos);
        nameos.close();
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (!myid) {
            ofstream gos((name + ".gene").c_str(), ios_base::app);
            GetModel()->PrintGeneEffects(gos);
            ofstream bos((name + ".cond").c_str(), ios_base::app);
            GetModel()->PrintCondEffects(bos);
            ofstream bgos((name + ".condgene").c_str(), ios_base::app);
            GetModel()->PrintDeviations(bgos);
        }
    }
};

int main(int argc, char *argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (nprocs <= 1) {
        cerr << "error: should run the program with at least 2 cores\n";
        exit(1);
    }

    MultiGeneSparseConditionOmegaChain *chain = 0;
    string name = "";

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneSparseConditionOmegaChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int ncond = 0;
        int nlevel = 1;
        int force = 1;
        int every = 1;
        int until = -1;
        int blmode = 1;
        int nucmode = 1;
        double pipos = 0.1;
        double pineg = 0.1;

        try {
            if (argc == 1) {
                throw(0);
            }

            int i = 1;
            while (i < argc) {
                string s = argv[i];

                if (s == "-d") {
                    i++;
                    datafile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-ncond") {
                    i++;
                    ncond = atoi(argv[i]);
                } else if (s == "-pi")  {
                    i++;
                    pipos = atof(argv[i]);
                    i++;
                    pineg = atof(argv[i]);
                } else if (s == "-nodev")   {
                    pipos = pineg = 0;
                } else if (s == "-nucrates") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shared") {
                        nucmode = 2;
                    } else if (tmp == "shrunken") {
                        nucmode = 1;
                    } else if ((tmp == "ind") || (tmp == "independent")) {
                        nucmode = 0;
                    } else {
                        cerr << "error: does not recongnize command after -nucrates\n";
                        exit(1);
                    }
                } else if (s == "-bl") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shared") {
                        blmode = 2;
                    } else if (tmp == "shrunken") {
                        blmode = 1;
                    } else if ((tmp == "ind") || (tmp == "independent")) {
                        blmode = 0;
                    } else {
                        cerr << "error: does not recongnize command after -bl\n";
                        exit(1);
                    }
                } else if ((s == "-x") || (s == "-extract")) {
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                } else {
                    if (i != (argc - 1)) {
                        throw(0);
                    }
                    name = argv[i];
                }
                i++;
            }
            if ((datafile == "") || (treefile == "") || (name == "")) {
                throw(0);
            }
        } catch (...) {
            cerr << "globom -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneSparseConditionOmegaChain(datafile, treefile, ncond, nlevel, pipos, pineg, blmode, nucmode, 
						every, until, name, force, myid, nprocs);
    }

    if (!myid) {
        cerr << "chain " << name << " started\n";
    }
    chain->Start();
    if (!myid) {
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }

    MPI_Finalize();
}
