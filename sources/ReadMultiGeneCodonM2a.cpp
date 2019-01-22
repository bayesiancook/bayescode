#include <cmath>
#include <fstream>
#include "MultiGeneCodonM2aModel.hpp"
#include "MultiGeneSample.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneCodonM2aSample : public MultiGeneSample {
  private:
    string modeltype, datafile, treefile;
    int writegenedata;
    int blmode, blsamplemode, nucmode, purommode, dposommode, purwmode, poswmode;
    double pihypermean, pihyperinvconc;
    double puromhypermean, puromhyperinvconc;
    double dposomhypermean, dposomhyperinvshape;
    double purwhypermean, purwhyperinvconc;
    double poswhypermean, poswhyperinvconc;
    int modalprior;

  public:
    string GetModelType() { return modeltype; }

    const MultiGeneCodonM2aModel *GetModel() const { return (MultiGeneCodonM2aModel *)model; }
    MultiGeneCodonM2aModel *GetModel() { return (MultiGeneCodonM2aModel *)model; }

    MultiGeneCodonM2aSample(string filename, int inburnin, int inevery, int inuntil, int myid,
                            int nprocs)
        : MultiGeneSample(filename, inburnin, inevery, inuntil, myid, nprocs) {
        Open();
    }

    void Open() {
        // open <name>.param
        ifstream is((name + ".param").c_str());

        // check that file exists
        if (!is) {
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }

        is >> modeltype;
        is >> datafile >> treefile;
        is >> writegenedata;
        is >> blmode >> blsamplemode >> nucmode >> dposommode >> purwmode >> poswmode;
        is >> pihypermean >> pihyperinvconc;
        is >> puromhypermean >> puromhyperinvconc;
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> purwhypermean >> purwhyperinvconc;
        is >> poswhypermean >> poswhyperinvconc;
        is >> modalprior;

        int check;
        is >> check;
        if (check)  {
            cerr << "Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        if (modeltype == "MULTIGENECODONM2A") {
            model = new MultiGeneCodonM2aModel(datafile, treefile, pihypermean, pihyperinvconc,
                                               myid, nprocs);
            GetModel()->SetAcrossGenesModes(blmode, nucmode, purommode, dposommode, purwmode,
                                            poswmode);
            GetModel()->SetBLSamplingMode(blsamplemode);
            GetModel()->SetMixtureHyperParameters(
                puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape,
                purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc);
            GetModel()->SetModalMixturePrior(modalprior);
        } else {
            cerr << "Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        GetModel()->FromStream(is);

        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    void MasterRead() {}

    void SlaveRead() {}
};

int main(int argc, char *argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int burnin = 0;
    int every = 1;
    int until = -1;
    string name;
    int ppred = 0;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if (s == "-ppred") {
                ppred = 1;
            } else if ((s == "-x") || (s == "-extract")) {
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                if (!IsInt(s)) {
                    throw(0);
                }
                burnin = atoi(argv[i]);
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                if (IsInt(s)) {
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    s = argv[i];
                    if (IsInt(s)) {
                        until = atoi(argv[i]);
                    } else {
                        i--;
                    }
                } else {
                    i--;
                }
            } else {
                if (i != (argc - 1)) {
                    throw(0);
                }
                name = argv[i];
            }
            i++;
        }
        if (name == "") {
            throw(0);
        }
    } catch (...) {
        cerr << "readglobom [-x <burnin> <every> <until>] <chainname> \n";
        cerr << '\n';
        exit(1);
    }

    MultiGeneCodonM2aSample *sample =
        new MultiGeneCodonM2aSample(name, burnin, every, until, myid, nprocs);

    if (ppred) {
        if (!myid) {
            sample->MasterPostPred();
        } else {
            sample->SlavePostPred();
        }
    } else {
        cerr << "default read function not yet implemented\n";
        /*
        if (! myid) {
            sample->MasterRead();
        }
        else    {
            sample->SlaveRead();
        }
        */
    }

    MPI_Finalize();
}
