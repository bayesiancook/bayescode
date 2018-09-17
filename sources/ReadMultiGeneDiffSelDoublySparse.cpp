#include <cmath>
#include <fstream>
#include "MultiGeneDiffSelDoublySparseModel.hpp"
#include "MultiGeneSample.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneDiffSelDoublySparseSample : public MultiGeneSample {
  private:
    string modeltype, datafile, treefile;
    int ncond;
    int nlevel;
    int codonmodel;
    double epsilon;
    double fitnessshape;
    int fitnesscentermode;
    int blmode, nucmode, shiftmode;
    double pihypermean, pihyperinvconc;
    double shiftprobmean, shiftprobinvconc;
    int burnin;
    int writegenedata;

  public:
    string GetModelType() { return modeltype; }

    const MultiGeneDiffSelDoublySparseModel *GetModel() const {
        return (MultiGeneDiffSelDoublySparseModel *)model;
    }
    MultiGeneDiffSelDoublySparseModel *GetModel() {
        return (MultiGeneDiffSelDoublySparseModel *)model;
    }

    MultiGeneDiffSelDoublySparseSample(string filename, int inburnin, int inevery, int inuntil,
                                       int myid, int nprocs)
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
        is >> ncond >> nlevel >> codonmodel;
        is >> epsilon >> fitnessshape;
        is >> fitnesscentermode;
        is >> blmode >> nucmode >> shiftmode;
        is >> pihypermean >> pihyperinvconc;
        is >> shiftprobmean >> shiftprobinvconc;

        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> burnin;
        is >> chainevery >> chainuntil >> chainsaveall >> writegenedata >> chainsize;

        if (modeltype == "MULTIGENEDIFFSELDSPARSE") {
            model = new MultiGeneDiffSelDoublySparseModel(
                datafile, treefile, ncond, nlevel, codonmodel, epsilon, fitnessshape, blmode,
                nucmode, shiftmode, pihypermean, pihyperinvconc, shiftprobmean, shiftprobinvconc,
                myid, nprocs);
        } else {
            cerr << "Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }

        if (chainsize < burnin) {
            GetModel()->SetWithToggles(0);
        } else {
            GetModel()->SetWithToggles(1);
        }
        GetModel()->SetFitnessCenterMode(fitnesscentermode);

        GetModel()->Allocate();
        GetModel()->FromStream(is);

        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    /*
    void MasterRead() {
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            MPI_Barrier(MPI_COMM_WORLD);
        }
        cerr << '\n';
    }

    void SlaveRead() {
        for (int i = 0; i < size; i++) {
            GetNextPoint();
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    */
};

int main(int argc, char *argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int blockcounts[2] = {1, 3};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
    MPI_Aint dtex, displacements[2];

    displacements[0] = (MPI_Aint)0;
    MPI_Type_extent(MPI_DOUBLE, &dtex);
    displacements[1] = dtex;
    MPI_Type_struct(2, blockcounts, displacements, types, &Propagate_arg);
    MPI_Type_commit(&Propagate_arg);

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

    MultiGeneDiffSelDoublySparseSample *sample =
        new MultiGeneDiffSelDoublySparseSample(name, burnin, every, until, myid, nprocs);

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
