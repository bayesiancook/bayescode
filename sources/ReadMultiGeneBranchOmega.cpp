#include <cmath>
#include <fstream>
#include "MultiGeneBranchOmegaModel.hpp"
#include "MultiGeneSample.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneBranchOmegaSample : public MultiGeneSample {
  private:
    string modeltype;
    string datafile;
    string treefile;

  public:
    string GetModelType() { return modeltype; }

    const MultiGeneBranchOmegaModel *GetModel() const { return (MultiGeneBranchOmegaModel *)model; }
    MultiGeneBranchOmegaModel *GetModel() { return (MultiGeneBranchOmegaModel *)model; }

    MultiGeneBranchOmegaSample(
        string filename, int inburnin, int inevery, int inuntil, int myid, int nprocs)
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
        if (modeltype == "MULTIGENEBRANCHOMEGA") {
            model = new MultiGeneBranchOmegaModel(datafile, treefile, myid, nprocs);
        } else {
            cerr << "error when opening file " << name << '\n';
            cerr << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        // read model (i.e. chain's last point) from <name>.param
        if (!myid) {
            model->FromStream(is);

            // open <name>.chain, and prepare stream and stream iterator
            OpenChainFile();
            // now, size is defined (it is the total number of points with which this
            // Sample object will make all its various posterior averages) all these
            // points can be accessed to (only once) by repeated calls to
            // GetNextPoint()
        }
    }

    int GetNbranch() const { return GetModel()->GetNbranch(); }

    // a very simple (and quite uninteresting) method for obtaining
    // the posterior mean and variance of the total length of the tree
    void Read() {
        cerr << size << " points to read\n";

        vector<vector<double>> pp(GetNgene(), vector<double>(GetNbranch(), 0));

        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->FastUpdate();
            for (int gene = 0; gene < GetNgene(); gene++) {
                for (int j = 0; j < GetNbranch(); j++) {
                    if (GetModel()->GetOmega(gene, j) > GetModel()->GetMeanOmega(gene, j)) {
                        pp[gene][j]++;
                    }
                }
            }
        }
        cerr << '\n';
        for (int gene = 0; gene < GetNgene(); gene++) {
            for (int j = 0; j < GetNbranch(); j++) { pp[gene][j] /= size; }
        }

        ofstream os((name + ".pp").c_str());
        for (int j = 0; j < GetNbranch(); j++) {
            os << j << '\t';
            for (int gene = 0; gene < GetNgene(); gene++) { os << pp[gene][j] << '\t'; }
            os << '\n';
        }
        cerr << "pp of positive departure in " << name << ".pp\n";
    }
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

    try {
        if (argc == 1) { throw(0); }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if ((s == "-x") || (s == "-extract")) {
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                if (!IsInt(s)) { throw(0); }
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
                if (i != (argc - 1)) { throw(0); }
                name = argv[i];
            }
            i++;
        }
        if (name == "") { throw(0); }
    } catch (...) {
        cerr << "readglobom [-x <burnin> <every> <until>] <chainname> \n";
        cerr << '\n';
        exit(1);
    }

    MultiGeneBranchOmegaSample *sample =
        new MultiGeneBranchOmegaSample(name, burnin, every, until, myid, nprocs);
    if (!myid) { sample->Read(); }

    MPI_Finalize();
}
