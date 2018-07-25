#include <cmath>
#include <fstream>
#include "MultiGeneSample.hpp"
#include "MultiGeneSingleOmegaModel.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneSingleOmegaSample : public MultiGeneSample {
  private:
    string modeltype;
    string datafile;
    string treefile;

  public:
    string GetModelType() { return modeltype; }

    const MultiGeneSingleOmegaModel *GetModel() const { return (MultiGeneSingleOmegaModel *)model; }
    MultiGeneSingleOmegaModel *GetModel() { return (MultiGeneSingleOmegaModel *)model; }

    MultiGeneSingleOmegaSample(string filename, int inburnin, int inevery, int inuntil, int myid,
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
        if (modeltype == "MULTIGENESINGLEOMEGA") {
            model = new MultiGeneSingleOmegaModel(datafile, treefile, myid, nprocs);
        } else {
            cerr << "error when opening file " << name << '\n';
            cerr << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        // read model (i.e. chain's last point) from <name>.param
        model->FromStream(is);
        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    void MasterRead() {
        cerr << size << " points to read\n";

        vector<double> meanom(GetNgene(), 0);
        vector<double> varom(GetNgene(), 0);
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            const vector<double> &om = GetModel()->GetOmegaArray();
            for (int gene = 0; gene < GetNgene(); gene++) {
                meanom[gene] += om[gene];
                varom[gene] += om[gene] * om[gene];
            }
        }
        cerr << '\n';
        for (int gene = 0; gene < GetNgene(); gene++) {
            meanom[gene] /= size;
            varom[gene] /= size;
            varom[gene] -= meanom[gene] * meanom[gene];
        }

        ofstream os((name + ".postmeanomega").c_str());
        for (int gene = 0; gene < GetNgene(); gene++) {
            os << GetModel()->GetLocalGeneName(gene) << '\t' << meanom[gene] << '\t'
               << sqrt(varom[gene]) << '\n';
        }
        cerr << "posterior mean omega per gene in " << name << ".postmeanomega\n";
    }

    void SlaveRead() {
        for (int i = 0; i < size; i++) {
            GetNextPoint();
        }
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

    MultiGeneSingleOmegaSample *sample =
        new MultiGeneSingleOmegaSample(name, burnin, every, until, myid, nprocs);

    if (ppred) {
        if (!myid) {
            sample->MasterPostPred();
        } else {
            sample->SlavePostPred();
        }
    } else {
        if (!myid) {
            sample->MasterRead();
        } else {
            sample->SlaveRead();
        }
    }

    MPI_Finalize();
}
