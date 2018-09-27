#include <cmath>
#include <fstream>
#include "MultiGeneSample.hpp"
#include "MultiGeneSparseConditionOmegaModel.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneSparseConditionOmegaSample : public MultiGeneSample {
  private:
    string modeltype, datafile, treefile;
    int ncond, nlevel;
    double pipos, pineg;
    int blmode, nucmode;

  public:
    string GetModelType() { return modeltype; }

    const MultiGeneSparseConditionOmegaModel *GetModel() const { return (MultiGeneSparseConditionOmegaModel *)model; }
    MultiGeneSparseConditionOmegaModel *GetModel() { return (MultiGeneSparseConditionOmegaModel *)model; }

    MultiGeneSparseConditionOmegaSample(string filename, int inburnin, int inevery, int inuntil, int myid,
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
        is >> ncond >> nlevel;
        is >> pipos >> pineg;
        is >> blmode >> nucmode;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "MULTIGENESPARSECONDOMEGA") {
                model = new MultiGeneSparseConditionOmegaModel(datafile, treefile, ncond, nlevel, pipos, pineg, myid, nprocs);
                GetModel()->SetAcrossGenesModes(blmode,nucmode);
        } else {
            cerr << "error when opening file " << name << '\n';
            cerr << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        model->FromStream(is);
        OpenChainFile();
    }

    void MasterRead() {
        cerr << size << " points to read\n";
        vector<vector<double>> pppos(GetModel()->GetNcond(),vector<double>(GetModel()->GetNgene(),0));
        vector<vector<double>> ppneg(GetModel()->GetNcond(),vector<double>(GetModel()->GetNgene(),0));
        vector<vector<double>> effect(GetModel()->GetNcond(),vector<double>(GetModel()->GetNgene(),0));
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->FastUpdate();
            for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
                for (int gene=0; gene<GetModel()->GetNgene(); gene++)    {
                    int tmp = GetModel()->GetAlloc(gene,cond);
                    if (tmp == 0)   {
                        ppneg[cond][gene]++;
                    }
                    else if (tmp == 2)  {
                        pppos[cond][gene]++;
                    }
                    double temp = log(GetModel()->GetOmega(gene,cond)) - GetModel()->GetMeanLogOmega(gene,cond);
                    effect[cond][gene] += temp;
                }
            }
        }
        cerr << '\n';

        ofstream os((name + ".postmeaneffects").c_str());
        os << "gene";
        for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
            os << '\t' << "cond" << cond << '\t' << "pp+" << '\t' << "pp-";
        }
        os << '\n';
        for (int gene=0; gene<GetModel()->GetNgene(); gene++)    {
            os << GetModel()->GetLocalGeneName(gene);
            for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
                effect[cond][gene] /= size;
                pppos[cond][gene] /= size;
                ppneg[cond][gene] /= size;
                os << '\t' << effect[cond][gene] << '\t' << pppos[cond][gene] << '\t' << ppneg[cond][gene];
            }
            os << '\n';
        }
        cerr << "post mean log deviations (and post probs) in " << name << ".postmeaneffects\n";
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
            } else if (s == "-ppred0")    {
                ppred = 2;
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

    MultiGeneSparseConditionOmegaSample *sample =
        new MultiGeneSparseConditionOmegaSample(name, burnin, every, until, myid, nprocs);

    if (ppred) {
        if (ppred == 2) {
            sample->GetModel()->SetPostPredMode(0);
        }
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
