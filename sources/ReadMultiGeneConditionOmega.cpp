#include <cmath>
#include <fstream>
#include "MultiGeneSample.hpp"
#include "MultiGeneConditionOmegaModel.hpp"
#include "RecursiveNewick.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneConditionOmegaSample : public MultiGeneSample {
  private:
    string modeltype, datafile, treefile;
    int ncond, nlevel;
    int blmode, nucmode, devmode;

  public:
    string GetModelType() { return modeltype; }

    const MultiGeneConditionOmegaModel *GetModel() const { return (MultiGeneConditionOmegaModel *)model; }
    MultiGeneConditionOmegaModel *GetModel() { return (MultiGeneConditionOmegaModel *)model; }

    MultiGeneConditionOmegaSample(string filename, int inburnin, int inevery, int inuntil, int myid,
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
        is >> blmode >> nucmode >> devmode;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "MULTIGENECONDOMEGA") {
                model = new MultiGeneConditionOmegaModel(datafile, treefile, ncond, nlevel, myid, nprocs);
                GetModel()->SetAcrossGenesModes(blmode,nucmode);
                GetModel()->SetDeviationMode(devmode);
        } else {
            cerr << "error when opening file " << name << '\n';
            cerr << modeltype << '\n';
            exit(1);
        }

        GetModel()->Allocate();
        model->FromStream(is);
        OpenChainFile();
    }

    void MasterRead(double cutoff) {
        cerr << size << " points to read\n";
        vector<vector<double>> pp(GetModel()->GetNcond(),vector<double>(GetModel()->GetNgene(),0));
        vector<vector<double>> effect(GetModel()->GetNcond(),vector<double>(GetModel()->GetNgene(),0));

        vector<double> condmean(GetModel()->GetNcond(), 0);

        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->FastUpdate();
            double meangene = GetModel()->GetMeanGeneEffect();
            for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
                double tmp = GetModel()->GetCondEffect(cond);
                condmean[cond] += tmp * meangene;
                for (int gene=0; gene<GetModel()->GetNgene(); gene++)    {
                    double temp = log(GetModel()->GetOmega(gene,cond) / GetModel()->GetMeanOmega(gene,cond));
                    if (temp > 0)   {
                        pp[cond][gene]++;
                    }
                    effect[cond][gene] += temp;
                }
            }
        }
        cerr << '\n';

        for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
            condmean[cond] /= size;
        }

        for (int gene=0; gene<GetModel()->GetNgene(); gene++)    {
            for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
                effect[cond][gene] /= size;
                pp[cond][gene] /= size;
            }
        }

        ofstream cos((name + ".postmeancondeffects").c_str());
        if (GetModel()->PerBranch())    {
            Tabulate(cos, GetModel()->GetTree(), condmean, false);
        }
        else    {
            for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
                cos << cond << '\t' << condmean[cond] << '\n';
            }
        }

        if (GetModel()->PerBranch())    {
            ofstream cos((name + ".postmeanleafbrancheffects").c_str());
            Tabulate(cos, GetModel()->GetTree(), condmean, true);
        }

        ofstream os((name + ".postmeaneffects").c_str());
        os << "gene";
        for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
            os << '\t' << "cond" << cond << '\t' << "pp";
        }
        os << '\n';
        for (int gene=0; gene<GetModel()->GetNgene(); gene++)    {
            int report = 0;
            for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
                if ((pp[cond][gene] > cutoff) || (pp[cond][gene] < (1-cutoff)))    {
                    report = 1;
                }
            }
            if (report) {
                os << GetModel()->GetLocalGeneName(gene);
                for (int cond=0; cond<GetModel()->GetNcond(); cond++)    {
                    os << '\t' << effect[cond][gene] << '\t' << pp[cond][gene];
                }
                os << '\n';
            }
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

    int burnin = 0;
    int every = 1;
    int until = -1;
    string name;
    int ppred = 0;

    double cutoff = 0.7;

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
            } else if (s == "-c")   {
                cutoff = atof(argv[i]);
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

    MultiGeneConditionOmegaSample *sample =
        new MultiGeneConditionOmegaSample(name, burnin, every, until, myid, nprocs);

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
            sample->MasterRead(cutoff);
        } else {
            sample->SlaveRead();
        }
    }

    MPI_Finalize();
}
