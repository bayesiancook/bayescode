
#include <cmath>
#include <fstream>
#include "MultiGeneSample.hpp"
#include "MultiGeneGTRModel.hpp"
#include "BranchToNewick.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneGTRSample : public MultiGeneSample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int blmode, nucmode;

  public:
    string GetModelType() { return modeltype; }

    const MultiGeneGTRModel *GetModel() const { return (MultiGeneGTRModel *)model; }
    MultiGeneGTRModel *GetModel() { return (MultiGeneGTRModel *)model; }

    MultiGeneGTRSample(string filename, int inburnin, int inevery, int inuntil, int myid,
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
        is >> blmode >> nucmode;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "MULTIGENEGTR") {
            model = new MultiGeneGTRModel(datafile, treefile, myid, nprocs);
            GetModel()->SetAcrossGenesModes(blmode,nucmode);
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
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
        }
        cerr << '\n';
    }

    void SlaveRead() {
        for (int i = 0; i < size; i++) {
            GetNextPoint();
        }
    }

    void MasterReadGeneNodePathSuffStat() {
        vector<RelativePathSuffStatNodeArray> array(GetModel()->GetNgene(), RelativePathSuffStatNodeArray(GetModel()->GetTree(), GetModel()->GetStateSpace()->GetNstate()));
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->MasterUpdate();
        }
        cerr << '\n';

        GetModel()->MasterReceiveGeneArray(array);

        RelativePathSuffStatNodeArray global_pathss(GetModel()->GetTree(), GetModel()->GetStateSpace()->GetNstate());
        global_pathss.Clear();
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            global_pathss.Add(array[gene]);
        }
        ofstream gos((name + ".globalnodepathsuffstat").c_str());
        gos << global_pathss << '\n';
        /*
        for (int j=0; j<GetModel()->GetTree().GetNnode(); j++)   {
            gos << global_pathss[j] << '\n';
        }
        */

        ofstream os((name + ".genenodepathsuffstat").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            os << GetModel()->GetLocalGeneName(gene) << '\t';
            os << array[gene];
            os << '\n';
        }
    }

    void SlaveReadGeneNodePathSuffStat() {
        vector<RelativePathSuffStatNodeArray> array(GetModel()->GetNgene(), RelativePathSuffStatNodeArray(GetModel()->GetTree(), GetModel()->GetStateSpace()->GetNstate()));
        for (int i = 0; i < size; i++) {
            GetNextPoint();
            GetModel()->SlaveUpdate();
            GetModel()->SlaveAddGeneNodePathSuffStat(array);
        }
        for (int i=0; i<GetModel()->GetLocalNgene(); i++) {
            array[i].Normalize(1.0/size);
        }
        GetModel()->SlaveSendGeneArray(array);
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
    int nodepathss = 0;

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
            } else if (s == "-nodepathss")	{
                nodepathss = 1;
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

    MultiGeneGTRSample *sample =
        new MultiGeneGTRSample(name, burnin, every, until, myid, nprocs);

    if (ppred) {
        if (!myid) {
            sample->MasterPostPred();
        } else {
            sample->SlavePostPred();
        }
    } else if (nodepathss)    {
        if (! myid) {
            sample->MasterReadGeneNodePathSuffStat();
        }
        else    {
            sample->SlaveReadGeneNodePathSuffStat();
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
