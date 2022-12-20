
#include <cmath>
#include <fstream>
#include "MultiGeneSample.hpp"
#include "MultiGeneSingleOmegaModel.hpp"
#include "BranchToNewick.hpp"
using namespace std;

MPI_Datatype Propagate_arg;

class MultiGeneSingleOmegaSample : public MultiGeneSample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int blmode, nucmode, omegamode;
    double omegahypermean, omegahyperinvshape;

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
        is >> blmode >> nucmode >> omegamode;
        is >> omegahypermean >> omegahyperinvshape;
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
            GetModel()->SetAcrossGenesModes(blmode,nucmode,omegamode);
            GetModel()->SetOmegaHyperParameters(omegahypermean,omegahyperinvshape);
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

        SimpleBranchArray<double> meanbl(GetModel()->GetTree(),0);
        vector<double> meanom(GetNgene(), 0);
        vector<double> varom(GetNgene(), 0);
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->AddLength(meanbl);
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

        for (int j=0; j<GetModel()->GetNbranch(); j++)  {
            meanbl[j] /= size;
        }
        ofstream blos((name + ".postmeansynrate.tre").c_str());
        BranchToNewick(blos, meanbl);
        cerr << "posterior mean branch lengths in " << name << ".postmeansynrate.tre\n";


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

    void MasterReaddSOmegaPathSuffStat() {
        vector<dSOmegaPathSuffStatBranchArray> array(GetModel()->GetNgene(), dSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        vector<GCConsdSOmegaPathSuffStatBranchArray> gcconsarray(GetModel()->GetNgene(), GCConsdSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->MasterUpdate();
        }
        cerr << '\n';

        GetModel()->MasterReceiveGeneArray(array);

        ofstream dsom_os((name + ".genedsomsuffstat").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            dsom_os << GetModel()->GetLocalGeneName(gene) << '\n';
            dsom_os << "counts_dS\t";
            array[gene].BranchToNewickSynCount(dsom_os);
            dsom_os << "counts_dS_norm\t";
            array[gene].BranchToNewickSynBeta(dsom_os);
            dsom_os << "counts_dN\t";
            array[gene].BranchToNewickNonSynCount(dsom_os);
            dsom_os << "counts_dN_norm\t";
            array[gene].BranchToNewickNonSynBeta(dsom_os);
        }

        /*
        ofstream syn_count_os((name + ".genewise.counts_dS.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            syn_count_os << GetModel()->GetLocalGeneName(gene) << '\t';
            array[gene].BranchToNewickSynCount(syn_count_os);
        }

        ofstream syn_beta_os((name + ".genewise.counts_dS_norm.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            syn_beta_os << GetModel()->GetLocalGeneName(gene) << '\t';
            array[gene].BranchToNewickSynBeta(syn_beta_os);
        }

        ofstream nonsyn_count_os((name + ".genewise.counts_dN.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            nonsyn_count_os << GetModel()->GetLocalGeneName(gene) << '\t';
            array[gene].BranchToNewickNonSynCount(nonsyn_count_os);
        }

        ofstream nonsyn_beta_os((name + ".genewise.counts_dN_norm.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            nonsyn_beta_os << GetModel()->GetLocalGeneName(gene) << '\t';
            array[gene].BranchToNewickNonSynBeta(nonsyn_beta_os);
        }
        */

        GetModel()->MasterReceiveGeneArray(gcconsarray);

        /*
        ofstream gccons_syn_count_os((name + ".gccons.genewise.counts_dS.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            gccons_syn_count_os << GetModel()->GetLocalGeneName(gene) << '\t';
            gcconsarray[gene].BranchToNewickSynCount(gccons_syn_count_os);
        }

        ofstream gccons_syn_beta_os((name + ".gccons.genewise.counts_dS_norm.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            gccons_syn_beta_os << GetModel()->GetLocalGeneName(gene) << '\t';
            gcconsarray[gene].BranchToNewickSynBeta(gccons_syn_beta_os);
        }

        ofstream gccons_nonsyn_count_os((name + ".gccons.genewise.counts_dN.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            gccons_nonsyn_count_os << GetModel()->GetLocalGeneName(gene) << '\t';
            gcconsarray[gene].BranchToNewickNonSynCount(gccons_nonsyn_count_os);
        }

        ofstream gccons_nonsyn_beta_os((name + ".gccons.genewise.counts_dN_norm.dnd").c_str());
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            gccons_nonsyn_beta_os << GetModel()->GetLocalGeneName(gene) << '\t';
            gcconsarray[gene].BranchToNewickNonSynBeta(gccons_nonsyn_beta_os);
        }
        */
    }

    void SlaveReaddSOmegaPathSuffStat() {
        vector<dSOmegaPathSuffStatBranchArray> array(GetModel()->GetLocalNgene(), dSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        vector<GCConsdSOmegaPathSuffStatBranchArray> gcconsarray(GetModel()->GetLocalNgene(), GCConsdSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        for (int i = 0; i < size; i++) {
            GetNextPoint();
            GetModel()->SlaveUpdate();
            GetModel()->SlaveAdddSOmegaPathSuffStat(array);
            GetModel()->SlaveAddGCConsdSOmegaPathSuffStat(gcconsarray);
        }
        for (int i=0; i<GetModel()->GetLocalNgene(); i++) {
            array[i].Normalize(1.0/size);
            gcconsarray[i].Normalize(1.0/size);
        }
        GetModel()->SlaveSendGeneArray(array);
        GetModel()->SlaveSendGeneArray(gcconsarray);
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
    int dsom = 0;

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
            } else if (s == "-dsomss")    {
                dsom = 1;
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
    } else if (dsom)    {
        if (! myid) {
            sample->MasterReaddSOmegaPathSuffStat();
        }
        else    {
            sample->SlaveReaddSOmegaPathSuffStat();
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
