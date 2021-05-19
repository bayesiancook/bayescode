
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

    void MasterReadNucRates() {
        cerr << size << " points to read\n";

        vector<vector<vector<double>>> meannucrates(GetNgene(), vector<vector<double>>(Nnuc, vector<double>(Nnuc, 0)));
        vector<vector<double>> meanstat(GetNgene(), vector<double>(Nnuc, 0));
        vector<double> meanom(GetNgene(), 0);

        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->AddNucRates(meannucrates, meanstat);
            const vector<double> &om = GetModel()->GetOmegaArray();
            for (int gene = 0; gene < GetNgene(); gene++) {
                meanom[gene] += om[gene];
            }
        }
        cerr << '\n';
        for (int gene = 0; gene < GetNgene(); gene++) {
            meanom[gene] /= size;
            for (int i=0; i<Nnuc; i++)  {
                for (int j=0; j<Nnuc; j++)  {
                    meannucrates[gene][i][j] /= size;
                }
            }
            for (int i=0; i<Nnuc; i++)  {
                meanstat[gene][i] /= size;
            }
        }

        ofstream os((name + ".postmeangenenucrates").c_str());
        os << "#gene";
        for (int i=0; i<Nnuc; i++)  {
            for (int j=0; j<Nnuc; j++)  {
                if (i != j) {
                    os << '\t' << DNAletters[i] << DNAletters[j];
                }
            }
        }
        for (int i=0; i<Nnuc; i++)  {
            os << '\t' << "pi" << DNAletters[i];
        }
        os << "\t%GC";
        os << "\tomega";
        os << '\n';
        for (int gene = 0; gene < GetNgene(); gene++) {
            os << GetModel()->GetLocalGeneName(gene);
            for (int i=0; i<Nnuc; i++)  {
                for (int j=0; j<Nnuc; j++)  {
                    if (i != j) {
                        os << '\t' << meannucrates[gene][i][j];
                    }
                }
            }
            for (int i=0; i<Nnuc; i++)  {
                os << '\t' << meanstat[gene][i];
            }
            os << '\t' << meanstat[gene][1] + meanstat[gene][2];
            os << '\t' << meanom[gene];
            os << '\n';
        }
        cerr << "posterior mean nuc rates across genes in " << name << ".postmeangenenucrates\n";
    }

    void SlaveReadNucRates() {
        for (int i = 0; i < size; i++) {
            GetNextPoint();
        }
    }

    void WritedSGCTree(ostream& os, const BranchSelector<double>& dstree, const BranchSelector<double>& gctree) {
        RecursiveWritedSGCTree(GetModel()->GetTree().GetRoot(), os, dstree, gctree);
    }

    void RecursiveWritedSGCTree(const Link* from, ostream& os, const BranchSelector<double>& dstree, const BranchSelector<double>& gctree)  {

        if (from->isLeaf()) {
            os << from->GetNode()->GetName();
            os << "_";
        }
        else    {
            os << "(";
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                RecursiveWritedSGCTree(link->Out(), os, dstree, gctree);
                if (link->Next() != from)   {
                    os << ",";
                }
            }
            os << ")";
        }
        if (! from->isRoot())	{
            os << gctree.GetVal(from->GetBranch()->GetIndex());
            os << ":";
            os << dstree.GetVal(from->GetBranch()->GetIndex());
        }
        else    {
            os << gctree.GetVal(from->Next()->GetBranch()->GetIndex());
        }
    }

    void WriteGCLeafTab(ostream& os, const BranchSelector<double>& gctree) {
        RecursiveWriteGCLeafTab(GetModel()->GetTree().GetRoot(), os, gctree);
    }

    void RecursiveWriteGCLeafTab(const Link* from, ostream& os, const BranchSelector<double>& gctree)  {

        if (from->isLeaf()) {
            os << from->GetNode()->GetName() << '\t' << gctree.GetVal(from->GetBranch()->GetIndex()) << '\n';
        }
        else    {
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                RecursiveWriteGCLeafTab(link->Out(), os, gctree);
            }
        }
    }

    void MasterReaddSOmegaPathSuffStat() {
        vector<dSOmegaPathSuffStatBranchArray> array(GetModel()->GetNgene(), dSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        vector<GCConsdSOmegaPathSuffStatBranchArray> gcconsarray(GetModel()->GetNgene(), GCConsdSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        vector<GCCodonPathSuffStatBranchArray> gcarray(GetModel()->GetNgene(), GCCodonPathSuffStatBranchArray(GetModel()->GetTree()));
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->MasterUpdate();
        }
        cerr << '\n';

        GetModel()->MasterReceiveGeneArray(array);
        ofstream os((name + ".genebranchdsomsuffstat").c_str());
        for (int i=0; i<GetModel()->GetNgene(); i++) {
            os << array[i] << '\n';
        }
        cerr << "gene dsom path suffstats in " << name << ".genebranchdsomsuffstat\n";

        dSOmegaPathSuffStatBranchArray globdsomss(GetModel()->GetTree());
        for (int i=0; i<GetModel()->GetNgene(); i++) {
            globdsomss.Add(array[i]);
        }

        ofstream gos((name + ".meanbranchdsomsuffstat").c_str());
        gos << globdsomss << '\n';
        cerr << "global dsom path suffstats in " << name << ".meanbranchdsomsuffstat\n";

        GetModel()->MasterReceiveGeneArray(gcconsarray);
        ofstream gcos((name + ".genebranchgcconsdsomsuffstat").c_str());
        for (int i=0; i<GetModel()->GetNgene(); i++) {
            gcos << gcconsarray[i] << '\n';
        }
        cerr << "GC-cons gene dsom path suffstats in " << name << ".genebranchgcconsdsomsuffstat\n";
        GCConsdSOmegaPathSuffStatBranchArray gcconsglobdsomss(GetModel()->GetTree());
        for (int i=0; i<GetModel()->GetNgene(); i++) {
            gcconsglobdsomss.Add(gcconsarray[i]);
        }
        ofstream gcgos((name + ".meanbranchgcconsdsomsuffstat").c_str());
        gcgos << gcconsglobdsomss << '\n';
        cerr << "global GC-cons dsom path suffstats in " << name << ".meanbranchgcconsdsomsuffstat\n";

        GetModel()->MasterReceiveGeneArray(gcarray);
        GCCodonPathSuffStatBranchArray globgcss(GetModel()->GetTree());
        for (int i=0; i<GetModel()->GetNgene(); i++) {
            globgcss.Add(gcarray[i]);
        }

        SimpleBranchArray<double> gctree(GetModel()->GetTree(),0);
        SimpleBranchArray<double> dstree(GetModel()->GetTree(),0);

        globdsomss.GetdS(dstree);
        globgcss.GetGC(gctree);

        ofstream gcnucos((name + ".meanbranchgc.tre").c_str());
        WritedSGCTree(gcnucos, dstree, gctree);
        cerr << "empirical GC in newick format in " << name << ".meanbranchgc.tre\n";

        ofstream tabgcnucos((name + ".leafbranchgc.tab").c_str());
        WriteGCLeafTab(tabgcnucos, gctree);
        cerr << "empirical GC tabulated for terminal branches in " << name << ".leafbranchgc.tab\n";

    }

    void SlaveReaddSOmegaPathSuffStat() {
        vector<dSOmegaPathSuffStatBranchArray> array(GetModel()->GetLocalNgene(), dSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        vector<GCConsdSOmegaPathSuffStatBranchArray> gcconsarray(GetModel()->GetLocalNgene(), GCConsdSOmegaPathSuffStatBranchArray(GetModel()->GetTree()));
        vector<GCCodonPathSuffStatBranchArray> gcarray(GetModel()->GetLocalNgene(), GCCodonPathSuffStatBranchArray(GetModel()->GetTree()));
        for (int i = 0; i < size; i++) {
            GetNextPoint();
            GetModel()->SlaveUpdate();
            GetModel()->SlaveAdddSOmegaPathSuffStat(array);
            GetModel()->SlaveAddGCConsdSOmegaPathSuffStat(gcconsarray);
            GetModel()->SlaveAddGCCodonPathSuffStat(gcarray);
        }
        for (int i=0; i<GetModel()->GetLocalNgene(); i++) {
            array[i].Normalize(1.0/size);
            gcconsarray[i].Normalize(1.0/size);
            gcarray[i].Normalize(1.0/size);
        }
        GetModel()->SlaveSendGeneArray(array);
        GetModel()->SlaveSendGeneArray(gcconsarray);
        GetModel()->SlaveSendGeneArray(gcarray);
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
    int nuc = 0;

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
            } else if (s == "-nuc") {
                nuc = 1;
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
    } else if (nuc) {
        if (! myid) {
            sample->MasterReadNucRates();
        }
        else    {
            sample->SlaveReadNucRates();
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
