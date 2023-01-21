#include <cmath>
#include <fstream>
#include <algorithm>
#include "Sample.hpp"
#include "GeneBranchStrandSymmetricCodonModel.hpp"
// #include "RecursiveNewick.hpp"

using namespace std;

/**
 * \brief An MCMC sample for GeneBranchStrandSymmetricCodonModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class GeneBranchStrandSymmetricSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    string taxonfile;
    int syn_devmode, om_devmode, nuc_devmode;

  public:
    string GetModelType() override { return modeltype; }

    GeneBranchStrandSymmetricCodonModel *GetModel() override { return (GeneBranchStrandSymmetricCodonModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    GeneBranchStrandSymmetricSample(string filename, int inburnin, int inevery, int inuntil)
        : Sample(filename, inburnin, inevery, inuntil) {
        Open();
    }

    void Open() override {
        // open <name>.param
        ifstream is((name + ".param").c_str());

        // check that file exists
        if (!is) {
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }

        // read model type, and other standard fields
        is >> modeltype;
        is >> datafile >> treefile >> taxonfile;
        is >> syn_devmode >> om_devmode >> nuc_devmode;
        int check;
        is >> check;
        if (check)  {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "GENEBRANCHCODON") {
            model = new GeneBranchStrandSymmetricCodonModel(datafile, treefile, taxonfile, syn_devmode, om_devmode, nuc_devmode);
        } else {
            cerr << "error when opening file " << name << '\n';
            exit(1);
        }

        // read model (i.e. chain's last point) from <name>.param
        model->FromStream(is);

        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    //! \brief computes the posterior mean estimate (and the posterior standard
    //! deviation) of omega
    void Read() {

        int Nbranch = GetModel()->GetNbranch();
        int Ngene = GetModel()->GetNgene();
        vector<double> mean_branchsyn_array(Nbranch,0);
        vector<double> mean_branchom_array(Nbranch,0);
        vector<double> mean_genesyn_array(Ngene,0);
        vector<double> mean_geneom_array(Ngene,0);

        vector<double> mean_branchgc_array(Nbranch,0);
        vector<double> mean_genegc_array(Ngene,0);

        vector<vector<double>> syn_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> syn_val(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_val(Ngene, vector<double>(Nbranch,0));

        vector<vector<double>> gc_val(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> gc_z(Ngene, vector<double>(Nbranch,0));

        double meanAC, meanAG, meanCA, meanCG, meanCT;
        double geneAC, geneAG, geneCA, geneCG, geneCT;
        double branchAC, branchAG, branchCA, branchCG, branchCT;
        double devAC, devAG, devCA, devCG, devCT;

        meanAC = meanAG = meanCA = meanCG = meanCT = 0;
        geneAC = geneAG = geneCA = geneCG = geneCT = 0;
        branchAC = branchAG = branchCA = branchCG = branchCT = 0;
        devAC = devAG = devCA = devCG = devCT = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->GetSynModel()->AddBranchArrayTo(mean_branchsyn_array);
            GetModel()->GetOmegaModel()->AddBranchArrayTo(mean_branchom_array);
            GetModel()->GetSynModel()->AddGeneArrayTo(mean_genesyn_array);
            GetModel()->GetOmegaModel()->AddGeneArrayTo(mean_geneom_array);

            GetModel()->GetSynModel()->AddValTo(syn_val);
            GetModel()->GetOmegaModel()->AddValTo(om_val);
            GetModel()->GetSynModel()->AddZscoreTo(syn_z);
            GetModel()->GetOmegaModel()->AddZscoreTo(om_z);

            GetModel()->AddGCBiasTo(gc_val, gc_z, 
                    mean_branchgc_array, mean_genegc_array);

            GetModel()->AddNucStats(meanAC, meanAG, meanCA, meanCG, meanCT,
                    geneAC, geneAG, geneCA, geneCG, geneCT,
                    branchAC, branchAG, branchCA, branchCG, branchCT,
                    devAC, devAG, devCA, devCG, devCT);
        }
        cerr << '\n';

        meanAC /= size;
        meanAG /= size;
        meanCA /= size;
        meanCG /= size;
        meanCT /= size;

        geneAC /= size;
        geneAG /= size;
        geneCA /= size;
        geneCG /= size;
        geneCT /= size;

        branchAC /= size;
        branchAG /= size;
        branchCA /= size;
        branchCG /= size;
        branchCT /= size;

        devAC /= size;
        devAG /= size;
        devCA /= size;
        devCG /= size;
        devCT /= size;

        double totAC = devAC + branchAC + geneAC + branchAC*geneAC;
        double pbranchAC = branchAC / totAC;
        double pgeneAC = geneAC / totAC;
        double pdevAC = devAC / totAC;

        double totAG = devAG + branchAG + geneAG + branchAG*geneAG;
        double pbranchAG = branchAG / totAG;
        double pgeneAG = geneAG / totAG;
        double pdevAG = devAG / totAG;

        double totCA = devCA + branchCA + geneCA + branchCA*geneCA;
        double pbranchCA = branchCA / totCA;
        double pgeneCA = geneCA / totCA;
        double pdevCA = devCA / totCA;

        double totCG = devCG + branchCG + geneCG + branchCG*geneCG;
        double pbranchCG = branchCG / totCG;
        double pgeneCG = geneCG / totCG;
        double pdevCG = devCG / totCG;

        double totCT = devCT + branchCT + geneCT + branchCT*geneCT;
        double pbranchCT = branchCT / totCT;
        double pgeneCT = geneCT / totCT;
        double pdevCT = devCT / totCT;

        ofstream nos((name + ".nucstats").c_str());
        nos << "          \tmean\t\tc-gene\tc-brnch\tc-dev\tc-tot\t\tp-gene\tp-brnch\tp-dev\n";
        nos << fixed << setw(5) << setprecision(2);
        nos << "A:T -> C:G\t" << meanAC << '\t' << '\t' << geneAC << '\t' << branchAC << '\t' << devAC << '\t';
        nos << totAC << '\t' << '\t' << pgeneAC << '\t' << pbranchAC << '\t' << pdevAC << '\n';
        nos << "A:T -> G:C\t" << meanAG << '\t' << '\t' << geneAG << '\t' << branchAG << '\t' << devAG << '\t';
        nos << totAG << '\t' << '\t' << pgeneAG << '\t' << pbranchAG << '\t' << pdevAG << '\n';
        nos << "C:G -> A:T\t" << meanCA << '\t' << '\t' << geneCA << '\t' << branchCA << '\t' << devCA << '\t';
        nos << totCA << '\t' << '\t' << pgeneCA << '\t' << pbranchCA << '\t' << pdevCA << '\n';
        nos << "C:G -> G:C\t" << meanCG << '\t' << '\t' << geneCG << '\t' << branchCG << '\t' << devCG << '\t';
        nos << totCG << '\t' << '\t' << pgeneCG << '\t' << pbranchCG << '\t' << pdevCG << '\n';
        nos << "C:G -> T:A\t" << meanCT << '\t' << '\t' << geneCT << '\t' << branchCT << '\t' << devCT << '\t';
        nos << totCT << '\t' << '\t' << pgeneCT << '\t' << pbranchCT << '\t' << pdevCT << '\n';
        cerr << "post mean nuc stats in " << name << ".nucstats\n";

        for (int i=0; i<Ngene; i++)   {
            for (int j=0; j<Nbranch; j++)   {
                syn_z[i][j] /= size;
                om_z[i][j] /= size;
                syn_val[i][j] /= size;
                om_val[i][j] /= size;
                gc_z[i][j] /= size;
                gc_val[i][j] /= size;
            }
        }

        for (int j=0; j<Nbranch; j++)   {
            mean_branchsyn_array[j] /= size;
            mean_branchom_array[j] /= size;
            mean_branchgc_array[j] /= size;
        }
        for (int i=0; i<Ngene; i++) {
            mean_genesyn_array[i] /= size;
            mean_geneom_array[i] /= size;
            mean_genegc_array[i] /= size;
        }

        ofstream gos((name + ".postmean.genedsomgc.tab").c_str());
        gos << "#genename\tgenesyn\tgeneom\tgenegcbias\n";
        for (int i=0; i<Ngene; i++) {
            gos << GetModel()->GetGeneName(i) << '\t' << mean_genesyn_array[i] << '\t' << mean_geneom_array[i] << '\t' << mean_genegc_array[i] << '\n';
        }
        cerr << "post mean gene dS, dN/dS and eq. GC in " << name << ".postmean.geneom.tab\n";

        ofstream os((name + ".postmean.leafdsomgcbias.tab").c_str());
        os << "#taxon\tbranchsyn\tbranchom\tbranchgcbias\n";
        Tabulate(os, GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array, mean_branchgc_array, true);

        ofstream allos((name + ".postmean.dsomgcbias.tab").c_str());
        allos << "#taxon1\ttaxon2\tbranchsyn\tbranchom\tbranchgcbias\n";
        Tabulate(allos, GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array, mean_branchgc_array, false);

        ofstream tos((name + ".dsom.tre").c_str());
        ToNewick(tos, *GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array);

        cerr << "post mean branch effects on dN/dS in " << name << ".postmean.leafdsomgcbias.tab and " << name << ".postmean.dsomgcbias.tab\n";
        cerr << "newick format in " << name << ".dsom.tre\n";

        ofstream devos((name + ".postmeandev.tab").c_str());
        devos << "#genename";
        devos << "\tterminal\tleft_taxon\tright_taxon";
        devos << "\tbranchsyn\tgenesyn\tsyn\tsynz";
        devos << "\tbranchom\tgeneom\tom\tomz";
        devos << "\tbranchgcbias\tgenegcbias\tgcbias\tgcbiasz";
        devos << "\n";

        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                devos << GetModel()->GetGeneName(i); 
                devos << '\t' << ((GetModel()->GetTree()->GetLeftTaxon(j) == GetModel()->GetTree()->GetRightTaxon(j)) ? 1 : 0);
                devos << '\t' << GetModel()->GetTree()->GetLeftTaxon(j);
                devos << '\t' << GetModel()->GetTree()->GetRightTaxon(j);
                devos << '\t' << mean_branchsyn_array[j] << '\t' << mean_genesyn_array[i] << '\t' << syn_val[i][j] << '\t' << syn_z[i][j]; 
                devos << '\t' << mean_branchom_array[j] << '\t' << mean_geneom_array[i] << '\t' << om_val[i][j] << '\t' << om_z[i][j]; 
                devos << '\t' << mean_branchgc_array[j] << '\t' << mean_genegc_array[i] << '\t' << gc_val[i][j] << '\t' << gc_z[i][j]; 
                devos << '\n';
            }
        }
    }

    void ReadSuffStat() {

        vector<dSOmegaPathSuffStatBranchArray> array(GetModel()->GetNgene(), dSOmegaPathSuffStatBranchArray(*GetModel()->GetTree()));
        vector<GCConsdSOmegaPathSuffStatBranchArray> gcconsarray(GetModel()->GetNgene(), GCConsdSOmegaPathSuffStatBranchArray(*GetModel()->GetTree()));

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->CollectdSOmPathSuffStat();
            GetModel()->CollectGCConsdSOmPathSuffStat();
            GetModel()->GetdSOmPathSuffStat().AddTo(array);
            GetModel()->GetGCConsdSOmPathSuffStat().AddTo(gcconsarray);
        }
        cerr << '\n';

        for (int i=0; i<GetModel()->GetNgene(); i++) {
            array[i].Normalize(1.0/size);
            gcconsarray[i].Normalize(1.0/size);
        }

        ofstream dsom_os((name + ".genedsomsuffstat").c_str());
        dsom_os << GetModel()->GetNgene() << '\n';
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            dsom_os << GetModel()->GetGeneName(gene) << '\n';
            dsom_os << "counts_dS\t";
            array[gene].BranchToNewickSynCount(dsom_os);
            dsom_os << "counts_dS_norm\t";
            array[gene].BranchToNewickSynBeta(dsom_os);
            dsom_os << "counts_dN\t";
            array[gene].BranchToNewickNonSynCount(dsom_os);
            dsom_os << "counts_dN_norm\t";
            array[gene].BranchToNewickNonSynBeta(dsom_os);
        }

        ofstream gcdsom_os((name + ".genegcconsdsomsuffstat").c_str());
        gcdsom_os << GetModel()->GetNgene() << '\n';
        for (int gene=0; gene<GetModel()->GetNgene(); gene++) {
            gcdsom_os << GetModel()->GetGeneName(gene) << '\n';
            gcdsom_os << "counts_dS\t";
            gcconsarray[gene].BranchToNewickSynCount(gcdsom_os);
            gcdsom_os << "counts_dS_norm\t";
            gcconsarray[gene].BranchToNewickSynBeta(gcdsom_os);
            gcdsom_os << "counts_dN\t";
            gcconsarray[gene].BranchToNewickNonSynCount(gcdsom_os);
            gcdsom_os << "counts_dN_norm\t";
            gcconsarray[gene].BranchToNewickNonSynBeta(gcdsom_os);
        }
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int dsomss = 0;

    string name;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if (s == "-dsomss") {
                dsomss = 1;
            }
            else if ((s == "-x") || (s == "-extract")) {
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                burnin = atoi(argv[i]);
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                every = atoi(argv[i]);
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                until = atoi(argv[i]);
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

    GeneBranchStrandSymmetricSample *sample = new GeneBranchStrandSymmetricSample(name, burnin, every, until);
    if (dsomss) {
        sample->ReadSuffStat();
    }
    else    {
        sample->Read();
    }
}
