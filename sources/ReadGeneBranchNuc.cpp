#include <cmath>
#include <fstream>
#include <algorithm>
#include "Sample.hpp"
#include "GeneBranchNucModel.hpp"
// #include "RecursiveNewick.hpp"

using namespace std;

/**
 * \brief An MCMC sample for GeneBranchNucModel
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
    int syn_devmode, nuc_devmode;

  public:
    string GetModelType() override { return modeltype; }

    GeneBranchNucModel *GetModel() override { return (GeneBranchNucModel *)model; }

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
        is >> syn_devmode >> nuc_devmode;
        int check;
        is >> check;
        if (check)  {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "GENEBRANCHNUC") {
            model = new GeneBranchNucModel(datafile, treefile, taxonfile, syn_devmode, nuc_devmode);
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
        vector<double> mean_branchsyn_array(Nbranch,0);
        vector<double> mean_nucss_branchsyn_array(Nbranch,0);

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
            GetModel()->AddEmpiricalBranchLengthsFromNucPathSuffStat(mean_nucss_branchsyn_array);
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

        double totlength = 0;
        double totlength2 = 0;
        for (int j=0; j<Nbranch; j++)   {
            mean_branchsyn_array[j] /= size;
            totlength2 += mean_branchsyn_array[j];
            mean_nucss_branchsyn_array[j] /= size;
            totlength += mean_nucss_branchsyn_array[j];
        }

        ofstream tos((name + ".postmeansyn.tre").c_str());
        ToNewick(tos, *GetModel()->GetTree(), mean_branchsyn_array, mean_nucss_branchsyn_array);
        // ToNewick(tos, *GetModel()->GetTree(), mean_nucss_branchsyn_array, mean_branchsyn_array);
        cerr << "post mean tree in " << name << ".postmeansyn.tre\n";
        cerr << "total tree length : " << totlength << '\t' << totlength2 << '\n';
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;

    string name;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if ((s == "-x") || (s == "-extract")) {
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
    sample->Read();
}
