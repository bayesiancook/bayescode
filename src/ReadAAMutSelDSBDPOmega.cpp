#include <cmath>
#include <fstream>
#include "AAMutSelDSBDPOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;

class ReadAAMutSelDSBDPOmegaArgParse {
    CmdLine &cmd;

public:
    ReadAAMutSelDSBDPOmegaArgParse(CmdLine &cmd) : cmd(cmd) {}

    ValueArg<int> every{
            "e", "every", "Number of iterations between two traces", false, 1, "int", cmd};
    ValueArg<int> until{"u", "until", "Maximum number of (saved) iterations (-1 means unlimited)",
                        false, -1, "int", cmd};
    ValueArg<int> burnin{"b", "burnin", "Number of iterations for burnin", false, 0, "int", cmd};
    SwitchArg ppred{"p", "ppred", "For each point of the chain (after burn-in), produces a data replicate simulated "
                                  "from the posterior predictive distribution", cmd};
    SwitchArg om{"o", "om", "Computes the mean predicted omega under mutation-selection balance", cmd};
    SwitchArg ss{"s", "ss", "Computes the mean posterior site-specific state equilibrium frequencies", cmd};
    SwitchArg stats{"s", "stats", "Computes the mean ", cmd};
    UnlabeledValueArg<std::string> chain_name{
            "chain_name", "Chain name (output file prefix)", true, "chain", "string", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"AAMUTSELDSBDPOMEGA", ' ', "0.1"};
    ReadAAMutSelDSBDPOmegaArgParse args(cmd);
    cmd.parse(argc, argv);

    int burnin = args.burnin.getValue();
    int every = args.every.getValue();
    int until = args.until.getValue();
    int size = (until - burnin) / every;
    std::string chain_name = args.chain_name.getValue();

    std::ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is);  // We're not interested in the ChainDriver of the param file
    AAMutSelDSBDPOmegaModel model{is};
    ChainReader cr{model, chain_name + ".chain"};

    cr.skip(burnin);
    if (args.ppred.getValue()) {
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model.PostPred("ppred_" + chain_name + "_" + std::to_string(i) + ".ali");
        }
        cerr << '\n';
    } else if (args.om.getValue()) {
        cerr << size << " points to read\n";

        double meandnds = 0;
        double vardnds = 0;

        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            double om = model.GetPredictedDNDS();
            meandnds += om;
            vardnds += om * om;
        }
        cerr << '\n';
        meandnds /= size;
        vardnds /= size;
        vardnds -= meandnds * meandnds;

        cout << "posterior mean omega : " << meandnds << '\t' << sqrt(vardnds) << '\n';
    } else if (args.ss.getValue()) {
        cerr << size << " points to read\n";

        double **sitestat = new double *[model.GetNsite()];
        for (int i = 0; i < model.GetNsite(); i++) {
            sitestat[i] = new double[model.GetNsite()];
            for (int k = 0; k < model.GetNsite(); k++) {
                sitestat[i][k] = 0;
            }
        }

        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model.GetNsite(); i++) {
                // TO DO !!!!!
                // double *p = model.GetProfile(i);
                for (int k = 0; k < model.GetNsite(); k++) {
                    // sitestat[i][k] += p[k];
                    sitestat[i][k] += 0;
                }
            }
        }
        cerr << '\n';

        ofstream os((chain_name + ".siteprofiles").c_str());
        os << model.GetNsite() << '\n';
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1;
            for (int k = 0; k < model.GetNsite(); k++) {
                sitestat[i][k] /= size;
                os << '\t' << sitestat[i][k];
            }
            os << '\n';
        }
        cerr << "mean site-specific profiles in " << chain_name << ".siteprofiles\n";
        cerr << '\n';
    }
}