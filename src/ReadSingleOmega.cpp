#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
#include "tclap/CmdLine.h"
#include "components/ChainReader.hpp"
#include "components/ChainDriver.hpp"
using namespace std;
using namespace TCLAP;

class ReadSingleOmegaArgParse {
    CmdLine &cmd;

public:
    ReadSingleOmegaArgParse(CmdLine &cmd) : cmd(cmd) {}
    ValueArg<int> every{"e",   "every", "Number of iterations between two traces", false, 1,
                        "int", cmd};
    ValueArg<int> until{"u",   "until", "Maximum number of (saved) iterations (-1 means unlimited)",
                        false, -1,      "int",
                        cmd};
    ValueArg<int> burnin{"b",   "burnin", "Number of iterations for burnin",
                        false, 0,      "int",
                        cmd};
    SwitchArg ppred{"p", "ppred", "Perform simulations under posterior distribution", cmd};
    UnlabeledValueArg<std::string> chain_name{
        "chain_name", "Chain name (output file prefix)", true, "chain", "string", cmd};
};


int main(int argc, char *argv[]) {
    CmdLine cmd{"ReadSingleOmega", ' ', "0.1"};
    ReadSingleOmegaArgParse args(cmd);
    cmd.parse(argc, argv);

    int burnin = args.burnin.getValue();
    int every = args.every.getValue();
    int until = args.until.getValue();
    int ppred = args.ppred.getValue();
    int size = (until - burnin) / every;
    std::string chain_name = args.chain_name.getValue();

    std::ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is); // We're not interested in the ChainDriver of the param file
    SingleOmegaModel model{is};
    ChainReader cr{model, chain_name + ".chain"};

    cr.skip(burnin);
    if (ppred) {
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model.PostPred("ppred_" + chain_name + "_" + std::to_string(i) + ".ali");
        }
        cerr << '\n';
    } else {
        cerr << size << " points to read\n";

        double meanomega = 0;
        double varomega = 0;

        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            double om = model.GetOmega();
            meanomega += om;
            varomega += om * om;
        }
        cerr << '\n';
        meanomega /= size;
        varomega /= size;
        varomega -= meanomega * meanomega;

        cout << "posterior mean omega : " << meanomega << '\t' << sqrt(varomega) << '\n';
    }
}
