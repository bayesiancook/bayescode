#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;

int main(int argc, char *argv[]) {
    CmdLine cmd{"ReadSingleOmega", ' ', "0.1"};
    ReadArgParse read_args(cmd);
    cmd.parse(argc, argv);

    std::string chain_name = read_args.GetChainName();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    std::ifstream is{chain_name + ".param"};
    ChainDriver *fake_read = nullptr;
    unique_ptr<SingleOmegaModel> model = nullptr;
    fake_read = new ChainDriver(is);
    is >> model;
    ChainReader cr(*model, chain_name + ".chain");

    cr.skip(read_args.GetBurnIn());
    if (read_args.GetPpred()) {
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model->PostPred("ppred_" + chain_name + "_" + std::to_string(i) + ".ali");
        }
        cerr << '\n';
    } else {
        stats_posterior<SingleOmegaModel>(*model, cr, every, size);
    }
}
