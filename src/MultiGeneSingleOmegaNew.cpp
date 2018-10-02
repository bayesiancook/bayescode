#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneSingleOmegaModel.hpp"
#include "InferenceAppArgParse.hpp"
#include "components/ChainDriver.hpp"

using namespace std;

class MultiGeneSingleOmegaArgParse : public BaseArgParse {
  public:
    MultiGeneSingleOmegaArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd), app(cmd) {}
    InferenceAppArgParse app;
    ValueArg<string> omega{"","omega","omega hyperparameter, expected syntax: <float>, <float>",false,"uninf","string",cmd};
    ValueArg<string> nucrates{"","nucrates","Possible values are, shared, shrunken, independent(ind)",false,"shrunken","string",cmd};
    ValueArg<string> bl{"","blmode","Possible values are, shared, shrunken, independent(ind)",false,"shrunken","string",cmd};

    int blmode() { return decode_mode("bl", bl.getValue()); }
    int nucmode() { return decode_mode("nucrates", bl.getValue()); }
    int omegamode(double& omegahypermean, double& omegahyperinvshape) {
        if (omega.getValue() != "uninf") {
            std::stringstream ss(omega.getValue());
            std::string item;
            std::getline(ss, item, ',');
            omegahypermean = atof(item.c_str());
            std::getline(ss, item, ',');
            omegahyperinvshape = atof(item.c_str());
            return 0;
        }
        else {
            omegahypermean = 1.0;
            omegahyperinvshape = 1.0;
            return 1;
        }
    }

private:
    int decode_mode(std::string option_name, std::string mode) {
        if (mode == "shared") {
            return 2;
        } else if (mode == "shrunken") {
            return 1;
        } else if ((mode == "ind") || (mode == "independent")) {
            return 0;
        } else {
            cerr << "error: does not recongnize command after -" << option_name << "\n";
            exit(1);
        }
    }
};


int main(int argc, char *argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    ChainCmdLine cmd{argc, argv, "MultiGeneSingleOmega", ' ', "0.1"};

    ChainDriver *chain_driver = nullptr;
    MultiGeneSingleOmegaModel *model = nullptr;

    if (cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = new ChainDriver(is);
        // model = new SingleOmegaModel(is);
    } else {
        MultiGeneSingleOmegaArgParse args(cmd);
        cmd.parse();
        chain_driver =
            new ChainDriver(cmd.chain_name(), args.app.every.getValue(), args.app.until.getValue());
        model = new MultiGeneSingleOmegaModel(args.app.alignment.getValue(),
                                              args.app.treefile.getValue(),
                                              myid, nprocs);

        double omegahypermean;
        double omegahyperinvshape;
        model->SetAcrossGenesModes(args.blmode(),
                                   args.nucmode(),
                                   args.omegamode(omegahypermean, omegahyperinvshape));
        model->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape);
        model->Allocate();
        model->Update();
    }

    MPI_Finalize();
}
