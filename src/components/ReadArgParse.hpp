#pragma once

#include <cassert>
#include <fstream>
#include "tclap/CmdLine.h"

class ReadArgParse {
  protected:
    TCLAP::CmdLine& cmd;

  public:
    explicit ReadArgParse(TCLAP::CmdLine& cmd) : cmd{cmd} {}

    TCLAP::ValueArg<int> every{
        "e", "every", "Number of iterations between two traces", false, 1, "int", cmd};
    TCLAP::ValueArg<int> until{"u", "until",
        "Maximum number of (saved) iterations (-1 means unlimited)", false, -1, "int", cmd};
    TCLAP::ValueArg<int> burnin{
        "b", "burnin", "Number of iterations for burnin", false, 0, "int", cmd};
    TCLAP::SwitchArg ppred{"p", "ppred",
        "For each point of the chain (after burn-in), produces a data replicate simulated "
        "from the posterior predictive distribution",
        cmd};
    TCLAP::UnlabeledValueArg<std::string> chain_name{
        "chain_name", "Chain name (output file prefix)", true, "chain", "string", cmd};

    std::string GetChainName() { return chain_name.getValue(); }

    bool GetPpred() { return ppred.getValue(); }

    int GetBurnIn() { return burnin.getValue(); }

    int GetEvery() { return every.getValue(); }

    int GetSize() { return (GetUntil() - GetBurnIn()) / GetEvery(); }

    int GetUntil() {
        int until_value{until.getValue()};
        if (until_value == -1) {
            std::ifstream trace{GetChainName() + ".trace"};
            std::string line;
            while (std::getline(trace, line)) {
                if (line != "\n") { until_value++; }
            }
        }
        assert(until_value > 0);
        return until_value;
    }
};
