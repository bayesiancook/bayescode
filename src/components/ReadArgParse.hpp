#pragma once

#include <cassert>
#include <fstream>
#include "tclap/CmdLine.h"

class ReadArgParse {
  protected:
    TCLAP::CmdLine& cmd;
    int until{-1};
    int burnin{-1};

  public:
    explicit ReadArgParse(TCLAP::CmdLine& cmd) : cmd{cmd} {}

    TCLAP::ValueArg<int> every{"e", "every",
        "Number of MCMC iterations between two saved point in the trace.", false, 1, "int", cmd};
    TCLAP::ValueArg<int> until_input{"u", "until",
        "Maximum number of (saved) iterations (-1 means unlimited).", false, -1, "int", cmd};
    TCLAP::ValueArg<int> burnin_input{
        "b", "burnin", "Number of MCMC iterations for the burn-in.", false, 0, "int", cmd};
    TCLAP::SwitchArg ppred{"p", "ppred",
        "For each point of the chain (after burn-in), produces a data replicate simulated "
        "from the posterior predictive distribution",
        cmd};
    TCLAP::SwitchArg trace{"", "trace",
        "Recompute the trace."
        "Trace is written in {chain_name}.trace.tsv by default (optionally "
        "use the --output argument to specify a different output path).",
        cmd};
    TCLAP::UnlabeledValueArg<std::string> chain_name{
        "chain_name", "Chain name (output file prefix).", true, "chain", "string", cmd};
    TCLAP::ValueArg<string> output{
        "o", "output", "Output file path (optional)", false, "", "string", cmd};

    std::string GetChainName() { return chain_name.getValue(); }

    bool GetPpred() { return ppred.getValue(); }

    int GetBurnIn() {
        if (burnin == -1) {
            if (int((GetUntil() - burnin_input.getValue()) / GetEvery()) < 1) {
                std::cerr << "Burn-in was set to " << burnin_input.getValue()
                          << " but there is only " << GetUntil() << "  points to read.\n";
                burnin = std::max(GetUntil() - 10 * GetEvery(), 0);
                std::cerr << "Setting the burn-in to " << burnin << std::endl;
            } else {
                burnin = burnin_input.getValue();
            }
        }
        return burnin;
    }

    int GetEvery() { return every.getValue(); }

    int GetSize() { return (GetUntil() - GetBurnIn()) / GetEvery(); }

    int GetUntil() {
        if (until == -1) {
            std::ifstream trace{GetChainName() + ".trace"};
            std::string line;
            while (std::getline(trace, line)) {
                if (line != "\n") { until++; }
            }
            if (until_input.getValue() != -1) { until = std::min(until, until_input.getValue()); }
            assert(until > 0);
        }
        return until;
    }


    std::string OutputFile(const std::string& default_suffix = "") {
        if (!output.getValue().empty()) {
            return output.getValue();
        } else {
            std::cout << "No output file name specified, using default: "
                      << GetChainName() + default_suffix << std::endl;
            return GetChainName() + default_suffix;
        }
    }
};
