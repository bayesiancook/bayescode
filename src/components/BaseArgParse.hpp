#pragma once

#include <fstream>
#include "tclap/CmdLine.h"

using namespace TCLAP;

class ChainCmdLine {
    CmdLine cmd;
    int argc;
    char** argv;
    bool already_parsed{false};
    UnlabeledValueArg<std::string> chain_name_arg{
        "chain_name", "Chain name (output file prefix)", true, "chain", "string", cmd};

  public:
    template <class... Args>
    ChainCmdLine(int argc, char* argv[], Args&&... args)
        : cmd(std::forward<Args>(args)...), argc(argc), argv(argv) {}

    void parse() {
        if (!resume_from_checkpoint() and !already_parsed) {
            cmd.parse(argc, argv);
            already_parsed = true;
        }
    }

    bool resume_from_checkpoint() const { return (argc == 2 && argv[1][0] != '-'); }

    std::ifstream checkpoint_file() { return std::ifstream(chain_name() + ".param"); }

    CmdLine& get() { return cmd; }

    std::string chain_name() {
        if (resume_from_checkpoint())
            return std::string(argv[1]);
        else if (already_parsed)
            return chain_name_arg.getValue();
        else {
            std::cerr << "Trying to get chain name before argument parsing\n";
            exit(1);
        }
    }
};

class BaseArgParse {
  protected:
    CmdLine& cmd;

  public:
    BaseArgParse(ChainCmdLine& cmd) : cmd(cmd.get()) {}
};
