#pragma once

#include "components/BaseArgParse.hpp"

class TreeAppArgParse : public BaseArgParse {
  public:
    explicit TreeAppArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}
    ValueArg<std::string> treefile{"t", "tree", "File path to the tree (NHX format).", true, "", "string", cmd};
    ValueArg<int> every{
        "e", "every", "Number of MCMC iterations between two saved point in the trace.", false, 1, "int", cmd};
    ValueArg<int> until{"u", "until", "Maximum number of (saved) iterations (-1 means unlimited).",
                        false, -1, "int", cmd};
    SwitchArg force{"f", "force", "Overwrite existing output files.", cmd};
};

class InferenceAppArgParse : public TreeAppArgParse {
  public:
    explicit InferenceAppArgParse(ChainCmdLine &cmd) : TreeAppArgParse(cmd) {}
    ValueArg<std::string> alignment{
        "a", "alignment", "File path to alignment (PHYLIP format).", true, "", "string", cmd};
};
