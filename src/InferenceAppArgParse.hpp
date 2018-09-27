#include "components/BaseArgParse.hpp"

class InferenceAppArgParse : public BaseArgParse {
  public:
    InferenceAppArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}
    ValueArg<string> alignment{
        "a", "alignment", "Alignment file (PHYLIP)", true, "", "string", cmd};
    ValueArg<string> treefile{"t", "tree", "Tree file (NHX)", true, "", "string", cmd};
    ValueArg<int> every{
        "e", "every", "Number of iterations between two traces", false, 1, "int", cmd};
    ValueArg<int> until{"u", "until", "Maximum number of (saved) iterations (-1 means unlimited)",
        false, -1, "int", cmd};
    SwitchArg force{"f", "force", "Overwrite existing output files", cmd};
};
