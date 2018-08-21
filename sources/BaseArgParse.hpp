#pragma once
#include "tclap/CmdLine.h"

using namespace TCLAP;

class BaseArgParse {
  protected:
    CmdLine& cmd;

  public:
    BaseArgParse(CmdLine& cmd) : cmd(cmd) {}
};
