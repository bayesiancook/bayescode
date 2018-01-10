#ifndef COMPONENT_DEFS_H
#define COMPONENT_DEFS_H

#include "tinycompo.hpp"

template <class Type>
struct Wrapper : tc::Component {
    Type data;
};

#endif
