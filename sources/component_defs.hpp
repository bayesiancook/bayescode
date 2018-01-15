#ifndef COMPONENT_DEFS_H
#define COMPONENT_DEFS_H

#include "tinycompo.hpp"

template <class Type>
struct Wrapper : tc::Component, public Type {
    template <class... Args>
    Wrapper(Args... args) : Type(args...) {}
};

#endif
