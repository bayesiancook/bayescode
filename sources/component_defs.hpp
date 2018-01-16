#ifndef COMPONENT_DEFS_H
#define COMPONENT_DEFS_H

#include "tinycompo.hpp"

// Wrapper for classes (with virtual destructors)
template <class Type>
struct Wrapper : tc::Component, public Type {
    template <class... Args>
    Wrapper(Args... args) : Type(args...) {}
};

// Wrapper for fundamental types or other classes
template <class Type>
struct FWrapper : tc::Component {
    Type data;

    template <class... Args>
    FWrapper(Args... args) : data(args...) {}

    Type& operator*() { return data; }
};

#endif
