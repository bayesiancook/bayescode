#ifndef COMPONENT_DEFS_H
#define COMPONENT_DEFS_H

#include "Random.hpp"
#include "tinycompo.hpp"

template <class Type>
struct AbstractWrapper {
    virtual Type operator*() = 0;
};

// Wrapper for classes (with virtual destructors)
template <class Type>
struct Wrapper : tc::Component, public Type {
    template <class... Args>
    Wrapper(Args... args) : Type(args...) {}

    Type& operator*() { return *this; }
};

// Wrapper for fundamental types or other classes
template <class Type>
struct FWrapper : tc::Component {
    Type data;

    template <class... Args>
    FWrapper(Args... args) : data(args...) {}

    Type& operator*() { return data; }

    std::string debug() const override {
        stringstream ss;
        ss << "FWrapper [" << data << "]";
        return ss.str();
    }
};

// Set the value to the value of a wrapped object (if it makes sense)
template <class Type>
struct SetValueTo {
    static void _connect(tc::Assembly& assembly, tc::PortAddress user, tc::Address provider) {
        assembly.at(user.address).set(user.prop, *assembly.at<AbstractWrapper<Type>>(provider));
    }
};

struct DirichletSample {
    static void _connect(tc::Assembly& assembly, tc::Address x, const std::vector<double>& center, double concentration) {
        Random::DirichletSample(*assembly.at<Wrapper<std::vector<double>>>(x), center, concentration);
    }
};

#endif
