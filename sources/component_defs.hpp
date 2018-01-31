#ifndef COMPONENT_DEFS_H
#define COMPONENT_DEFS_H

#include "Random.hpp"
#include "tinycompo.hpp"

//==========================================================
// INTERFACES
//==========================================================

class Go : public tc::Component {
    virtual void go() = 0;

  public:
    Go() { port("go", &Go::go); }
};

template <class Type>
struct AbstractWrapper {
    virtual Type operator*() = 0;
};

struct Lifecycle {
    virtual void Init() = 0;
    virtual void EndMove() = 0;
    virtual void End() = 0;
};

struct AbstractTraceFile {
    virtual void write_header() = 0;
    virtual void write_line() = 0;
};

//==========================================================
// WRAPPERS
//==========================================================

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

//==========================================================
// FILE HANDLERS
//==========================================================
template <class Type>
class TraceFile : public AbstractTraceFile, public tc::Component {
    ofstream fs;
    void (Type::*trace_method)(ostream&) const;
    void (Type::*header_method)(ostream&) const;
    Type* target{nullptr};

  public:
    TraceFile(string filename, void (Type::*header_method)(ostream&) const, void (Type::*trace_method)(ostream&) const,
              bool erase_contents = true)
        : trace_method(trace_method), header_method(header_method) {
        port("target", &TraceFile::target);
        if (erase_contents) {  // TODO warning/error if file is not empty
            fs.open(filename, ios_base::trunc);
        } else {
            fs.open(filename, ios_base::app);
        }
    }

    void write_header() override { (target->*header_method)(fs); }

    void write_line() override { (target->*trace_method)(fs); }
};

class RunToggle : public tc::Component {
    fstream fs;

    void set(int i) {
        fs << i;
        fs.seekg(fs.beg);
    }

  public:
    RunToggle(string chainname) {
        fs.open(chainname + ".run", ios_base::in | ios_base::out | ios_base::trunc);
        set(1);
    }

    void toggle() {
        if (fs.peek() == '1') {
            set(0);
        } else {
            set(1);
        }
    }

    bool check() { return fs.peek() == '1'; }
};

//==========================================================
// CONNECTORS
//==========================================================

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
