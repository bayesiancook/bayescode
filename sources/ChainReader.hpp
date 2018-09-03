#pragma once

#include "Tracer.hpp"

class ChainReader {
    Tracer tracer;
    std::ifstream is;

public:
    template<class M>
    ChainReader(M& model, std::string filename) : tracer(model, &M::declare_model),
                                                  is(filename) {
        tracer.ignore_header();
    }

    void next() { tracer.read_line(); }
};
