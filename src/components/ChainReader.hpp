#pragma once

#include "Tracer.hpp"

class ChainReader {
    Tracer tracer;
    std::ifstream is;

  public:
    template <class M>
    ChainReader(M& model, std::string filename) : tracer(model), is(filename) {
        tracer.ignore_header(is);
    }

    void next() { tracer.read_line(is); }
    void skip(int n) {
        for (int i = 0; i < n; i++) next();
    }
};
