#pragma once

#include "global/defs.hpp"

namespace param {
    enum mode_t { fixed, shared, shrunk, independent, invalid };

    istream &operator>>(istream &is, param::mode_t &mode) {
        string input;
        is >> input;
        if (input == "fixed") {
            mode = fixed;
        } else if (input == "shrunk") {
            mode = shrunk;
        } else if (input == "shared") {
            mode = shared;
        } else if (input == "independent") {
            mode = independent;
        } else {
            FAIL("Could not recognize '{}' as a valid parameter mode.", input);
        }
        return is;
    }

    template <class ValueT, class HyperT>
    struct Config {
        ValueT value;
        HyperT hyper;
        mode_t mode{invalid};
    };

}  // namespace param