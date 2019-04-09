#pragma once

#include "ParamConfig.hpp"

template <class ValueT, class HyperValueT>
istream &operator>>(istream &is, param::Config<ValueT, HyperValueT> &c) {
    auto split = [](string input, char sep) {
        vector<string> output;
        auto it_begin = input.begin();
        auto it_end = std::find(input.begin(), input.end(), sep);
        do {
            output.emplace_back(it_begin, it_end);
            it_begin = it_end + 1;
            it_end = std::find(it_begin, input.end(), sep);
        } while (it_begin < input.end());
        return output;
    };

    string declaration;
    is >> declaration;

    auto fail = [&declaration](string msg) {
        FAIL("Parameter mode declaration '{}' invalid.\n\t{}", declaration, msg);
    };

    auto colon_split = split(declaration, ':');
    if (colon_split.size() == 1) {  // case 'mode'
        string mode_str = colon_split.at(0);
        if (mode_str == "shared" or mode_str == "shrunk") {  // 'shared'
            std::stringstream(mode_str) >> c.mode;
        } else {
            fail("Shared and shrunk are the only modes which can be used without parameters.");
        }
    } else if (colon_split.size() == 2) {  // case 'mode:things...'
        string mode_str = colon_split.at(0);
        auto slash_split = split(colon_split.at(1), '/');
        if (slash_split.size() == 1) {  // case 'mode:value(s)'
            string value_str = slash_split.at(0);
            if (mode_str == "shared" or mode_str == "fixed") {
                std::stringstream(mode_str) >> c.mode;
                std::stringstream(value_str) >> c.value;
            } else {
                fail("Shared and fixed are the only modes which can have a parameter value.");
            }
        }
    } else {
        fail("Found two or more occurrences of ':' in declaration.");
    }


    return is;
}