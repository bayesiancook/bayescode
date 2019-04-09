#pragma once

#include "param_classes.hpp"

template <class ValueType, class HyperValueType>
istream &operator>>(istream &is, MultiGeneParameter<ValueType, HyperValueType> &c) {
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
    auto colon_split = split(declaration, ':');
    if (colon_split.size() == 1) {
        if (colon_split.at(0) == "shared") {}
    }


    return is;
}