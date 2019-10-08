#pragma once

#include "lib/BranchArray.hpp"

class LegacyArrayProxy : public BranchSelector<double> {
    std::vector<double>& data_ref;
    const Tree& tree_ref;

  public:
    LegacyArrayProxy(std::vector<double>& data_ref, const Tree& tree_ref)
        : data_ref(data_ref), tree_ref(tree_ref) {}

    virtual const Tree& GetTree() const override { return tree_ref; }
    virtual const double& GetVal(int index) const override { return data_ref[index]; }
};