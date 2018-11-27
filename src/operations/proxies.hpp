#pragma once

#include <vector>
#include "interfaces.hpp"

// does not own its contents
class ForAll : public Proxy {
    std::vector<Proxy*> pointers;

  public:
    ForAll() = default;

    template <class... Pointers>
    ForAll(Pointers&&... pointers) : pointers({pointers...}) {}

    void add(Proxy* pointer) { pointers.push_back(pointer); }

    void acquire() final {
        for (auto pointer : pointers) { pointer->acquire(); }
    }

    void release() final {
        for (auto pointer : pointers) { pointer->release(); }
    }
};