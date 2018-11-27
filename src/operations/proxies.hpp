#pragma once

#include <vector>
#include "interfaces.hpp"

/*
====================================================================================================
  ForAll class
  Propagates acquire and release to all proxies registered to it
==================================================================================================*/
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

template <class... Args>
std::unique_ptr<Proxy> make_forall(Args&&... args) {
    return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new ForAll(std::forward<Args>(args)...)));
}

/*
====================================================================================================
  Group class
  Same as forall but owns its contents in the form of unique pointers
==================================================================================================*/
class Group : public Proxy {
    std::vector<std::unique_ptr<Proxy>> operations;

  public:
    Group() = default;

    template <class... Operations>
    Group(std::unique_ptr<Proxy>&& operation, Operations&&... operations)
        : Group(std::forward<Operations>(operations)...) {
        /* -- */
        this->operations.emplace_back(std::move(operation));
    }

    void add(std::unique_ptr<Proxy> ptr) { operations.push_back(std::move(ptr)); }

    void acquire() final {
        for (auto&& operation : operations) { operation->acquire(); }
    }

    void release() final {
        for (auto&& operation : operations) { operation->release(); }
    }
};

template <class... Args>
std::unique_ptr<Proxy> make_group(Args&&... args) {
    return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new Group(std::forward<Args>(args)...)));
}
