#pragma once

#include <functional>
#include <memory>
#include <vector>
#include "interfaces.hpp"
#include "mpi_components/Process.hpp"

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
        if (operation.get() != nullptr) {
            this->operations.insert(this->operations.begin(), std::move(operation));
        }
    }

    void add(std::unique_ptr<Proxy> ptr) {
        if (ptr.get() != nullptr) { operations.push_back(std::move(ptr)); }
    }

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

/*
====================================================================================================
  Operation class
==================================================================================================*/
class Operation : public Proxy {
    std::function<void()> f_acquire{[]() {}};
    std::function<void()> f_release{[]() {}};

  public:
    template <class Acquire, class Release>
    Operation(Acquire f_acquire, Release f_release) : f_acquire(f_acquire), f_release(f_release) {}

    void acquire() final { f_acquire(); }
    void release() final { f_release(); }
};

template <class Acquire, class Release>
std::unique_ptr<Proxy> make_operation(Acquire f_acquire, Release f_release) {
    return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new Operation(f_acquire, f_release)));
}

template <class Acquire>
std::unique_ptr<Proxy> make_acquire_operation(Acquire f_acquire) {
    return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new Operation(f_acquire, []() {})));
}

template <class Release>
std::unique_ptr<Proxy> make_release_operation(Release f_release) {
    return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new Operation([]() {}, f_release)));
}

/*
====================================================================================================
  ForInContainer class
==================================================================================================*/
template <class Container>
class ForInContainer : public Proxy {
    using Element = typename Container::value_type;
    Container& container;
    std::function<Proxy&(Element&)> get_proxy;

  public:
    template <class GetProxy>
    ForInContainer(Container& container, GetProxy get_proxy)
        : container(container), get_proxy(get_proxy) {}

    ForInContainer(Container& container)
        : container(container), get_proxy([](Element& e) -> Proxy& { return e; }) {}

    void acquire() final {
        for (auto& element : container) { get_proxy(element).acquire(); }
    }

    void release() final {
        for (auto& element : container) { get_proxy(element).release(); }
    }
};

template <class Container, class GetProxy>
std::unique_ptr<Proxy> make_for_in_container(Container& container, GetProxy get_proxy) {
    return std::unique_ptr<Proxy>(
        dynamic_cast<Proxy*>(new ForInContainer<Container>(container, get_proxy)));
}

template <class Container>
std::unique_ptr<Proxy> make_for_in_container(Container& container) {
    return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(new ForInContainer<Container>(container)));
}

std::unique_ptr<Proxy> slave_only(std::unique_ptr<Proxy>&& ptr) {
    return MPI::p->rank ? std::move(ptr) : nullptr;
}

std::unique_ptr<Proxy> master_only(std::unique_ptr<Proxy>&& ptr) {
    return (!MPI::p->rank) ? std::move(ptr) : nullptr;
}