#pragma once

#include <set>
#include <string>

template <class T>
class RegistrarBase {
    std::set<std::string> filter;

  public:
    RegistrarBase(std::set<std::string> filter) : filter(filter) {}

    template <class... Args>
    void add(std::string name, Args&&... args) {
        if (filter.find(name) != filter.end()) {
            static_cast<T*>(this)->register_element(name, std::forward<Args>(args)...);
        }
    }

    template <class M>
    void register_from_method(M* ptr, void (M::*f)(RegistrarBase<T>&)) {
        (ptr->*f)(*this);
    }
};