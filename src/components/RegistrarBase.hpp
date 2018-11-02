#pragma once

#include <set>
#include <string>

template <class T>
class RegistrarBase {
    mutable std::set<std::string> _filter;

  public:
    RegistrarBase() {}

    template <class... Args>
    void add(std::string name, Args&&... args) {
        if ((_filter.size() == 0) or (_filter.find(name) != _filter.end())) {
            static_cast<T*>(this)->register_element(name, std::forward<Args>(args)...);
        }
    }

    template <class M>
    void register_from_method(
        M& ref, void (M::*f)(RegistrarBase<T>&), std::set<std::string> filter = {}) {
        _filter = filter;
        (ref.*f)(*this);
        _filter.clear();
    }
};