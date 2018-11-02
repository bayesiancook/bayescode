#pragma once

#include <set>
#include <string>

class Partition;

template <class T>
class RegistrarBase {
    mutable std::set<std::string> _filter;

  public:
    void register_element(std::string, double&) {
        std::cerr << "Error: register_element(std::string, double&) not implemented!\n";
        exit(1);
    }
    void register_element(std::string, std::vector<double>&) {
        std::cerr
            << "Error: register_element(std::string, std::vector<double>&) not implemented!\n";
        exit(1);
    }
    void register_element(std::string, std::vector<double>&, const Partition&) {
        std::cerr << "Error: register_element(std::string, std::vector<double>&, const Partition&) "
                     "not implemented!\n";
        exit(1);
    }

    RegistrarBase() = default;

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