#pragma once
#include <iostream>
#include "BranchArray.hpp"
#include "MultinomialAllocationVector.hpp"
#include "StickBreakingProcess.hpp"

class Tracer {
    std::vector<std::function<void(std::ostream&)>> header_to_stream;
    std::vector<std::function<void(std::ostream&)>> data_to_stream;
    std::vector<std::function<void(std::istream&)>> set_from_stream;

  public:
    template <class T>
    Tracer(T& x, void (T::*f)(Tracer&)) {
        (x.*f)(*this);
    }

    void write_header(std::ostream& os) const {
        int n = header_to_stream.size();
        if (n > 0) {
            header_to_stream.at(0)(os);
            for (int i = 1; i < n; i++) {
                os << "\t";
                header_to_stream.at(i)(os);
            }
        }
    }

    void write_line(std::ostream& os) const {
        int n = data_to_stream.size();
        if (n > 0) {
            os << "\n";
            data_to_stream.at(0)(os);
            for (int i = 1; i < n; i++) {
                os << "\t";
                data_to_stream.at(i)(os);
            }
        }
    }

    void ignore_header(std::istream& is) const {
        std::string s;
        std::getline(is, s);
    }

    void read_line(std::istream& is) const {
        for (auto& f : set_from_stream) f(is);
    }

    void add(std::string name, double& d) {
        header_to_stream.push_back([name](std::ostream& os) { os << name; });
        data_to_stream.push_back([&d](std::ostream& os) { os << d; });
        set_from_stream.push_back([&d](std::istream& is) { is >> d; });
    }

    void add(std::string name, std::vector<double>& v) {
        header_to_stream.push_back([&v, name](std::ostream& os) {
            int n = v.size();
            if (n > 0) {
                os << name << "[0]";
                for (int i = 1; i < n; i++) os << "\t" << name << "[" << i << "]";
            }
        });
        data_to_stream.push_back([&v](std::ostream& os) {
            int n = v.size();
            if (n > 0) {
                os << v.at(0);
                for (int i = 1; i < n; i++) os << "\t" << v.at(i);
            }
        });
        set_from_stream.push_back([&v](std::istream& is) {
            for (auto& e : v) is >> e;
        });
    }

    template <class T>
    void add(std::string name, SimpleArray<T>& v) {
        header_to_stream.push_back([&v, name](std::ostream& os) {
            int n = v.GetSize();
            if (n > 0) {
                os << name << "[0]";
                for (int i = 1; i < n; i++) os << "\t" << name << "[" << i << "]";
            }
        });
        data_to_stream.push_back([&v](std::ostream& os) { os << v; });
        set_from_stream.push_back([&v](std::istream& is) { is >> v; });
    }


    void add(std::string name, StickBreakingProcess& sbp) {
        add(name + "_array", dynamic_cast<SimpleArray<double>&>(sbp));

        auto& kappa = sbp.GetKappa();
        add(name + "_kappa", kappa);

        auto& beta_variates = sbp.GetBetaVariates();
        add(name + "_betavariates", beta_variates);
    }

    template <class T>
    void add(std::string name, BranchArray<T>& v) {
        header_to_stream.push_back([&v, name](std::ostream& os) {
            int n = v.GetNbranch();
            if (n > 0) {
                os << name << "[0]";
                for (int i = 1; i < n; i++) os << "\t" << name << "[" << i << "]";
            }
        });
        data_to_stream.push_back([&v](std::ostream& os) {
            int n = v.GetNbranch();
            if (n > 0) {
                os << v.GetVal(0);
                for (int i = 1; i < n; i++) os << "\t" << v.GetVal(i);
            }
        });
        set_from_stream.push_back([&v](std::istream& is) {
            for (int i = 0; i < v.GetNbranch(); i++) is >> v[i];
        });
    }

    template <class T>
    void add(std::string name, Array<T>& v) {
        for(int i = 0; i < v.GetSize(); i++) {
            add(name + "[" + std::to_string(i) + "]", v[i]);
        }
    }

    template <class T>
    void add(std::string name, T* o, double (T::*f)() const) {
        header_to_stream.push_back([name](std::ostream& os) { os << name; });
        data_to_stream.push_back([o, f](std::ostream& os) { os << (o->*f)(); });
        set_from_stream.push_back([](std::istream& is) {
            double d;
            is >> d;
        });
    }

    void add(std::string name, std::function<double()> f) {
        header_to_stream.push_back([name](std::ostream& os) { os << name; });
        data_to_stream.push_back([f](std::ostream& os) { os << f(); });
        set_from_stream.push_back([](std::istream& is) {
            double d;
            is >> d;
        });
    }
};
