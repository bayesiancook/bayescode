#pragma once
#include <iostream>
#include "BidimArray.hpp"
#include "BranchArray.hpp"
#include "MultinomialAllocationVector.hpp"
#include "StickBreakingProcess.hpp"
#include "components/RegistrarBase.hpp"
#include "mpi_components/partition.hpp"

class Tracer : public RegistrarBase<Tracer> {
    std::vector<std::function<void(std::ostream&)>> header_to_stream;
    std::vector<std::function<void(std::ostream&)>> data_to_stream;
    std::vector<std::function<void(std::istream&)>> set_from_stream;

    friend RegistrarBase<Tracer>;
    using RegistrarBase<Tracer>::register_element;

  public:
    template <class T>
    Tracer(T& x, void (T::*f)(Tracer&)) {
        (x.*f)(*this);
    }

    void write_header(std::ostream& os) const {
        size_t n = header_to_stream.size();
        if (n > 0) {
            header_to_stream.at(0)(os);
            for (size_t i = 1; i < n; i++) {
                os << "\t";
                header_to_stream.at(i)(os);
            }
        }
    }

    unsigned nbr_header_fields() const { return static_cast<unsigned>(header_to_stream.size()); }

    std::vector<double> line_values() const {
        std::stringstream ss_line;
        write_line(ss_line);
        std::vector<double> values{};
        std::string str_value;
        while (getline(ss_line, str_value, '\t')) { values.push_back(stod(str_value)); }
        return values;
    }

    void write_line(std::ostream& os) const {
        size_t n = data_to_stream.size();
        if (n > 0) {
            os << "\n";
            data_to_stream.at(0)(os);
            for (size_t i = 1; i < n; i++) {
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

    void register_element(std::string const& name, double& d) {
        header_to_stream.emplace_back([name](std::ostream& os) { os << name; });
        data_to_stream.emplace_back([&d](std::ostream& os) { os << d; });
        set_from_stream.emplace_back([&d](std::istream& is) { is >> d; });
    }

    void register_element(std::string const& name, int& d) {
        header_to_stream.emplace_back([name](std::ostream& os) { os << name; });
        data_to_stream.emplace_back([&d](std::ostream& os) { os << d; });
        set_from_stream.emplace_back([&d](std::istream& is) { is >> d; });
    }

    template <class T>
    void register_element(std::string const& name, std::vector<T>& v,
        Partition partition = Partition(IndexSet(), 0)) {
        /* -- */
        header_to_stream.emplace_back([&v, name](std::ostream& os) {
            size_t n = v.size();
            if (n > 0) {
                os << name << "[0]";
                for (size_t i = 1; i < n; i++) os << "\t" << name << "[" << i << "]";
            }
        });
        data_to_stream.emplace_back([&v](std::ostream& os) {
            size_t n = v.size();
            if (n > 0) {
                os << v.at(0);
                for (size_t i = 1; i < n; i++) os << "\t" << v.at(i);
            }
        });
        set_from_stream.emplace_back([&v](std::istream& is) {
            for (auto& e : v) is >> e;
        });
    }

    void register_element(std::string const& name, StickBreakingProcess& sbp) {
        add(name + "_array", dynamic_cast<SimpleArray<double>&>(sbp));

        auto& beta_variates = sbp.GetBetaVariates();
        add(name + "_betavariates", beta_variates);
    }

    template <class T>
    void register_element(std::string const& name, BranchArray<T>& v) {
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
    void register_element(std::string const& name, Array<T>& v) {
        for (int i = 0; i < v.GetSize(); i++) { add(name + "[" + std::to_string(i) + "]", v[i]); }
    }

    template <class T>
    void register_element(std::string const& name, T* o, double (T::*f)() const) {
        header_to_stream.push_back([name](std::ostream& os) { os << name; });
        data_to_stream.push_back([o, f](std::ostream& os) { os << (o->*f)(); });
        set_from_stream.push_back([](std::istream& is) {
            double d;
            is >> d;
        });
    }

    void register_element(std::string const& name, std::function<double()> const& f) {
        header_to_stream.emplace_back([name](std::ostream& os) { os << name; });
        data_to_stream.emplace_back([f](std::ostream& os) { os << f(); });
        set_from_stream.emplace_back([](std::istream& is) {
            double d;
            is >> d;
        });
    }
};
