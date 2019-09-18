#pragma once
#include <functional>
#include <iostream>
#include "model_decl_utils.hpp"
#include "mpi_components/partition.hpp"
#include "tags/decl_utils.hpp"
#include "traits.hpp"

// template <class T, class... Args>
// void call_interface_helper(std::true_type, T& target, Args&&... args) {
//     target.declare_interface(std::forward<Args>(args)...);
// }

// template <class T, class... Args>
// void call_interface_helper(std::false_type, T& target, Args&&... args) {
//     external_interface<T>::declare_interface(target, std::forward<Args>(args)...);
// }

// template <class T, class... Args>
// void call_interface(T& target, Args&&... args) {
//     call_interface_helper(has_interface<T>{}, target, std::forward<Args>(args)...);
// }


class Tracer {
    std::vector<std::function<void(std::ostream&)>> header_to_stream;
    std::vector<std::function<void(std::ostream&)>> data_to_stream;
    std::vector<std::function<void(std::istream&)>> set_from_stream;

  public:
    template <class Provider, class Test = processing::HasTag<ModelNode>>
    Tracer(Provider& p, Test test = processing::HasTag<ModelNode>()) {
        using namespace processing;
        using must_be_unrolled = Or<HasTrait<has_either_interface>, HasTrait<is_nontrivial_vector>>;
        using recursive_processing = RecursiveUnroll<must_be_unrolled, FullNameEnd>;
        using toplevel_filter = Filter<Test, recursive_processing>;
        auto prinfo = make_processing_info<toplevel_filter>(*this);
        call_interface(prinfo, p);
        // p.declare_interface(prinfo);
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

    size_t nbr_header_fields() const { return header_to_stream.size(); }

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

    void process_declaration(std::string name, double& d) {
        header_to_stream.emplace_back([name](std::ostream& os) { os << name; });
        data_to_stream.emplace_back([&d](std::ostream& os) { os << d; });
        set_from_stream.emplace_back([&d](std::istream& is) { is >> d; });
    }

    void process_declaration(std::string name, int& d) {
        header_to_stream.emplace_back([name](std::ostream& os) { os << name; });
        data_to_stream.emplace_back([&d](std::ostream& os) { os << d; });
        set_from_stream.emplace_back([&d](std::istream& is) { is >> d; });
    }

    template <class T>
    void process_declaration(
        std::string name, std::vector<T>& v, Partition partition = Partition(IndexSet(), 0)) {
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

    void process_declaration(std::string name, std::function<double()> const& f) {
        header_to_stream.emplace_back([name](std::ostream& os) { os << name; });
        data_to_stream.emplace_back([f](std::ostream& os) { os << f(); });
        set_from_stream.emplace_back([](std::istream& is) {
            double d;  // ignoring
            is >> d;
        });
    }
};
