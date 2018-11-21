#pragma once

#include "components/RegistrarBase.hpp"
#include "interfaces.hpp"
#include "partition.hpp"

/*==================================================================================================
  GatherMaster
==================================================================================================*/
template <typename T>
class GatherMaster : public Proxy, public RegistrarBase<GatherMaster<T>> {
    using buf_it = typename std::vector<T>::iterator;

    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<int> displs, revcounts;
    std::vector<std::function<void(std::vector<buf_it>&)>> readers;

    friend RegistrarBase<GatherMaster<T>>;
    using RegistrarBase<GatherMaster<T>>::register_element;

    void register_element(std::string, std::vector<T>& target, const Partition& partition) {
        for (int i = 1; i < p.size; i++) {
            size_t p_size = partition.partition_size(i);
            revcounts.at(i) += p_size;
            for (int j = i + 1; j < p.size; j++) { displs.at(j) += p_size; }
        }
        for (size_t i = 0; i < partition.size_all(); i++) { buf.push_back(-1); }

        // FIXME not very efficient
        readers.push_back([&target, partition, this](std::vector<buf_it>& its) {
            for (int i = 1; i < p.size; i++) {
                buf_it end = its.at(i) + partition.partition_size(i);
                std::vector<T> tmp(its.at(i), end);
                target.insert(target.end(), tmp.begin(), tmp.end());
                its.at(i) = end;
            }
        });
    }

    void read_buffer() {
        std::vector<buf_it> its{buf.begin()};
        for (int i = 1; i < p.size; i++) { its.push_back(buf.begin() + displs.at(i)); }
        for (auto reader : readers) { reader(its); }
    }

  public:
    GatherMaster(Process& p = *MPI::p)
        : p(p), datatype(get_datatype<T>()), displs(p.size, 0), revcounts(p.size, 0) {}

    void acquire() final {
        MPI_Gatherv(NULL, 0, datatype, buf.data(), revcounts.data(), displs.data(), datatype, 0,
            MPI_COMM_WORLD);
        read_buffer();
    }
};

/*==================================================================================================
  GatherSlave
==================================================================================================*/
template <typename T>
class GatherSlave : public Proxy, public RegistrarBase<GatherSlave<T>> {
    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<void()>> writers;

    friend RegistrarBase<GatherSlave<T>>;
    using RegistrarBase<GatherSlave<T>>::register_element;

    void register_element(std::string, std::vector<T>& target, Partition&) {
        writers.push_back(
            [&target, this]() { buf.insert(buf.end(), target.begin(), target.end()); });
    }

    void write_buffer() {
        buf.clear();
        for (auto writer : writers) { writer(); }
    }

  public:
    GatherSlave(Process& p = *MPI::p) : p(p), datatype(get_datatype<T>()) {}

    void release() final {
        write_buffer();
        MPI_Gatherv(
            buf.data(), buf.size(), datatype, NULL, NULL, NULL, datatype, 0, MPI_COMM_WORLD);
    }
};


/*==================================================================================================
  Gather functions
  Functions that are meant to be called globally and that will create either a master or slave
  component depending on the process
==================================================================================================*/
template <class Model, class T = double>
std::unique_ptr<Proxy> gather(Model& m, void (Model::*f_master)(RegistrarBase<GatherMaster<T>>&),
    void (Model::*f_slave)(RegistrarBase<GatherSlave<T>>&),
    std::set<std::string> filter = std::set<std::string>{}) {
    std::unique_ptr<Proxy> result{nullptr};
    if (!MPI::p->rank) {
        auto component = new GatherMaster<T>();
        component->register_from_method(m, f_master, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    } else {
        auto component = new GatherSlave<T>();
        component->register_from_method(m, f_slave, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    }
    return result;
}

template <class Model, class T = double>
std::unique_ptr<Proxy> gather_model(
    Model& m, std::set<std::string> filter = std::set<std::string>{}) {
    return gather<Model, T>(m, &Model::declare_model, &Model::declare_model, filter);
}
