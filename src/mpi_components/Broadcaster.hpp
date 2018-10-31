#pragma once

#include <functional>
#include "Process.hpp"
#include "components/RegistrarBase.hpp"
#include "interfaces.hpp"

/*==================================================================================================
  Broadcaster
  An object responsible for broadcasting the values of specified fields from one process to others
==================================================================================================*/
template <typename T>
class BroadcasterMaster : public Proxy, public RegistrarBase<BroadcasterMaster<T>> {
    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<void()>> writers;


  public:
    BroadcasterMaster(std::set<std::string> filter, Process& p = *MPI::p)
        : RegistrarBase<BroadcasterMaster<T>>(filter), p(p), datatype(get_datatype<T>()) {}

    void register_element(std::string, T& target) {
        writers.push_back([&target, this]() { buf.push_back(target); });
    }

    void write_buffer() {
        buf.clear();
        for (auto writer : writers) { writer(); }
    }

    void release() final {
        write_buffer();
        MPI_Bcast(buf.data(), buf.size(), datatype, p.rank, MPI_COMM_WORLD);
    }
};

template <typename T>
class BroadcasterSlave : public Proxy, public RegistrarBase<BroadcasterSlave<T>> {
    using buf_it = typename std::vector<T>::iterator;

    Process& p;
    int origin;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<buf_it(buf_it)>> readers;

  public:
    BroadcasterSlave(std::set<std::string> filter, int origin = 0, Process& p = *MPI::p)
        : RegistrarBase<BroadcasterSlave<T>>(filter),
          p(p),
          origin(0),
          datatype(get_datatype<T>()) {}

    void register_element(std::string, T& target) {
        readers.push_back([&target](buf_it it) {
            target = *it;
            return it + 1;
        });
        buf.push_back(-1);
    }

    void read_buffer() {
        buf_it it = buf.begin();
        for (auto reader : readers) { it = reader(it); }
    }

    void acquire() final {
        MPI_Bcast(buf.data(), buf.size(), datatype, origin, MPI_COMM_WORLD);
        read_buffer();
    }
};

template <class Model, class T = double>
std::unique_ptr<Proxy> broadcast(Model& m,
    void (Model::*f_master)(RegistrarBase<BroadcasterMaster<T>>&),
    void (Model::*f_slave)(RegistrarBase<BroadcasterSlave<T>>&), std::set<std::string> filter) {
    std::unique_ptr<Proxy> result{nullptr};
    if (!MPI::p->rank) {
        auto component = new BroadcasterMaster<T>(filter);
        component->register_from_method(m, f_master);
        result.reset(dynamic_cast<Proxy*>(component));
    } else {
        auto component = new BroadcasterSlave<T>(filter);
        component->register_from_method(m, f_slave);
        result.reset(dynamic_cast<Proxy*>(component));
    }
    return result;
}

template <class Model, class T = double>
std::unique_ptr<Proxy> broadcast_model(Model& m, std::set<std::string> filter) {
    return broadcast(m, &Model::declare_model, &Model::declare_model, filter);
}