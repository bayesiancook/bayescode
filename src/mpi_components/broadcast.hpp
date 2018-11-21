#pragma once

#include <functional>
#include "Process.hpp"
#include "components/RegistrarBase.hpp"
#include "interfaces.hpp"

/*==================================================================================================
  BroadcasterMaster
  An object responsible for broadcasting the values of specified fields to other processes
  is meant to communicate with one BroadcasterSlave per other process
==================================================================================================*/
template <typename T>
class BroadcasterMaster : public Proxy, public RegistrarBase<BroadcasterMaster<T>> {
    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<void()>> writers;

    friend RegistrarBase<BroadcasterMaster<T>>;
    using RegistrarBase<BroadcasterMaster<T>>::register_element;

    void register_element(std::string, T& target) {
        writers.push_back([&target, this]() { buf.push_back(target); });
    }

    void register_element(std::string, std::vector<T>& target) {
        writers.push_back(
            [&target, this]() { buf.insert(buf.end(), target.begin(), target.end()); });
    }

    void write_buffer() {
        buf.clear();
        for (auto writer : writers) { writer(); }
    }

  public:
    BroadcasterMaster(Process& p = *MPI::p) : p(p), datatype(get_datatype<T>()) {}

    void release() final {
        write_buffer();
        MPI_Bcast(buf.data(), buf.size(), datatype, p.rank, MPI_COMM_WORLD);
    }
};

/*==================================================================================================
  BroadcasterSlave
==================================================================================================*/
template <typename T>
class BroadcasterSlave : public Proxy, public RegistrarBase<BroadcasterSlave<T>> {
    using buf_it = typename std::vector<T>::iterator;

    Process& p;
    int origin;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<void(buf_it&)>> readers;

    friend RegistrarBase<BroadcasterSlave<T>>;
    using RegistrarBase<BroadcasterSlave<T>>::register_element;

    void register_element(std::string, T& target) {
        readers.push_back([&target](buf_it& it) {
            target = *it;
            it++;
        });
        buf.emplace_back();
    }

    void register_element(std::string, std::vector<T>& target) {
        readers.push_back([&target](buf_it& it) {
            target = std::vector<T>(it, it + target.size());
            it += target.size();
        });
        buf.insert(buf.end(), std::vector<T>(target.size()));
    }

    void read_buffer() {
        buf_it it = buf.begin();
        for (auto reader : readers) { reader(it); }
    }

  public:
    BroadcasterSlave(int origin = 0, Process& p = *MPI::p)
        : p(p), origin(0), datatype(get_datatype<T>()) {}

    void acquire() final {
        MPI_Bcast(buf.data(), buf.size(), datatype, origin, MPI_COMM_WORLD);
        read_buffer();
    }
};

/*==================================================================================================
  Broadcast functions
  Functions that are meant to be called globally and that will create either a master or slave
  component depending on the process
==================================================================================================*/
template <class Model, class T = double>
std::unique_ptr<Proxy> broadcast(Model& m,
    void (Model::*f_master)(RegistrarBase<BroadcasterMaster<T>>&),
    void (Model::*f_slave)(RegistrarBase<BroadcasterSlave<T>>&),
    std::set<std::string> filter = std::set<std::string>{}) {
    std::unique_ptr<Proxy> result{nullptr};
    if (!MPI::p->rank) {
        auto component = new BroadcasterMaster<T>();
        component->register_from_method(m, f_master, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    } else {
        auto component = new BroadcasterSlave<T>();
        component->register_from_method(m, f_slave, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    }
    return result;
}

template <class Model, class T = double>
std::unique_ptr<Proxy> broadcast_model(
    Model& m, std::set<std::string> filter = std::set<std::string>{}) {
    return broadcast<Model, T>(m, &Model::declare_model, &Model::declare_model, filter);
}
