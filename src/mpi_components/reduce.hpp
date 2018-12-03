#pragma once

#include <functional>
#include "Process.hpp"
#include "components/RegistrarBase.hpp"
#include "utils.hpp"

/*==================================================================================================
  ReducerMaster
  An object responsible for Reducing the values of specified fields to other processes
  is meant to communicate with one ReducerSlave per other process
==================================================================================================*/
class ReducerMaster : public Proxy {
    BufferManager manager;

    std::vector<int> zeroes_int;
    std::vector<double> zeroes_double;

  public:
    ReducerMaster() = default;

    template <class... Variables>
    void add(Variables&&... vars) {
        manager.add(std::forward<Variables>(vars)...);
    }

    void reduce_ints() {
        if (manager.nb_ints() > 0) {
            if (zeroes_int.size() != manager.nb_ints()) {
                zeroes_int = std::vector<int>(manager.nb_ints(), 0);
            }
            BufferManager tmp_manager = manager.int_manager();
            auto buf = tmp_manager.receive_buffer();
            auto nb_elems = tmp_manager.nb_ints();
            MPI_Reduce(zeroes_int.data(), buf, nb_elems, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            tmp_manager.receive();
        }
    }

    void reduce_doubles() {
        if (manager.nb_doubles() > 0) {
            if (zeroes_double.size() != manager.nb_doubles()) {
                zeroes_double = std::vector<double>(manager.nb_doubles(), 0);
            }
            BufferManager tmp_manager = manager.double_manager();
            auto buf = tmp_manager.receive_buffer();
            auto nb_elems = tmp_manager.nb_doubles();
            MPI_Reduce(zeroes_double.data(), buf, nb_elems, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            tmp_manager.receive();
        }
    }

    void acquire() final {
        assert(manager.buffer_size() > 0);
        reduce_ints();
        reduce_doubles();
    }
};

/*==================================================================================================
  ReducerSlave
==================================================================================================*/
class ReducerSlave : public Proxy {
    BufferManager manager;

  public:
    ReducerSlave() = default;

    template <class... Variables>
    void add(Variables&&... vars) {
        manager.add(std::forward<Variables>(vars)...);
    }

    void reduce_ints() {
        if (manager.nb_ints() > 0) {
            BufferManager tmp_manager = manager.int_manager();
            auto buf = tmp_manager.send_buffer();
            auto nb_elems = tmp_manager.nb_ints();
            MPI_Reduce(buf, nullptr, nb_elems, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }

    void reduce_doubles() {
        if (manager.nb_doubles() > 0) {
            BufferManager tmp_manager = manager.double_manager();
            auto buf = tmp_manager.send_buffer();
            auto nb_elems = tmp_manager.nb_doubles();
            MPI_Reduce(buf, nullptr, nb_elems, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }

    void release() final {
        assert(manager.buffer_size() > 0);
        reduce_ints();
        reduce_doubles();
    }
};

/*==================================================================================================
  Reduce functions
  Functions that are meant to be called globally and that will create either a master or slave
  component depending on the process
==================================================================================================*/
template <class... Variables>
std::unique_ptr<Proxy> reduce(Variables&&... vars) {
    if (!MPI::p->rank) {  // master
        auto reducer = new ReducerMaster();
        reducer->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(reducer));
    } else {  // slave
        auto reducer = new ReducerSlave();
        reducer->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(reducer));
    }
}