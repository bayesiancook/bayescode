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

    void acquire() final {
        MPI_Request request_int;
        MPI_Request request_double;
        manager.receive_buffer();  // so the buffer is ready
        if (manager.buffer_int_size() > 0) {
            if (zeroes_int.size() != manager.buffer_int_size()) {
                zeroes_int = std::vector<int>(manager.buffer_int_size(), 0);
            }
            auto buf_int = manager.receive_int_buffer();
            MPI_Ireduce(zeroes_int.data(), buf_int, manager.buffer_int_size(), MPI_INT, MPI_SUM,
                MPI::p->rank, MPI_COMM_WORLD, &request_int);
        }
        if (manager.buffer_double_size() > 0) {
            if (zeroes_double.size() != manager.buffer_double_size()) {
                zeroes_double = std::vector<double>(manager.buffer_double_size(), 0);
            }
            auto buf_double = manager.receive_double_buffer();
            MPI_Ireduce(zeroes_double.data(), buf_double, manager.buffer_double_size(), MPI_DOUBLE,
                MPI_SUM, MPI::p->rank, MPI_COMM_WORLD, &request_double);
        }
        if (manager.buffer_int_size() > 0) { MPI_Wait(&request_int, MPI_STATUS_IGNORE); }
        if (manager.buffer_double_size() > 0) { MPI_Wait(&request_double, MPI_STATUS_IGNORE); }
        manager.receive();
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

    void release() final {
        MPI_Request request_int;
        MPI_Request request_double;
        manager.send_buffer();  // so the buffer is ready
        if (manager.buffer_int_size() > 0) {
            auto buf_int = manager.send_int_buffer();
            MPI_Ireduce(buf_int, NULL, manager.buffer_int_size(), MPI_INT, MPI_SUM, 0,
                MPI_COMM_WORLD, &request_int);
        }
        if (manager.buffer_double_size() > 0) {
            auto buf_double = manager.send_double_buffer();
            MPI_Ireduce(buf_double, NULL, manager.buffer_double_size(), MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD, &request_double);
        }
        if (manager.buffer_int_size() > 0) { MPI_Wait(&request_int, MPI_STATUS_IGNORE); }
        if (manager.buffer_double_size() > 0) { MPI_Wait(&request_double, MPI_STATUS_IGNORE); }
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