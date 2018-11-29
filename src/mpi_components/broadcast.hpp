#pragma once

#include <functional>
#include "Process.hpp"
#include "components/RegistrarBase.hpp"
#include "serialization.hpp"
#include "utils.hpp"

/*==================================================================================================
  BroadcasterMaster
  An object responsible for broadcasting the values of specified fields to other processes
  is meant to communicate with one BroadcasterSlave per other process
==================================================================================================*/
class BroadcasterMaster : public Proxy, public RegistrarBase<BroadcasterMaster> {
    BufferManager buf;

  public:
    BroadcasterMaster() { assert(MPI::p->rank == 0); }

    using RegistrarBase<BroadcasterMaster>::register_element;

    template <class T>
    void register_element(std::string, T& target) {
        buf.add(target);
    }

    void release() final {
        assert(buf.buffer_size() > 0);
        auto ready_buf = buf.send_buffer();
        MPI_Bcast(ready_buf, buf.buffer_size(), MPI_PACKED, 0, MPI_COMM_WORLD);
    }
};

/*==================================================================================================
  BroadcasterSlave
==================================================================================================*/
class BroadcasterSlave : public Proxy, public RegistrarBase<BroadcasterSlave> {
    BufferManager buf;

  public:
    BroadcasterSlave() {}

    using RegistrarBase<BroadcasterSlave>::register_element;

    template <class T>
    void register_element(std::string, T& target) {
        buf.add(target);
    }

    void acquire() final {
        assert(buf.buffer_size() > 0);
        auto rcv_buf = buf.receive_buffer();
        MPI_Bcast(rcv_buf, buf.buffer_size(), MPI_PACKED, 0, MPI_COMM_WORLD);
        buf.receive();
    }
};

/*==================================================================================================
  Broadcast functions
  Functions that are meant to be called globally and that will create either a master or slave
  component depending on the process
==================================================================================================*/
template <class T, class Model, class... Args>
std::unique_ptr<Proxy> broadcast(Model& m, filter_t filter) {
    return instantiate_and_declare_from_model<BroadcasterMaster, BroadcasterSlave, Model>(
        m, filter);
}
