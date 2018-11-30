#pragma once

#include <functional>
#include "BufferManager.hpp"
#include "Process.hpp"

/*==================================================================================================
  BroadcasterMaster
  An object responsible for broadcasting the values of specified fields to other processes
  is meant to communicate with one BroadcasterSlave per other process
==================================================================================================*/
class BroadcasterMaster : public Proxy {
    BufferManager buf;

  public:
    BroadcasterMaster() { assert(MPI::p->rank == 0); }

    template <class... Variables>
    void add(Variables&&... vars) {
        buf.add(std::forward<Variables>(vars)...);
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
class BroadcasterSlave : public Proxy {
    BufferManager buf;

  public:
    BroadcasterSlave() {}

    template <class... Variables>
    void add(Variables&&... vars) {
        buf.add(std::forward<Variables>(vars)...);
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
template <class... Variables>
std::unique_ptr<Proxy> broadcast(Variables&&... vars) {
    if (!MPI::p->rank) {  // master
        auto broadcaster = new BroadcasterMaster();
        broadcaster->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(broadcaster));
    } else {  // slave
        auto broadcaster = new BroadcasterSlave();
        broadcaster->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(broadcaster));
    }
}