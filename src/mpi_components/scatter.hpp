#pragma once

#include "PartitionedBufferManager.hpp"
#include "partition.hpp"

/*==================================================================================================
  ScatterMaster
==================================================================================================*/
class ScatterMaster : public Proxy {
    Partition partition;
    PartitionedBufferManager manager;

  public:
    ScatterMaster(Partition partition) : partition(partition), manager(partition) {}

    template <class... Variables>
    void add(Variables&&... vars) {
        manager.add(std::forward<Variables>(vars)...);
    }

    void release() final {
        assert(manager.buffer_size() > 0);
        auto buf = manager.send_buffer();
        MPI_Scatterv(buf, manager.revcounts(), manager.displs(), MPI_PACKED, nullptr, 0, MPI_PACKED,
            0, MPI_COMM_WORLD);
    }
};

/*==================================================================================================
  ScatterSlave
==================================================================================================*/
class ScatterSlave : public Proxy {
    BufferManager manager;
    Partition partition;

  public:
    ScatterSlave(Partition partition) : partition(partition) {}

    template <class... Variables>
    void add(Variables&&... vars) {
        manager.add(std::forward<Variables>(vars)...);
    }

    void acquire() final {
        assert(manager.buffer_size() > 0);
        auto buf = manager.receive_buffer();
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_PACKED, buf, manager.buffer_size(), MPI_PACKED,
            0, MPI_COMM_WORLD);
        manager.receive();
    }
};


/*==================================================================================================
  Scatter functions
  Functions that are meant to be called globally and that will create either a master or slave
  component depending on the process
==================================================================================================*/
template <class... Variables>
std::unique_ptr<Proxy> scatter(Partition partition, Variables&&... vars) {
    if (!MPI::p->rank) {  // master
        auto scatterer = new ScatterMaster(partition);
        scatterer->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(scatterer));
    } else {  // slave
        auto scatterer = new ScatterSlave(partition);
        scatterer->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(scatterer));
    }
}
