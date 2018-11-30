#pragma once

#include "PartitionedBufferManager.hpp"
#include "partition.hpp"

/*==================================================================================================
  GatherMaster
==================================================================================================*/
class GatherMaster : public Proxy {
    Partition partition;
    PartitionedBufferManager manager;

  public:
    GatherMaster(Partition partition) : partition(partition), manager(partition) {}

    template <class T>
    void add(T& target) {
        static_assert(is_partitionable<T>::value, "GatherMaster: target is not partitionable!");
        assert(target.size() % partition.size_all() == 0);
        manager.add(target);
    }

    void acquire() final {
        assert(manager.buffer_size() > 0);
        auto buf = manager.receive_buffer();
        MPI_Gatherv(NULL, 0, MPI_PACKED, buf, manager.revcounts(), manager.displs(), MPI_PACKED, 0,
            MPI_COMM_WORLD);
        manager.receive();
    }
};

/*==================================================================================================
  GatherSlave
==================================================================================================*/
class GatherSlave : public Proxy {
    BufferManager manager;
    Partition partition;

  public:
    GatherSlave(Partition partition) : partition(partition) {}

    template <class T>
    void add(T& target) {
        static_assert(is_partitionable<T>::value, "GatherSlave: target is not partitionable!");
        assert(target.size() % partition.my_partition_size() == 0);
        manager.add(target);
    }

    void release() final {
        assert(manager.buffer_size() > 0);
        auto buf = manager.send_buffer();
        MPI_Gatherv(buf, manager.buffer_size(), MPI_PACKED, NULL, NULL, NULL, MPI_PACKED, 0,
            MPI_COMM_WORLD);
    }
};


/*==================================================================================================
  Gather functions
  Functions that are meant to be called globally and that will create either a master or slave
  component depending on the process
==================================================================================================*/
template <class... Variables>
std::unique_ptr<Proxy> gather(Partition partition, Variables&&... vars) {
    if (!MPI::p->rank) {  // master
        auto gatherer = new GatherMaster(partition);
        gatherer->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(gatherer));
    } else {  // slave
        auto gatherer = new GatherSlave(partition);
        gatherer->add(std::forward<Variables>(vars)...);
        return std::unique_ptr<Proxy>(dynamic_cast<Proxy*>(gatherer));
    }
}
