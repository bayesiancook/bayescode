#pragma once

#include "BufferManager.hpp"
#include "partition.hpp"

class PartitionedBufferManager {
    Partition partition;
    std::vector<BufferManager> temp_managers;

    BufferManager manager;
    std::vector<int> _revcounts, _displs;
    bool manager_ready{false};

    void check_manager() {
        if (!manager_ready) {
            for (size_t subset = partition.first(); subset < partition.max_index(); subset++) {
                manager.merge(temp_managers.at(subset));
                _revcounts.at(subset) = temp_managers.at(subset).buffer_size();
                for (size_t i = subset + 1; i < partition.max_index(); i++) {
                    _displs.at(i) += temp_managers.at(subset).buffer_size();
                }
            }
            manager_ready = true;
        }
    }

  public:
    PartitionedBufferManager(Partition partition)
        : partition(partition),
          temp_managers(partition.max_index()),
          _revcounts(partition.max_index()),
          _displs(partition.max_index(), 0) {}

    template <class T>
    void add(T& x) {
        custom_add_dispatch(x, has_custom_serialization<T>());
    }

  private:
    template <class T>
    void custom_add_dispatch(T& x, std::false_type) {
        static_assert(is_partitionable<T>::value,
            "PartitionedBufferManager::add: type T is not partitionable");
        assert(x.size() % partition.size_all() == 0);
        assert(x.size() != 0);
        assert(!manager_ready);

        size_t multiplicity = x.size() / partition.size_all();
        int i = 0;
        for (size_t subset = partition.first(); subset < partition.max_index(); subset++) {
            size_t subset_size = partition.partition_size(subset);
            size_t nb_elements = subset_size * multiplicity;
            temp_managers.at(subset).add_subset(x, i, nb_elements);
            i += nb_elements;
        }
    }

    template <class T>
    void custom_add_dispatch(T& x, std::true_type) {
        x.serialization_interface(*this);
    }

  public:
    template <class Var, class... Vars>
    void add(Var& var, Vars&&... vars) {
        add(var);
        add(std::forward<Vars>(vars)...);
    }

    void* send_buffer() {
        check_manager();
        return manager.send_buffer();
    }

    void* receive_buffer() {
        check_manager();
        return manager.receive_buffer();
    }

    void receive() {
        check_manager();
        manager.receive();
    }

    size_t buffer_size() {
        check_manager();
        return manager.buffer_size();
    }

    int* revcounts() {
        check_manager();
        return _revcounts.data();
    }

    int* displs() {
        check_manager();
        return _displs.data();
    }
};
