#pragma once

#include <iostream>
#include <memory>
#include <vector>

class AbstractMonitor {
  public:
    virtual void print(std::ostream& os) const = 0;
    virtual ~AbstractMonitor() = default;
};

class MonitorManager {
    std::vector<std::unique_ptr<AbstractMonitor>> monitors;

  public:
    using monitor_index_t = size_t;

    MonitorManager() = default;

    template <class T, class... Args>
    monitor_index_t new_monitor(Args&&... args) {
        static_assert(std::is_base_of<AbstractMonitor, T>::value,
            "Monitor does not inherit from AbstractMonitor");
        monitors.emplace_back(dynamic_cast<AbstractMonitor*>(new T(std::forward<Args>(args)...)));
        return monitors.size() - 1;
    }

    void print(std::ostream& os) const {
        for (auto& monitor : monitors) {
            monitor->print(os);
            os << "\n";
        }
    }

    template <class M, class F, class... Args>
    void run_and_monitor(monitor_index_t i, F f, Args&&... args) {
        assert(dynamic_cast<M*>(monitors.at(i).get()) != nullptr);
        auto result = f(std::forward<Args>(args)...);
        auto& monitor_ref = dynamic_cast<M&>(*monitors.at(i));
        monitor_ref.update(result);
    }
};