#pragma once

#include <iostream>
#include <map>
#include <memory>

class AbstractMonitor {
  public:
    virtual void print(std::ostream& os) const = 0;
    virtual ~AbstractMonitor() = default;
};

class MonitorManager {
    std::map<std::string, std::unique_ptr<AbstractMonitor>> monitors;

  public:
    MonitorManager() = default;

    template <class T, class... Args>
    void new_monitor(std::string name, Args&&... args) {
        static_assert(std::is_base_of<AbstractMonitor, T>::value,
            "Monitor does not inherit from AbstractMonitor");
        monitors.emplace(name, dynamic_cast<AbstractMonitor*>(new T(std::forward<Args>(args)...)));
    }

    void print(std::ostream& os) const {
        for (auto& monitor : monitors) {
            os << monitor.first << ": ";
            monitor.second->print(os);
            os << "\n";
        }
    }

    template <class M, class F, class... Args>
    void run_and_monitor(std::string name, F f, Args&&... args) {
        assert(dynamic_cast<M*>(monitors.at(name).get()) != nullptr);
        auto result = f(std::forward<Args>(args)...);
        auto& monitor_ref = dynamic_cast<M&>(*monitors.at(name));
        monitor_ref.update(result);
    }
};

// global monitor
extern std::unique_ptr<MonitorManager> gm;
