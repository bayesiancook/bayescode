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
    MonitorManager() = default;

    template <class T, class... Args>
    void new_monitor(Args&&... args) {
        static_assert(std::is_base_of<AbstractMonitor, T>::value, "YOLO");
        monitors.emplace_back(dynamic_cast<AbstractMonitor*>(new T(std::forward<Args>(args)...)));
    }

    void print(std::ostream& os) const {
        for (auto& monitor : monitors) {
            monitor->print(os);
            os << "\n";
        }
    }
};