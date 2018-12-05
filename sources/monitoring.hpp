#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <sstream>

class AbstractMonitor {
  public:
    virtual void print(std::ostream& os) const = 0;
    virtual ~AbstractMonitor() = default;
};

class MonitorManager {
  public:
    std::map<std::string, std::unique_ptr<AbstractMonitor>> monitors;

    MonitorManager() = default;

    template <class T, class... Args>
    void new_monitor(std::string name, Args&&... args) {
        static_assert(std::is_base_of<AbstractMonitor, T>::value,
                      "Monitor does not inherit from AbstractMonitor");
        monitors.emplace(std::piecewise_construct, std::forward_as_tuple(name),
                         std::forward_as_tuple(
                             dynamic_cast<AbstractMonitor*>(new T(std::forward<Args>(args)...))));
    }

    void print(std::ostream& os) const {
        for (auto& monitor : monitors) {
            os << monitor.first << ": ";
            monitor.second->print(os);
            os << "\n";
        }
    }

    template <class M, class F, class... Args>
    void run_and_monitor(std::string name, F&& f, Args&&... args) {
        if (monitors.find(name) == monitors.end()) {
            new_monitor<M>(name);  // HACKISH (no constructor params)
        }
        auto result = f(std::forward<Args>(args)...);
        auto& monitor_ref = dynamic_cast<M&>(*monitors.at(name));
        monitor_ref.update(result);
    }
};

// global monitor
extern std::unique_ptr<MonitorManager> gm;

template <class T>
class MeanMonitor : public AbstractMonitor {
    T sum{0};
    int count{0};

    T tmp_sum{0};
    int tmp_count{0};

  public:
    void print(std::ostream& os) const final { os << sum / static_cast<double>(count); }

    void update(T x) {
        sum += x;
        tmp_sum += x;
        count++;
        tmp_count++;
    }

    void tmp_reset() {
        tmp_sum = 0;
        tmp_count = 0;
    }

    double tmp_mean() const {
        return tmp_sum / tmp_count;
    }
};
