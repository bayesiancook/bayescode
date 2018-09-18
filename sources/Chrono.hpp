#ifndef CHRONO_H
#define CHRONO_H

#include <sys/time.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>

/**
 * \brief A chronometer
 *
 * Measures time in milliseconds; also implements a counter, allowing to both
 * measure total time and count e.g. number of cycles (and then return mean time
 * per cycle).
 */

class Chrono {
  public:
    Chrono() { Reset(); }
    //! reset chrono
    void Reset();
    //! start chrono
    void Start();
    //! stop chrono
    void Stop();

    //! increment the counter of cycles
    inline int operator++() { return N++; }

    //! get total time
    inline double GetTime() const { return TotalTime; }

    //! get time per cycle
    inline double GetTimePerCount() const { return TotalTime / N; }

    //! get total number of cycles
    inline int GetCount() const { return N; }

  private:
    // this is in milli seconds
    double sec1;
    double sec2;
    double milli1;
    double milli2;
    double TotalTime;
    int N;
};

class MeasureTime : public std::stringstream {
  public:
    MeasureTime() : stopped(false) { start(); }

    void start() {
        counter = std::chrono::high_resolution_clock::now();
        stopped = false;
    }

    void stop() {
        stopped = true;
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - counter);
    }

    template <int i>
    void print(std::string message) {
        if (!stopped) { stop(); }
        std::string left(2 * i, ' ');
        std::cout << left << "* " << message << str() << "Time: " << duration.count()  //  << "ms."
                  << std::endl;
        str("");
        start();
    }

    template <int i>
    void print() {
        print<i>("");
    }

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> counter;
    std::chrono::milliseconds duration;
    bool stopped;
    std::stringstream ss;
};

#endif  // CHRONO_H
