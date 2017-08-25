#ifndef HIGHRESTIMER_HH
#define HIGHRESTIMER_HH

#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

namespace HighResTimer {
using namespace std;

typedef std::chrono::high_resolution_clock::time_point time_point;

class HighResTimer {
  string name;
  atomic<int64_t> cnt;
  atomic<double> wsum, wsum2, wmin, wmax; // Nanoseconds
  atomic<double> bsum, bsum2, bmin, bmax; // Bytes

  void init();
  void add(time_point t0, time_point t1, int64_t nbytes);

  // Timers are registered, thus they can't be copied or moved
  HighResTimer() = delete;
  HighResTimer(const HighResTimer &) = delete;
  HighResTimer(HighResTimer &&) = delete;
  HighResTimer &operator=(const HighResTimer &) = delete;
  HighResTimer &operator=(HighResTimer &&) = delete;

public:
  HighResTimer(const string &name);
  ~HighResTimer();

  class HighResClock {
    HighResTimer *timer;
    time_point t0;
    friend class HighResTimer;
    HighResClock() = delete;
    HighResClock(HighResTimer *timer) : timer(timer) {
      t0 = std::chrono::high_resolution_clock::now();
    }

  public:
    void stop(int64_t nbytes) {
      time_point t1 = std::chrono::high_resolution_clock::now();
      timer->add(t0, t1, nbytes);
      timer = nullptr;
    }
  };
  friend class HighResClock;

  HighResClock start() __attribute__((__warn_unused_result__)) {
    return HighResClock(this);
  }

  double get_cnt() const; // count

  double get_wsum() const; // [ns]
  double get_wavg() const;
  double get_wsdv() const;
  double get_wmin() const;
  double get_wmax() const;

  double get_bsum() const; // [GB]
  double get_bavg() const;
  double get_bsdv() const;
  double get_bmin() const;
  double get_bmax() const;

  double get_bw() const; // [B/s]

  ostream &output(ostream &os) const;
  friend ostream &operator<<(ostream &os, const HighResTimer &timer) {
    return timer.output(os);
  }
};
}

#endif // #ifndef HIGHRESTIMER_HH
