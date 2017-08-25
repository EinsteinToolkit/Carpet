#include "HighResTimer.hh"

#include <cctk.h>
#include <cctk_Arguments.h>

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace HighResTimer {
using namespace std;

struct HighResTimerSet {
  map<string, HighResTimer *> timers;

  void add(const string &timername, HighResTimer *timer);
  void remove(const string &timername);

  void output(ostream &os) const;
};
HighResTimerSet timer_set;

void HighResTimerSet::add(const string &timername, HighResTimer *timer) {
  if (timers.count(timername) > 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "HighResTimer \"%s\" already registered",
               timername.c_str());
  timers[timername] = timer;
}

void HighResTimerSet::remove(const string &timername) {
  auto count = timers.erase(timername);
  if (count == 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "HighResTimer \"%s\" not registered",
               timername.c_str());
}

void HighResTimerSet::output(ostream &os) const {
  for (const auto &tnt : timers) {
    // const auto &timername = tnt.first;
    const auto &timer = tnt.second;
    os << "    " << *timer << "\n";
  }
}

extern "C" void HighResTimer_OutputAllTimers(CCTK_ARGUMENTS) {
  CCTK_INFO("High resolution timers:");
  timer_set.output(cout);
}

////////////////////////////////////////////////////////////////////////////////

HighResTimer::HighResTimer(const string &name) : name(name) {
  init();
  timer_set.add(name, this);
}

HighResTimer::~HighResTimer() { timer_set.remove(name); }

void HighResTimer::init() {
  cnt = 0;
  wsum = 0.0;
  wsum2 = 0.0;
  wmin = HUGE_VAL;
  wmax = 0.0;
  bsum = 0.0;
  bsum2 = 0.0;
  bmin = HUGE_VAL;
  bmax = 0.0;
}

template <typename F, typename T> T atomic_update(atomic<T> &var, const F &f) {
  T oldval(var);
  while (!var.compare_exchange_weak(oldval, f(oldval), memory_order_acq_rel))
    ;
  return oldval;
}

void HighResTimer::add(time_point t0, time_point t1, int64_t nbytes0) {
  double wtime =
      std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
  double nbytes = nbytes0;
  ++cnt;
  atomic_update(wsum, [&](double wsum) { return wsum + wtime; });
  atomic_update(wsum2, [&](double wsum2) { return wsum2 + pow(wtime, 2); });
  atomic_update(wmin, [&](double wmin) { return fmin(wmin, wtime); });
  atomic_update(wmax, [&](double wmax) { return fmax(wmax, wtime); });
  atomic_update(bsum, [&](double bsum) { return bsum + nbytes; });
  atomic_update(bsum2, [&](double bsum2) { return bsum2 + pow(nbytes, 2); });
  atomic_update(bmin, [&](double bmin) { return fmin(bmin, nbytes); });
  atomic_update(bmax, [&](double bmax) { return fmax(bmax, nbytes); });
}

double HighResTimer::get_cnt() const { return cnt; }

double HighResTimer::get_wsum() const { return wsum; }
double HighResTimer::get_wavg() const { return cnt == 0 ? 0.0 : wsum / cnt; }
double HighResTimer::get_wsdv() const {
  return cnt == 0 ? 0.0 : sqrt(fdim(wsum2 / cnt, pow(wsum / cnt, 2)));
}
double HighResTimer::get_wmin() const { return cnt == 0 ? 0.0 : double(wmin); }
double HighResTimer::get_wmax() const { return wmax; }

double HighResTimer::get_bsum() const { return bsum; }
double HighResTimer::get_bavg() const { return cnt == 0 ? 0.0 : bsum / cnt; }
double HighResTimer::get_bsdv() const {
  return cnt == 0 ? 0.0 : sqrt(fdim(bsum2 / cnt, pow(bsum / cnt, 2)));
}
double HighResTimer::get_bmin() const { return cnt == 0 ? 0.0 : double(bmin); }
double HighResTimer::get_bmax() const { return bmax; }

double HighResTimer::get_bw() const { return cnt == 0 ? 0.0 : bsum / wsum; }

ostream &HighResTimer::output(ostream &os) const {
  auto oldprec = os.precision(3);
  os << name << ": cnt=" << llrint(get_cnt())
     << " wsum=" << llrint(get_wsum() * 1.0e-6) * 1.0e-3
     << " wavg=" << llrint(get_wavg() * 1.0e-6) * 1.0e-3
     << " wsdv=" << llrint(get_wsdv() * 1.0e-6) * 1.0e-3
     << " wmin=" << llrint(get_wmin() * 1.0e-6) * 1.0e-3
     << " wmax=" << llrint(get_wmax() * 1.0e-6) * 1.0e-3 << " [s]"
     << " bsum=" << llrint(get_bsum() * 1.0e-6) * 1.0e-3
     << " bavg=" << llrint(get_bavg() * 1.0e-6) * 1.0e-3
     << " bsdv=" << llrint(get_bsdv() * 1.0e-6) * 1.0e-3
     << " bmin=" << llrint(get_bmin() * 1.0e-6) * 1.0e-3
     << " bmax=" << llrint(get_bmax() * 1.0e-6) * 1.0e-3 << " [GB]"
     << " wsum=" << llrint(get_wsum()) << " wavg=" << llrint(get_wavg())
     << " wsdv=" << llrint(get_wsdv()) << " wmin=" << llrint(get_wmin())
     << " wmax=" << llrint(get_wmax()) << " [ns]"
     << " bsum=" << llrint(get_bsum()) << " bavg=" << llrint(get_bavg())
     << " bsdv=" << llrint(get_bsdv()) << " bmin=" << llrint(get_bmin())
     << " bmax=" << llrint(get_bmax()) << " [B]";
  os.precision(oldprec);
  return os;
}
}
