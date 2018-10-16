#include <Timer.hh>
#include <cassert>
#include <cctk.h>
#include <unistd.h>
using namespace std;

extern "C" int TestTimers2_TestTimers() {
  Timers::Timer timer("Test");
  timer.start();
  unsigned int remaining = sleep(1);
  double val = timer.getTime();
  timer.stop();
  if (remaining == 0 && val < 0.5) {
    CCTK_WARN(CCTK_WARN_ALERT, "Timers don't work");
  }

  return 0;
}
