#include <cassert>
#include <cctk.h>
#include <unistd.h>
using namespace std;

extern "C" int TestTimers2_TestCactusTimers() {
  int handle = CCTK_TimerCreate("TestTimer");
  assert(handle >= 0);
  CCTK_TimerStartI(handle);
  unsigned int remaining = sleep(1);
  CCTK_TimerStopI(handle);

  cTimerData *timer = CCTK_TimerCreateData();
  assert(timer);
  CCTK_TimerI(handle, timer);

  const cTimerVal *tv = CCTK_GetClockValue("gettimeofday", timer);
  assert(tv);

  double val = CCTK_TimerClockSeconds(tv);
  if (remaining == 0 && val < 0.5) {
    CCTK_WARN(CCTK_WARN_ALERT, "Cactus timers don't work");
  }

  return 0;
}
