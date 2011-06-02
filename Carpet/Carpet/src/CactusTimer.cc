#include <cassert>
#include <cstdio>
#include <cstring>
#include <list>
#include <iomanip>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_String.h>

#if HAVE_UNISTD_H
#  include <fcntl.h>
#  include <unistd.h>
#endif

#include <defs.hh>

#include "CactusTimer.hh"
#include "TimerSet.hh"

namespace Carpet
{
  using namespace std;

  // Create a new Cactus timer with the given name
  CactusTimer::CactusTimer (const string &timername)
    : running (false)
  {
//    cout << "CactusTimer::CactusTimer(): name = " << timername << endl;
    handle = CCTK_TimerCreate (timername.c_str());
    assert (handle >= 0);

    timerSet.add (this);
  }

  // Destroy a timer
  CactusTimer::~CactusTimer ()
  {
    timerSet.remove (this);
    check (not CCTK_TimerDestroyI (handle));
  }

  // Start the timer
  void CactusTimer::start ()
  {
    msgStart ();
//    cout << "CactusTimer::start this = " << this << endl;
    running = true;
    CCTK_TimerStartI (handle);
  }

  // Stop the timer
  void CactusTimer::stop ()
  {
    CCTK_TimerStopI (handle);
    running = false;
    msgStop ();
  }

  // Reset the timer
  void CactusTimer::reset ()
  {
    CCTK_TimerResetI (handle);
  }

  // Timer name
  string CactusTimer::name () const
  {
    char const * const timername = CCTK_TimerName (handle);
    assert (timername);
    return string(timername);
  }

  double CactusTimer::getTime()
  {
    DECLARE_CCTK_PARAMETERS;

    bool const was_running = running;
    if (was_running) stop();

    static cTimerData * timer = 0;
    if (not timer) timer = CCTK_TimerCreateData ();
    assert (timer);
    CCTK_TimerI (handle, timer);

    double val = 0; // All these timers will be returned as doubles

    const cTimerVal  *tv = CCTK_GetClockValue(timer_xml_clock, timer);

    if (tv != NULL)
    {
      val = CCTK_TimerClockSeconds(tv);
    }
    else
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Clock \"%s\" not found for timer #%d \"%s\"",
                 timer_xml_clock, handle, CCTK_TimerName(handle));
      val = -1;
    }

    // for (int i=0; i<timer->n_vals; ++i) {
    //   switch (timer->vals[i].type) {
    //   case val_int:
    //     val = timer->vals[i].val.i;
    //     break;
    //   case val_long:
    //     val = timer->vals[i].val.l;
    //     break;
    //   case val_double:
    //     val = timer->vals[i].val.d;
    //     break;
    //   case val_none:
    //     val = 0; // ??
    //     break;
    //   default:
    //     assert (0);
    //   }
    // }

    if (was_running) start();

    return val;
  }


  // Print timer data
  void CactusTimer::printData ()
  {
    bool const was_running = running;
    if (was_running) stop();

#if 0
    check (not CCTK_TimerPrintDataI (handle, -1)); // -1 means: all clocks
#endif

    static cTimerData * timer = 0;
    if (not timer) timer = CCTK_TimerCreateData ();
    assert (timer);
    CCTK_TimerI (handle, timer);

    static bool firsttime = true;
    if (firsttime) {
      printf ("# 1: timer name");
      for (int i=0; i<timer->n_vals; ++i) {
        printf (" %d: %s [%s]",
                i+2, timer->vals[i].heading, timer->vals[i].units);
      }
      printf ("\n");
      firsttime = false;
    }

    printf ("%s:", name().c_str());
    for (int i=0; i<timer->n_vals; ++i) {
      switch (timer->vals[i].type) {
      case val_int:
        printf (" %d", timer->vals[i].val.i);
        break;
      case val_long:
        printf (" %ld", timer->vals[i].val.l);
        break;
      case val_double:
        printf (" %g", timer->vals[i].val.d);
        break;
      case val_none:
        break;
      default:
        assert (0);
      }
    }
    printf ("\n");

    if (was_running) start();
  }

  // Output (debug) messages that a timer is starting or stopping
  void CactusTimer::msgStart ()  const
  {
    DECLARE_CCTK_PARAMETERS;
    if (timers_verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "Timer \"%s\" starting", name().c_str());
    }
  }

  void CactusTimer::msgStop () const
  {
    DECLARE_CCTK_PARAMETERS;
    if (timers_verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "Timer \"%s\" stopping", name().c_str());
    }
  }

  ostream& CactusTimer::serialise(ostream &os)
  {
    os << scientific << setprecision(19) << getTime() << " " << name();
    return os;
  }

}  // namespace Carpet
