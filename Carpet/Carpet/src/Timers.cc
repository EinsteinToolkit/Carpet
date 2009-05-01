#include <cassert>
#include <cstdio>
#include <cstring>
#include <list>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_String.h"

#if HAVE_UNISTD_H
#  include <fcntl.h>
#  include <unistd.h>
#endif

#include <defs.hh>

#include "Timers.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  // A global timer set
  static
  TimerSet timerSet;
  
  
  
  // Add a timer
  void
  TimerSet::add (Timer * const timer)
  {
    timers.push_back (timer);
  }
  
  
  
  // Remove a timer
  void
  TimerSet::remove (Timer * const timer)
  {
    timers.remove (timer);
  }
  
  
  
  // Print all timer names
  void
  TimerSet::printNames ()
    const
  {
    printf ("Timer names:\n");
    int n = 0;
    for (list <Timer *>::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++ itimer)
    {
      printf ("   [%4d] %s\n", n, (* itimer)->name());
      ++ n;
    }
  }
  
  
    
  // Print all timer data
  void
  TimerSet::printData ()
  {
    for (list <Timer *>::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++ itimer)
    {
      (* itimer)->printData ();
    }
    printf ("\n");
  }
  
  
    
  // Print all timer data
  void
  TimerSet::writeData (cGH const * const cctkGH,
                       char const * const filename)
  {
    int const oldfd = redirect (cctkGH, filename);
#if 0
    printf ("********************************************************************************\n");
#endif
    printf ("# Carpet timing information at iteration %d time %g:\n",
            cctkGH->cctk_iteration, (double) cctkGH->cctk_time);
    timerSet.printData ();
    unredirect (oldfd);
  }
  
  
  
  // If filename is not empty, then redirect stdout to a file
  int
  TimerSet::redirect (cGH const * const cctkGH,
                      char const * const filename)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (filename, "")) {
      return -1;
    }
    
#ifndef HAVE_UNISTD_H
    CCTK_WARN (1, "Cannot redirect timer output to a file; the operating system does not support this");
    return -1;
#else
    
    int const myproc = CCTK_MyProc (cctkGH);
    char fullname [10000];
    Util_snprintf (fullname, sizeof fullname,
                   "%s/%s.%04d.txt", out_dir, filename, myproc);
    
    int flags = O_WRONLY | O_CREAT | O_APPEND; // append
    static bool first_time = true;
    if (first_time) {
      first_time = false;
      if (IO_TruncateOutputFiles (cctkGH)) {
        flags = O_WRONLY | O_CREAT | O_TRUNC; // truncate
      }
    }
    
    // Temporarily redirect stdout
    fflush (stdout);
    int const oldfd = dup (1);  // fd 1 is stdout
    int const mode = 0644;      // rw-r--r--, or a+r u+w
    int const fdfile = open (fullname, flags, mode);
    if (fdfile < 0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open timer output file \"%s\"", fullname);
      close (oldfd);
      return -1;
    }
    close (1);
    dup (fdfile);               // dup to 1, i.e., stdout again
    close (fdfile);
    return oldfd;
#endif
  }
  
  
  
  // Redirect stdout back
  void
  TimerSet::unredirect (int const oldfd)
  {
    if (oldfd < 0) return;
    
#ifdef HAVE_UNISTD_H
    fflush (stdout);
    close (1);
    dup (oldfd);
    close (oldfd);
#endif
  }
  
  
  
  // Create a new Cactus timer with the given name
  Timer::Timer (char const * const timername)
    : running (false)
  {
    assert (timername);
    handle = CCTK_TimerCreate (timername);
    assert (handle >= 0);
    
    timerSet.add (this);
  }
  
  
  
  // Destroy a timer
  Timer::~Timer ()
  {
    timerSet.remove (this);
    check (not CCTK_TimerDestroyI (handle));
  }
  
  
  
  // Timer name
  char const *
  Timer::name ()
    const
  {
    char const * const timername = CCTK_TimerName (handle);
    assert (timername);
    return timername;
  }
  
  
  
  // Print timer data
  void
  Timer::printData ()
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
    
    printf ("%s:", name());
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
  
} // namespace Carpet
