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

#include "Timers.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  // A global timer set
  TimerSet & timerSet ()
  {
    static TimerSet timerSet_;
    return timerSet_;
  }
  
  
  
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
  TimerSet::print ()
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
  TimerSet::printData (cGH const * const cctkGH,
                       char const * const filename)
  {
    redirect (cctkGH, filename);
    printf ("********************************************************************************\n");
    printf ("Carpet timing information at iteration %d time %g:\n",
            cctkGH->cctk_iteration, (double) cctkGH->cctk_time);
    for (list <Timer *>::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++ itimer)
    {
      (* itimer)->printData ();
    }
    printf ("********************************************************************************\n");
    unredirect ();
  }
  
  
  
  // If filename is not empty, then redirect stdout to a file
  void
  TimerSet::redirect (cGH const * const cctkGH,
                      char const * const filename)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (filename, "")) {
      fdsave = -1;
      return;
    }
    
#ifndef HAVE_UNISTD_H
    CCTK_WARN (1, "Cannot redirect timer output to a file; the operating system does not support this");
    return;
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
    fdsave = dup (1);          // fd 1 is stdout
    int const mode = 0644;     // rw-r--r--, or a+r u+w
    int const fdfile = open (fullname, flags, mode);
    if (fdfile < 0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open timer output file \"%s\"", fullname);
      close (fdsave);
      fdsave = -1;
      return;
    }
    close (1);
    dup (fdfile);             // dup to 1, i.e., stdout again
    close (fdfile);
#endif
  }
  
  
  
  // Redirect stdout back
  void
  TimerSet::unredirect ()
  {
    if (fdsave < 0) return;
    
#ifdef HAVE_UNISTD_H
    fflush (stdout);
    close (1);
    dup (fdsave);
    close (fdsave);
#endif
  }
  
  
  
  // Create a new Cactus timer with the give name, which belongs to a
  // certain timer set
  Timer::Timer (TimerSet & timerSet_, char const * const name)
    : running (false),
      timerSet (timerSet_)
  {
    assert (name);
    handle = CCTK_TimerCreate (name);
    assert (handle >= 0);
    
    timerSet.add (this);
  }
  
  
  
  // Destroy a timer
  Timer::~Timer ()
  {
    timerSet.remove (this);
    int const ierr = CCTK_TimerDestroyI (handle);
    assert (not ierr);
  }
  
  
  
  // Timer name
  char const *
  Timer::name ()
    const
  {
    char const * const name_ = CCTK_TimerName (handle);
    assert (name_);
    return name_;
  }
  
  
  
  // Print timer data
  void
  Timer::printData ()
  {
    bool const was_running = running;
    if (was_running) stop();
    int const ierr = CCTK_TimerPrintDataI (handle, -1); // -1 means: all clocks
    assert (not ierr);
    if (was_running) start();
  }
  
} // namespace Carpet
