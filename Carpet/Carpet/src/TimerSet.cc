#include <cassert>
#include <cstdio>
#include <cstring>
#include <list>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_String.h>

#if HAVE_UNISTD_H
#  include <fcntl.h>
#  include <unistd.h>
#endif

#include <defs.hh>

#include <Timers.hh>
#include <CactusTimer.hh>
#include <TimerNode.hh>
#include <TimerSet.hh>

namespace Carpet {

  using namespace std;

  // A global timer set
  TimerSet timerSet;

  // Add a timer
  void
  TimerSet::add (CactusTimer * const timer)
  {
    timers.push_back (timer);
  }

  // Remove a timer
  void
  TimerSet::remove (CactusTimer * const timer)
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
    for (list <CactusTimer *>::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++ itimer)
    {
      printf ("   [%4d] %s\n", n, (* itimer)->name().c_str());
      ++ n;
    }
  }

  // Print all timer data
  void
  TimerSet::printData ()
  {
    for (list <CactusTimer *>::const_iterator
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
    // close (1);
    // int const fd = dup (fdfile); // dup to 1, i.e., stdout again
    int const fd = dup2 (fdfile, 1); // dup to 1, i.e., stdout again
    assert (fd == 1);
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
    // close (1);
    // int const fd = dup (oldfd);
    int const fd = dup2 (oldfd, 1);
    if (not (fd == 1)) {
      fprintf(stderr, "oldfd=%d fd=%d\n", oldfd, fd);
    }
    assert (fd == 1);
    close (oldfd);
#endif
  }


  /// Reduce each timers in the set across all processes and update
  /// each timer with the reduction information.
  void TimerSet::reduce()
  {
    // Collect timer names that each process has

    // Construct union of all timer names, sort canonically and assign
    // integer identifiers

    // For each timer, identify which processes have that timer

    // Reduce the timer across all those processes (return to root proc only)

    serialise(cout);
  }

  ostream& TimerSet::serialise(ostream &os)
  {
    for (list <CactusTimer *>::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++ itimer)
    {
      (*itimer)->serialise(os);
      os << endl;
    }
    return os;
  }

/*


Each process has a list of (string,real) pairs.  I want to return a
list of these where the reals have been reduced using a reduction
operator.  Not all processes have the same strings present.
 */

} // namespace Carpet
