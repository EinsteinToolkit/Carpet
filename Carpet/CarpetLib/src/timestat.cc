#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <unistd.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"



namespace CarpetLib {
  
  using namespace std;
  
  
  
  static
  bool have_cputick = false;
  
  static
  double cputick = 0.0;
  
  static
  void
  calculate_cputick ()
  {
    // Make a few warm-up measurements
    getticks ();
    getticks ();
    getticks ();
    
#if 0
    // Use usleep to calibrate the timer
    for (int i=0; i<10; ++i) {
      useconds_t const waittime = 1000 * 1000;
      ticks const rstart = getticks ();
      int const ierr = usleep (waittime);
      ticks const rend = getticks ();
      cputick = waittime / 1.0e6 / elapsed (rend, rstart);
      if (not ierr) goto done;
    }
    CCTK_WARN (1, "Could not determine timer resolution");
  done:;
#endif
    
#if 1
    // Use MPI_Wtime to calibrate the timer
    ticks const rstart = getticks ();
    double const wstart = MPI_Wtime ();
    while (MPI_Wtime() < wstart + 1.0) {
      // do nothing, just wait
    }
    ticks const rend = getticks ();
    double const wend = MPI_Wtime ();
    cputick = (wend - wstart) / elapsed (rend, rstart);
#endif
    
    have_cputick = true;
  }
  
  
  
  // Call a timer
  static
  ticks
  call_timer ()
  {
    return getticks();
  }
  
  
  
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
  
  
  
  // Output all timer names
  void
  TimerSet::outputNames (ostream & os)
    const
  {
    os << "Timer names:" << eol;
    int n = 0;
    for (list <Timer *>::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++ itimer)
    {
      os << "   [" << setw (4) << setfill ('0') << n << "] "
         << (* itimer)->name() << eol;
      ++ n;
    }
  }
  
  
  
  // Output all timer data
  void
  TimerSet::outputData (ostream & os)
    const
  {
    for (list <Timer *>::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++ itimer)
    {
      os << * (* itimer);
    }
  }
  
  
  
  // Create a new timer with the given name
  Timer::Timer (char const * const timername_)
    : timername (timername_)
  {
    assert (timername_);
    if (not have_cputick) calculate_cputick ();
    resetstats ();
    timerSet.add (this);
  }
  
  
  
  // Destroy a timer
  Timer::~Timer ()
  {
    timerSet.remove (this);
  }
  
  
  
  // Reset the statistics
  void
  Timer::resetstats ()
  {
    wtime  = 0.0;
    wtime2 = 0.0;
    wmin   = 0.0;
    wmax   = 0.0;
    
    bytes  = 0.0;
    bytes2 = 0.0;
    bmin   = 0.0;
    bmax   = 0.0;
    
    count  = 0.0;
    
    running = false;
  }
  
  
  
  // Add statistics of a timing operation
  void
  Timer::addstat (double const t,
                  double const b)
  {
    wtime  += t;
    wtime2 += pow (t, 2);
    wmin   = min (wmin, t);
    wmax   = max (wmax, t);
    
    bytes  += b;
    bytes2 += pow (b, 2);
    bmin   = min (bmin, b);
    bmax   = max (bmax, b);
    
    ++ count;
  }
  
  
  
  // Start the timer
  void
  Timer::start ()
  {
    DECLARE_CCTK_PARAMETERS;
    assert (not running);
    running = true;
    starttime = call_timer ();
  }
  
  
  
  // Stop the timer
  void
  Timer::stop (double const b)
  {
    DECLARE_CCTK_PARAMETERS;
    assert (running);
    running = false;
    ticks const endtime = call_timer ();
    addstat (elapsed (endtime, starttime), b);
  }
  
  
  
  // Reset the timer
  void
  Timer::reset ()
  {
    resetstats ();
  }
  
  
  
  // Timer name
  string
  Timer::name ()
    const
  {
    return timername;
  }
  
  
  
  // Output timer data
  void
  Timer::outputData (ostream & os)
    const
  {
    double avg, stddev, bavg, bstddev;
    if (count == 0.0) {
      avg     = 0.0;
      stddev  = 0.0;
      bavg    = 0.0;
      bstddev = 0.0;
    } else {
      avg     = wtime / count;
      stddev  = sqrt (max (0.0, wtime2 / count - pow (avg, 2)));
      bavg    = bytes / count;
      bstddev = sqrt (max (0.0, bytes2 / count - pow (bavg, 2)));
    }
    
    os << timername << ":"
       << " cnt: " << count
       << "   time: sum: " << wtime
       << " avg: " << avg
       << " stddev: " << stddev
       << " min: " << wmin
       << " max: " << wmax
       << "   bytes: sum: " << bytes
       << " avg: " << bavg
       << " stddev: " << bstddev
       << " min: " << bmin
       << " max: " << bmax
       << eol;
  }
  
  
  
  extern "C" {
    void
    CarpetLib_printtimestats (CCTK_ARGUMENTS);
  }
  
  void
  CarpetLib_printtimestats (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if ((print_timestats_every > 0 and
         cctk_iteration % print_timestats_every == 0))
    {
      ostringstream filenamebuf;
      filenamebuf << out_dir << "/" << timestat_file
                  << "." << setw(4) << setfill('0') << dist::rank()
                  << ".txt";
      string const filename = filenamebuf.str();
      
      ofstream file;
      static bool do_truncate = true;
      if (do_truncate) {
        if (not IO_TruncateOutputFiles (cctkGH)) {
          do_truncate = false;
        }
      }
      if (do_truncate) {
        do_truncate = false;
        file.open (filename.c_str(), ios::out | ios::trunc);
      } else {
        file.open (filename.c_str(), ios::out | ios::app);
      }
      
      static bool do_print_info = true;
      if (do_print_info) {
        do_print_info = false;
        if (CCTK_IsFunctionAliased ("UniqueBuildID")) {
          char const * const build_id =
            static_cast <char const *> (UniqueBuildID (cctkGH));
          file << "Build ID: " << build_id << eol;
        }
        if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
          char const * const job_id =
            static_cast <char const *> (UniqueSimulationID (cctkGH));
          file << "Simulation ID: " << job_id << eol;
        }
        file << "Running on " << dist::size() << " processors" << eol;
      } // if do_print_info
      
      file << "********************************************************************************" << eol
           << "CarpetLib timing information at iteration " << cctkGH->cctk_iteration << " time " << cctkGH->cctk_time << ":" << eol
           << timerSet;
      
      file.close ();
      
    } // if print_timestats
    
  }
  
} // namespace CarpetLib
