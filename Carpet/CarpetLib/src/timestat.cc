#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <unistd.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif

#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"



namespace CarpetLib {
  
  using namespace std;
  
  
  
  // The timer type
  enum timer_type {
    timer_unset, timer_MPI_Wtime, timer_rdtsc, timer_rtc, timer_none
  };
  static timer_type thetimer = timer_unset;
  
  
  
  // A faster timing routine for i386 processors:
  // Read the Intel CPU time stamp counter
  static
  double rdtsc_cputick = 1.0;
  
  static inline
  double
  rdtsc ()
  {
    // Test for i386 (should use autoconf instead)
#if defined(__i386__) || defined(__x86_64__)
    // Serialise using cpuid
    // (This is strictly necessary only on some systems)
    {
      unsigned long eax;
      asm volatile ("movl %%ebx, %%esi; cpuid; movl %%esi, %%ebx" :
                    "=a" (eax) : "a" (0) : "edx", "ecx", "esi");
    }
    // Using "=A" does not work on x86_64, where this uses rax instead
    //    of (edx,eax)
    // unsigned long long val;
    // __asm__ __volatile__("rdtsc" : "=A" (val) : );
    unsigned long eax, edx;
    asm volatile ("rdtsc" : "=a" (eax), "=d" (edx));
    return rdtsc_cputick * ((unsigned long long) edx << 32 | eax);
#else
    static bool have_warned = false;
    if (not have_warned) {
      CCTK_WARN (1, "The Intel rdtsc timer is not supported");
      have_warned = true;
    }
    return 0;
#endif
  }
  
  static
  void
  init_rdtsc ()
  {
    // Make three warmup measurements
    double const rdummy1 = rdtsc ();
    double const rdummy2 = rdtsc ();
    double const rdummy3 = rdtsc ();
    double const rstart = rdtsc ();
    double const wstart = MPI_Wtime ();
    int const ierr = usleep (1000 * 1000);
    double const rend = rdtsc ();
    double const wend = MPI_Wtime ();
    if (ierr) {
      CCTK_WARN (1, "Could not determine a reliable rdtsc timer resolution");
    }
    rdtsc_cputick *= (wend - wstart) / (rend - rstart);
  }
  
  
  
  // A faster timing routine for AIX
  static inline
  double
  rtc ()
  {
    // Test for AIX (should use autoconf instead)
#ifdef __IBMC__
    timebasestruct_t val;
    read_real_time (& val, TIMEBASE_SZ);
    time_base_to_time (& val, TIMEBASE_SZ);
    return val.tb_high + 1.0e-9 * val.tb_low;
#else
    static bool have_warned = false;
    if (not have_warned) {
      CCTK_WARN (1, "The AIX rtc timer not supported");
      have_warned = true;
    }
    return 0;
#endif
  }
  
  
  
  // Call a timer
  static
  double
  call_timer ()
  {
    DECLARE_CCTK_PARAMETERS;
    if (thetimer == timer_unset) {
      if (CCTK_EQUALS (timestat_timer, "MPI_Wtime")) {
        thetimer = timer_MPI_Wtime;
      } else if (CCTK_EQUALS (timestat_timer, "rdtsc")) {
        thetimer = timer_rdtsc;
        init_rdtsc ();
      } else if (CCTK_EQUALS (timestat_timer, "rtc")) {
        thetimer = timer_rtc;
      } else if (CCTK_EQUALS (timestat_timer, "none")) {
        thetimer = timer_none;
      } else {
        assert (0);
      }
    }
    switch (thetimer) {
    case timer_MPI_Wtime:
      return MPI_Wtime ();
    case timer_rdtsc:
      return rdtsc ();
    case timer_rtc:
      return rtc ();
    case timer_none:
      return 0.0;
    default:
      assert (0);
    }
    abort ();
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
    if (thetimer == timer_none) return;
    assert (not running);
    running = true;
    starttime = call_timer ();
  }
  
  
  
  // Stop the timer
  void
  Timer::stop (double const b)
  {
    DECLARE_CCTK_PARAMETERS;
    if (thetimer == timer_none) return;
    assert (running);
    running = false;
    double const endtime = call_timer ();
    addstat (endtime - starttime, b);
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
