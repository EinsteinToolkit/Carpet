#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include <sys/time.h>
#include <unistd.h>

#ifdef CCTK_MPI
#  include <mpi.h>
#else
#  include "nompi.h"
#endif

#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"



namespace CarpetLib {
  
  using namespace std;
  
  
  
  static
  bool have_cputick = false;
  
  // CPU tick time in seconts
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
  TimerSet* timerSet = NULL;
  
  
  
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
    if (not timerSet) timerSet = new TimerSet;
    timerSet->add (this);
  }
  
  
  
  // Destroy a timer
  Timer::~Timer ()
  {
    assert (timerSet);
    timerSet->remove (this);
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
    if (use_ipm_timing_regions) {
      MPI_Pcontrol (+1, timername.c_str());
    }
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
    if (use_ipm_timing_regions) {
      MPI_Pcontrol (-1, timername.c_str());
    }
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
       << "   time: sum: " << cputick * wtime
       << " avg: " << cputick * avg
       << " stddev: " << cputick * stddev
       << " min: " << cputick * wmin
       << " max: " << cputick * wmax
       << "   bytes: sum: " << bytes
       << " avg: " << bavg
       << " stddev: " << bstddev
       << " min: " << bmin
       << " max: " << bmax
       << eol;
  }
  
  
  
  // Fortran wrappers
  extern "C" {
    
    // In Fortran, a timer should be declared and used like this:
    //
    //    ! Save the timer handle, and initialise it to zero; this
    //    ! ensures that the timer is created only once:
    //    CCTK_POINTER, save :: timer = 0
    //    call Timer_create (timer, "Name")
    //
    //    ! Start the timer:
    //    call Timer_start (timer)
    //
    //    ! Stop the timer, and pass the number of bytes:
    //    CCTK_REAL :: bytes
    //    bytes = ...
    //    call Timer_stop (timer, bytes)
    
    void CCTK_FCALL
    CCTK_FNAME(Timer_create) (CCTK_POINTER * timer, ONE_FORTSTRING_ARG)
    {
      if (*timer != 0) return;   // create the timer only once
      ONE_FORTSTRING_CREATE (timername);
      *timer = new Timer(timername);
      free (timername);
    }
    
    void CCTK_FCALL
    CCTK_FNAME(Timer_destroy) (CCTK_POINTER * timer)
    {
      if (*timer == 0) return;   // delete the timer only if it has been created
      delete (Timer*)*timer;
      *timer = 0;
    }
    
    void CCTK_FCALL
    CCTK_FNAME(Timer_start) (CCTK_POINTER * timer)
    {
      assert (*timer != 0);
      ((Timer*)*timer)->start();
    }
    
    void CCTK_FCALL
    CCTK_FNAME(Timer_stop) (CCTK_POINTER * timer, CCTK_REAL const * b)
    {
      assert (*timer != 0);
      ((Timer*)*timer)->stop(*b);
    }
    
    void CCTK_FCALL
    CCTK_FNAME(Timer_reset) (CCTK_POINTER * timer)
    {
      assert (*timer != 0);
      ((Timer*)*timer)->reset();
    }
    
  } // extern "C"
  
  
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
          char const * const sim_id =
            static_cast <char const *> (UniqueSimulationID (cctkGH));
          file << "Simulation ID: " << sim_id << eol;
        }
        file << "Running with " << dist::size() << " processes and " << dist::total_num_threads() << " threads" << eol;
      } // if do_print_info
      
      if (not timerSet) timerSet = new TimerSet;
      file << "********************************************************************************" << eol
           << "CarpetLib timing information at iteration " << cctkGH->cctk_iteration << " time " << cctkGH->cctk_time << ":" << eol
           << *timerSet;
      
      file.close ();
      
    } // if print_timestats
    
  }
  
  
  
  struct t_cycleclock {
    double total;
    double total_squared;
    double min_total;
    double max_total;
    double count;
    ticks last;
    
    t_cycleclock ()
    {
      reset();
    }
    
    ~t_cycleclock ()
    {
    }
    
    void start ()
    {
      last = getticks();
    }
    
    void stop ()
    {
      ticks const current = getticks();
      double const difference = elapsed (current, last);
      total += difference;
      total_squared += pow (difference, 2);
      min_total = min_total == 0.0 ? difference : min (min_total, difference);
      max_total = max (min_total, difference);
      count += 1;
    }
    
    void reset ()
    {
      total         = 0.0;
      total_squared = 0.0;
      min_total     = 0.0;      // numeric_limits<double>::max();
      max_total     = 0.0;
      count         = 0.0;
      // last          = 0.0;
    }
    
  };
  
  
  
  void * cycleclock_create (int const timernum)
  {
    return new t_cycleclock;
  }
  
  void cycleclock_destroy (int const timernum, void * const data)
  {
    if (data) {
      delete static_cast<t_cycleclock*> (data);
    }
  }
  
  void cycleclock_start (int const timernum, void * const data)
  {
    static_cast<t_cycleclock*> (data) -> start();
  }
  
  void cycleclock_stop (int const timernum, void * const data)
  {
    static_cast<t_cycleclock*> (data) -> stop();
  }
  
  void cycleclock_reset (int const timernum, void * const data)
  {
    static_cast<t_cycleclock*> (data) -> reset();
  }
  
  void cycleclock_get (int const timernum, void * const data_,
                       cTimerVal * const vals)
  {
    t_cycleclock const & data = * static_cast<t_cycleclock const *> (data_);
    
    // Total time
    vals[0].type       = val_double;
    vals[0].heading    = "cycle";
    vals[0].units      = "secs";
    vals[0].val.d      = data.total;
    vals[0].seconds    = cputick * vals[0].val.d;
    vals[0].resolution = cputick;
    
    // Average
    vals[1].type       = val_double;
    vals[1].heading    = "cycle[avg]";
    vals[1].units      = "secs";
    vals[1].val.d      = data.count == 0.0 ? 0.0 : data.total / data.count;
    vals[1].seconds    = cputick * vals[1].val.d;
    vals[1].resolution = cputick;
    
    // Standard deviation
    vals[2].type       = val_double;
    vals[2].heading    = "cycle[stddev]";
    vals[2].units      = "secs";
    vals[2].val.d      = (data.count == 0.0
                          ? 0.0
                          : sqrt (abs (data.total_squared * data.count -
                                       pow (data.total, 2)) / data.count));
    vals[2].seconds    = cputick * vals[2].val.d;
    vals[2].resolution = cputick;
    
    // Minimum
    vals[3].type       = val_double;
    vals[3].heading    = "cycle[min]";
    vals[3].units      = "secs";
    vals[3].val.d      = data.min_total;
    vals[3].seconds    = cputick * vals[3].val.d;
    vals[3].resolution = cputick;
    
    // Maximum
    vals[4].type       = val_double;
    vals[4].heading    = "cycle[max]";
    vals[4].units      = "secs";
    vals[4].val.d      = data.max_total;
    vals[4].seconds    = cputick * vals[4].val.d;
    vals[4].resolution = cputick;
  }
  
  void cycleclock_set (int const timernum, void * const data_,
                       cTimerVal * const vals)
  {
    t_cycleclock & data = * static_cast<t_cycleclock * restrict> (data_);
    
    data.reset();               // punt
    data.total = vals[0].val.d;
  }
  
  extern "C" {
    int CarpetLib_registercycleclock (void);
  }
  
  int CarpetLib_registercycleclock (void)
  {
    if (not have_cputick) calculate_cputick ();
    
    cClockFuncs functions;
    functions.n_vals  = 5;
    functions.create  = cycleclock_create;
    functions.destroy = cycleclock_destroy;
    functions.start   = cycleclock_start;
    functions.stop    = cycleclock_stop;
    functions.reset   = cycleclock_reset;
    functions.get     = cycleclock_get;
    functions.set     = cycleclock_set;
    
    CCTK_ClockRegister("cycle", &functions);
    
    return 0;
  }
  
} // namespace CarpetLib
