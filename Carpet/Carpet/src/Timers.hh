#ifndef TIMERS_HH
#define TIMERS_HH

#include <iostream>
#include <list>

#include <cctk.h>



namespace Carpet {

  using namespace std;

/**

This class allows the user to instrument their code with named timers
which can then be later queried to determine the amount of time spent
in the code between "start" and "end" calls.  The sequence of start
and end calls of different timers determines a dynamical hierarchical
tree structure implemented by the TimerNode class.

To use this class, create a timer object with a particular name:

  Timer timer("MyTimer")

Now wrap the code to be timed with start() and stop() calls:

  timer.start()

  some code

  timer.stop()

You can start and stop a timer multiple times.  The timer will be
created as a child of whatever timer is current (i.e. has been started
and not stopped) at the time of the first start() call.  Any timers
which are started between the start() and stop(), whether or not they
are in the same file, will be stored as children of this timer.

Timer objects must be started and stopped in a non-overlapping manner.
Specifically, a timer cannot be stopped if it is not the most recently
started timer.  Doing so will generate an error.

Timer objects can be allocated as "static" or not - it does not matter.

*/
  class TimerTree;

  class Timer {

  public:

    Timer (const string &name);
    Timer (const string &name, TimerTree *tree);
    ~Timer ();

    void instantiate ();
    void start ();
    void stop ();
    string name () const;
    double getTime();

  private:

    string d_name;
    TimerTree *d_tree;
  };

  // Macros for using timers in a convenient manner
  
#define TIMING_BEGIN(name)                      \
  do {                                          \
    static Carpet::Timer timer (name);          \
    timer.start();                              \
    {
#define TIMING_END                              \
    }                                           \
    timer.stop();                               \
  } while (0)
  
} // namespace Carpet

#endif // TIMERS_HH
