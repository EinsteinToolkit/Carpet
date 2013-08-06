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
#include <TimerNode.hh>
#include "variables.hh"


namespace Carpet {
  
  using namespace std;

/*********************************************************************
 Timer
 *********************************************************************/

  /// Create a timer with a given name, but do not start it, and do
  /// not associate it with a point in the timer hierarchy.
  Timer::Timer (const string &name_p) : d_name(name_p)
  {
    d_tree = &main_timer_tree;
  }

  Timer::Timer (const string &name_p, TimerTree *tree) : d_name(name_p), d_tree(tree)
  {
  }

  /// Destroy a timer
  Timer::~Timer ()
  {
  }

  /// Insert the timer into the tree of timers as a child of the most
  /// recently started timer that has not been stopped. Don't start
  /// the timer. This routine ensures a timer is created even if it is
  /// never started.
  void Timer::instantiate ()
  {
    TimerNode *current_timer = d_tree->current;
    assert(current_timer);
    current_timer->getChildTimer(name())->instantiate();
  }

  /// Start the timer and insert it into the tree of timers as a child
  /// of the most recently started timer that has not been stopped.
  void Timer::start ()
  {
    TimerNode *current_timer = d_tree->current;
    if (not d_tree->root) return; // do nothing if there is no root
    assert(current_timer);
    current_timer->getChildTimer(name())->start();
  }

  /// Stop the timer - it must be the most recently started timer
  void Timer::stop ()
  {
    TimerNode *current = d_tree->current;
    if (not d_tree->root) return; // do nothing if there is no root
    if (current->getName() != name())
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Trying to stop enclosing timer '%s' before enclosed time '%s'",
                  name().c_str(), current->getName().c_str());
    current->stop();
  }

  /// Return the name of the timer
  string Timer::name () const
  {
    return d_name;
  }

  /// Return the current time of the timer as a double
  double Timer::getTime ()
  {
    return d_tree->current->getTime();
  }
} // namespace Carpet
