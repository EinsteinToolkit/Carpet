#include <cassert>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <Timer.hh>
#include <TimerTree.hh>



namespace Timers {
  
  using namespace std;
  
  /*********************************************************************
   Timer
  *********************************************************************/
  
  
  
  TimerTree main_timer_tree;
  TimerTree mode_timer_tree;
  
  
  
  /// Create a timer with a given name, but do not start it, and do
  /// not associate it with a point in the timer hierarchy.
  Timer::Timer(string name_p, int tree):
    d_name(name_p), d_tree(tree==0 ? &main_timer_tree : &mode_timer_tree)
  {
  }
  
  /// Destroy the timer
  Timer::~Timer()
  {
  }
  
  /// Insert the timer into the tree of timers as a child of the most
  /// recently started timer that has not been stopped. Don't start
  /// the timer. This routine ensures a timer is created even if it is
  /// never started.
  void Timer::instantiate()
  {
    TimerNode *current_timer = d_tree->current;
    if (not d_tree->root) return; // do nothing if there is no root
    assert(current_timer);
    current_timer->getChildTimer(name())->instantiate();
  }
  
  /// Start the timer and insert it into the tree of timers as a child
  /// of the most recently started timer that has not been stopped.
  void Timer::start()
  {
    TimerNode *current_timer = d_tree->current;
    if (not d_tree->root) return; // do nothing if there is no root
    assert(current_timer);
    current_timer->getChildTimer(name())->start();
  }
  
  /// Stop the timer - it must be the most recently started timer
  void Timer::stop()
  {
    TimerNode *current = d_tree->current;
    if (not d_tree->root) return; // do nothing if there is no root
    if (current->getName() != name())
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Trying to stop enclosing timer '%s' before enclosed time '%s'",
                  name().c_str(), current->getName().c_str());
    current->stop();
  }
  
  /// Return the name of the timer
  string Timer::name() const
  {
    return d_name;
  }
  
  /// Return the current time of the timer as a double
  double Timer::getTime()
  {
    return d_tree->current->getTime();
  }
  
  void Timer::outputTree(string name)
  {
    DECLARE_CCTK_PARAMETERS;
    
    TimerNode *tt = main_timer_tree.root->getChildTimer(name.c_str());
    double total_avg, total_max;
    tt->getGlobalTime(total_avg, total_max);
    tt->print
      (cout, total_max, 0, threshold_percentage, output_precision);
    mode_timer_tree.root->getGlobalTime(total_avg, total_max);
    mode_timer_tree.root->print
      (cout, total_max, 0, threshold_percentage, output_precision);
  }
  
  void Timer::outputTreeXML()
  {
    DECLARE_CCTK_PARAMETERS;
    
    main_timer_tree.root->outputXML(out_dir, CCTK_MyProc(0));
  }
  
  
  
  extern "C"
  int Timer_Startup()
  {
    // This must happen before any Timer objects are created
    main_timer_tree.root = new TimerNode(&main_timer_tree, "main");
    main_timer_tree.current = 0; // No timer has been started yet
    main_timer_tree.root->start();
    
    mode_timer_tree.root = new TimerNode(&mode_timer_tree, "meta mode");
    mode_timer_tree.current = 0; // No timer has been started yet
    mode_timer_tree.root->start();
    
    return 0;
  }
  
  extern "C"
  int Timer_Shutdown()
  {
    // main_timer_tree.root->stop();
    // mode_timer_tree.root->stop();
    
    // Delete timer tree
    delete main_timer_tree.root; main_timer_tree.root = 0;
    delete mode_timer_tree.root; mode_timer_tree.root = 0;
    
    return 0;
  }
  
} // namespace Timers
