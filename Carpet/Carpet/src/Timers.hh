#include <list>

#include "cctk.h"



namespace Carpet {
  
  using namespace std;
  
  
  
  class Timer;
  
  
  
  // A set of timers
  class TimerSet {
    
    list <Timer *> timers;
    
  public:
    
    // Add a timer
    void
    add (Timer * const timer);
    
    // Remove a timer
    void
    remove (Timer * const timer);
    
    // Print all timer data
    void
    printData (cGH const * const cctkGH,
               char const * const filename);
    
  private:
    
    int fdsave;
    
    // If filename is not empty, then redirect stdout to a file
    void
    redirect (cGH const * const cctkGH,
              char const * const filename);
    
    // Redirect stdout back
    void
    unredirect ();
    
  }; // class TimerSet
  
  // A global timer set
  extern TimerSet timerSet;
  
  
  
  // A timer, which is a wrapper around a Cactus timer
  class Timer {
    
    int handle;
    bool running;
    TimerSet & timerSet;
    
  public:
    
    // Create a new Cactus timer with the give name, which belongs to
    // a certain timer set
    Timer (TimerSet & timerSet_,
           char const * const timername);
    
    // Destroy a timer
    ~Timer ();
    
    // Start the timer
    void
    start ()
    {
      running = true;
      CCTK_TimerStartI (handle);
    }
    
    // Stop the timer
    void
    stop ()
    {
      CCTK_TimerStopI (handle);
      running = false;
    }
    
    // Reset the timer
    void
    reset ()
    {
      CCTK_TimerResetI (handle);
    }
    
    // Print timer data
    void
    printData ();
    
  };
  
} // namespace Carpet
