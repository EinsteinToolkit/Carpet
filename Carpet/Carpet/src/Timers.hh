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
    add (Timer * timer);
    
    // Remove a timer
    void
    remove (Timer * timer);
    
    // Print all timer names
    void
    printNames ()
      const;
    
    // Print timer data
    void
    printData ();
    
    // Write all timer data
    static
    void
    writeData (cGH const * cctkGH,
               char const * filename);
    
  private:
    
    // If filename is not empty, then redirect stdout to a file
    static
    int
    redirect (cGH const * cctkGH,
              char const * filename);
    
    // Redirect stdout back
    static
    void
    unredirect (int oldfd);
    
  }; // class TimerSet
  
  
  
  // A timer, which is a wrapper around a Cactus timer
  class Timer {
    
    int handle;
    bool running;
    
  public:
    
    // Create a new Cactus timer with the given name
    Timer (char const * timername);
    
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
    
    // Timer name
    char const *
    name ()
      const;
    
    // Print timer data
    void
    printData ();
    
  };
  
} // namespace Carpet
