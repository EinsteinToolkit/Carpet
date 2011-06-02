#ifndef TIMERSET_HH
#define TIMERSET_HH

#include <iostream>
#include <list>

#include <cctk.h>
#include "CactusTimer.hh"

//class Carpet::TimerSet;

//ostream& operator <<(ostream &os, const Carpet::TimerSet &obj);

namespace Carpet {

  class TimerSet;
  extern TimerSet timerSet;

  using namespace std;

  // A set of timers
  class TimerSet {

    list <CactusTimer *> timers;

  public:

    // Add a timer
    void
    add (CactusTimer * timer);

    // Remove a timer
    void
    remove (CactusTimer * timer);

    // Print all timer names
    void
    printNames ()
      const;

    // Print timer data
    void
    printData ();

    // Write all timer data
    static void writeData (cGH const * cctkGH, char const * filename);

    void reduce();

    ostream& serialise(ostream &os);

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

} // namespace Carpet

#endif // TIMERSET_HH
