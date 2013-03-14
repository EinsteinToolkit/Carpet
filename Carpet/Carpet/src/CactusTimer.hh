#ifndef CACTUSTIMER_HH
#define CACTUSTIMER_HH

#include <cctk.h>

#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

namespace Carpet {

  using namespace std;

/** The CactusTimer class wraps the Cactus timer mechanism.  All times
    are returned as doubles for now. */

  class CactusTimer
  {
    int handle;
    bool running;

  public:

    /// Create a new Cactus timer with the given name
    CactusTimer (const string &timername);

    /// Destroy a timer
    ~CactusTimer ();

    /// Start the timer
    void start ();

    /// Stop the timer
    void stop ();

    /// Reset the timer
    void reset ();

    /// Timer name
    string name () const;

    /// Return the current time of the timer as a double
    double getTime();

    /// Return all clock names
    vector<pair<string,string> > getAllTimerNames() const;

    /// Return all clock values of the timer as double
    vector<double> getAllTimerValues();

    /// Print timer data
    void printData ();

    ostream& serialise(ostream &os);

  private:

    // Output (debug) messages that a timer is starting or stopping
    void
    msgStart ()
      const;

    void
    msgStop ()
      const;

  };
}  // namespace Carpet

#endif // CACTUSTIMER_HH
