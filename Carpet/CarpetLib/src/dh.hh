// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dh.hh,v 1.17 2004/03/03 16:20:19 schnetter Exp $

#ifndef DH_HH
#define DH_HH

#include <assert.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "gh.hh"
#include "vect.hh"

using namespace std;



// Forward declaration
template<int D> class ggf;
template<int D> class dh;

// Output
template<int D>
ostream& operator<< (ostream& os, const dh<D>& d);



// A data hierarchy (grid hierarchy plus ghost zones)
template<int D>
class dh {
  
  // Types
  typedef vect<int,D>    ivect;
  typedef bbox<int,D>    ibbox;
  typedef bboxset<int,D> ibset;
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect; // vector of lists
  
public:
  
  // in here, the term "boundary" means both ghost zones and
  // refinement boundaries, but does not refer to outer (physical)
  // boundaries.
  
  // ghost zones, refinement boundaries, and outer boundaries are not
  // used as sources for synchronisation.  this design choice might
  // not be good.
  
  struct dboxes {
    ibbox exterior;		// whole region (including boundaries)
    
    ibbox interior;		// interior (without boundaries)
    iblist send_mg_fine;
    iblist send_mg_coarse;
    iblist recv_mg_fine;
    iblist recv_mg_coarse;
    iblistvect send_ref_fine;
    iblistvect send_ref_coarse;
    iblistvect recv_ref_fine;
    iblistvect recv_ref_coarse;
    iblistvect send_sync;	// send while syncing
    iblistvect send_ref_bnd_fine; // sent to finer grids
    
    ibset boundaries;		// boundaries
    iblistvect recv_sync;	// received while syncing
    iblistvect recv_ref_bnd_coarse; // received from coarser grids
    ibset sync_not;		// not received while syncing (outer boundary of that level)
    ibset recv_not;		// not received while syncing or prolongating (globally outer boundary)
  };
  
private:
  
  struct dbases {
    ibbox exterior;		// whole region (including boundaries)
    ibbox interior;		// interior (without boundaries)
    ibset boundaries;		// boundaries
  };
  
  typedef vector<dboxes> mboxes; // ... for each multigrid level
  typedef vector<mboxes> cboxes; // ... for each component
  typedef vector<cboxes> rboxes; // ... for each refinement level
  
  typedef vector<dbases> mbases; // ... for each multigrid level
  typedef vector<mbases> rbases; // ... for each refinement level
  
public:				// should be readonly
  
  // Fields
  gh<D>& h;			// hierarchy
  ivect lghosts, ughosts;	// ghost zones
  
  int prolongation_order_space;	// order of spatial prolongation operator
  int buffer_width;             // buffer inside refined grids
  
  rboxes boxes;
  rbases bases;
  
  list<ggf<D>*> gfs;            // list of all grid functions
  
public:
  
  // Constructors
  dh (gh<D>& h, const ivect& lghosts, const ivect& ughosts,
      int prolongation_order_space, int buffer_width);
  
  // Destructors
  virtual ~dh ();
  
  // Helpers
  int prolongation_stencil_size () const;
  
  // Modifiers
  void recompose (const int initialise_from, const bool do_prolongate);
  
  // Grid function management
  void add (ggf<D>* f);
  void remove (ggf<D>* f);
  
  // Output
  virtual void output (ostream& os) const;
};



template<int D>
inline ostream& operator<< (ostream& os, const dh<D>& d) {
  d.output(os);
  return os;
}



#endif // DH_HH
