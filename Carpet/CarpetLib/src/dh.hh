#ifndef DH_HH
#define DH_HH

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "gh.hh"
#include "region.hh"
#include "vect.hh"

using namespace std;

#define CARPET_HAVE_BUFFER_WIDTH

// Forward declaration
class ggf;
class dh;


// A data hierarchy (grid hierarchy plus ghost zones)
class dh {
  
  // Types
public:
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect; // vector of lists
  
  typedef vector <pseudoregion_t> pvect;
  typedef vector <sendrecv_pseudoregion_t> srpvect;
  
  
  
  struct dboxes {
    
    // Region description:
    
    ibbox exterior;             // whole region (including boundaries)
    
    b2vect is_outer_boundary;
    ibset outer_boundaries;     // outer boundary
    ibbox communicated;         // exterior without outer boundary
    
    ibset boundaries;           // ghost zones
    ibbox owned;                // evolved in time
    
    ibset buffers;              // buffer zones
    ibset active;               // owned minus buffers
    
    ibset sync;                 // filled by synchronisation
    ibset bndref;               // filled by boundary prolongation
    
    // For Cactus: (these are like boundary or owned, but include the
    // outer boundary)
    ibset ghosts;               // ghost zones, as seen from Cactus
    ibbox interior;             // interior (without ghost zones)
    
    size_t memory () const;
    ostream & output (ostream & os) const;
  };
  
  struct fast_dboxes {
    
    // Communication schedule:
    
    srpvect fast_mg_rest_sendrecv;
    srpvect fast_mg_prol_sendrecv;
    srpvect fast_ref_prol_sendrecv;
    srpvect fast_ref_rest_sendrecv;
    srpvect fast_sync_sendrecv;
    srpvect fast_ref_bnd_prol_sendrecv;
    
    // Regridding schedule:
    
    srpvect fast_old2new_sync_sendrecv;
    srpvect fast_old2new_ref_prol_sendrecv;
    
    size_t memory () const;
    ostream & output (ostream & os) const;
  };
  
private:
  
  typedef vector<dboxes> cboxes; // ... for each component
  typedef vector<cboxes> rboxes; // ... for each refinement level
  typedef vector<rboxes> mboxes; // ... for each multigrid level
  
  typedef vector<fast_dboxes> fast_cboxes; // ... for each component
  typedef vector<fast_cboxes> fast_rboxes; // ... for each refinement level
  typedef vector<fast_rboxes> fast_mboxes; // ... for each multigrid level
  
  
  
  void
  setup_bboxes ();
  
public:                         // should be readonly
  
  // Fields
  gh & h;                       // hierarchy
  i2vect ghost_width;           // number of ghost zones
  i2vect buffer_width;          // number of buffer zones
  
  int prolongation_order_space; // order of spatial prolongation operator
  
  mboxes boxes;                 // grid hierarchy
  mboxes oldboxes;              // old grid hierarchy, used during regridding
  fast_mboxes fast_boxes;       // grid hierarchy
  fast_mboxes fast_oldboxes;
  
  list<ggf*> gfs;               // list of all grid functions
  
public:
  
  // Constructors
  dh (gh & h,
      i2vect const & ghosts, i2vect const & buffers,
      int prolongation_order_space);
  
  // Destructors
  ~dh ();
  
  // Helpers
  int prolongation_stencil_size () const;
  
  // Modifiers
  void regrid ();
  void recompose (int rl, bool do_prolongate);
  
private:
  int this_proc (int rl, int c) const;
  bool on_this_proc (int rl, int c) const;
  bool on_this_proc (int rl, int c, int cc) const;
  int this_oldproc (int rl, int c) const;
  bool on_this_oldproc (int rl, int c) const;
  
public:
  // Grid function management
  void add (ggf * f);
  void remove (ggf * f);
  
  // Output
  size_t memory () const;
  ostream & output (ostream & os) const;
};



inline size_t memoryof (dh::dboxes const & b)
{
  return b.memory ();
}

inline size_t memoryof (dh::fast_dboxes const & b)
{
  return b.memory ();
}

inline size_t memoryof (dh const & d)
{
  return d.memory ();
}

inline ostream & operator<< (ostream & os, dh::dboxes const & b)
{
  return b.output (os);
}

inline ostream & operator<< (ostream & os, dh::fast_dboxes const & b)
{
  return b.output (os);
}

inline ostream & operator<< (ostream & os, dh const & d)
{
  return d.output (os);
}



#endif // DH_HH
