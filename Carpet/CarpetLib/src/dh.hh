#ifndef DH_HH
#define DH_HH

#include <cassert>
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
    
    // Communication schedule:
    
    // ref_prol_recv[cc] and ref_rest_send[cc] determine what needs to
    // be sent from and received from cc for prolongation to this box
    
    iblist mg_rest_recv;
    iblist mg_rest_send;
    iblist mg_prol_recv;
    iblist mg_prol_send;
    iblistvect ref_prol_recv;
    iblistvect ref_prol_send;
    iblistvect ref_rest_recv;
    iblistvect ref_rest_send;
    iblistvect sync_recv;
    iblistvect sync_send;
    iblistvect ref_bnd_prol_recv;
    iblistvect ref_bnd_prol_send;
    
    pvect fast_mg_rest_recv;
    pvect fast_mg_rest_send;
    pvect fast_mg_prol_recv;
    pvect fast_mg_prol_send;
    pvect fast_ref_prol_recv;
    pvect fast_ref_prol_send;
    pvect fast_ref_rest_recv;
    pvect fast_ref_rest_send;
    pvect fast_sync_recv;
    pvect fast_sync_send;
    pvect fast_ref_bnd_prol_recv;
    pvect fast_ref_bnd_prol_send;
    
    srpvect fast_mg_rest_sendrecv;
    srpvect fast_mg_prol_sendrecv;
    srpvect fast_ref_prol_sendrecv;
    srpvect fast_ref_rest_sendrecv;
    srpvect fast_sync_sendrecv;
    srpvect fast_ref_bnd_prol_sendrecv;
    
    // Regridding schedule:
    
    iblistvect old2new_sync_recv;
    iblistvect old2new_sync_send;
    iblistvect old2new_ref_prol_recv;
    iblistvect old2new_ref_prol_send;
    
    pvect fast_old2new_sync_recv;
    pvect fast_old2new_sync_send;
    pvect fast_old2new_ref_prol_recv;
    pvect fast_old2new_ref_prol_send;
    
    srpvect fast_old2new_sync_sendrecv;
    srpvect fast_old2new_ref_prol_sendrecv;
    
    ostream & output (ostream & os) const;
  };
  
private:
  
  typedef vector<dboxes> cboxes; // ... for each component
  typedef vector<cboxes> rboxes; // ... for each refinement level
  typedef vector<rboxes> mboxes; // ... for each multigrid level
  
  
  
  void
  setup_bboxes ();
  
  static
  void
  optimise_field (dboxes & b,
                  iblistvect const dboxes::* field,
                  pvect dboxes::* fast_field);
  static
  void
  optimise_field (dboxes & b,
                  int proc,
                  iblist const dboxes::* field,
                  pvect dboxes::* fast_field);
  
  static
  void
  optimise2_field (cboxes & bs,
                   pvect const dboxes::* fast_field_recv,
                   pvect const dboxes::* fast_field_send,
                   srpvect dboxes::* fast_field_sendrecv);
  
public:                         // should be readonly
  
  // Fields
  gh & h;                       // hierarchy
  i2vect ghost_width;           // number of ghost zones
  i2vect buffer_width;          // number of buffer zones
  
  int prolongation_order_space; // order of spatial prolongation operator
  
  mboxes boxes;                 // grid hierarchy
  mboxes oldboxes;              // old grid hierarchy, used during regridding
  
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
  
  // Grid function management
  void add (ggf * f);
  void remove (ggf * f);
  
  // Output
  ostream & output (ostream & os) const;
};



inline ostream & operator<< (ostream & os, dh::dboxes const & b)
{
  return b.output (os);
}

inline ostream & operator<< (ostream & os, dh const & d)
{
  return d.output (os);
}



#endif // DH_HH
