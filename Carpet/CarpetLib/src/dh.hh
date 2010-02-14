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

#define CARPET_HAVE_BUFFER_WIDTHS

// Forward declaration
class ggf;
class dh;



// A data hierarchy (grid hierarchy plus ghost zones)
class dh {
  
  static list<dh*> alldh;
  list<dh*>::iterator alldhi;
  
  // Types
public:
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect; // vector of lists
  
  typedef vector <pseudoregion_t> pvect;
  typedef vector <sendrecv_pseudoregion_t> srpvect;
  
  
  
  struct dboxes {
    
    // Region description:
    
    ibbox exterior;             // whole region (including boundaries)
    ibbox owned;                // evolved in time
    ibbox interior;             // interior (without ghost zones)
    
    // Region statistics:
    typedef ibbox::size_type size_type;
    size_type exterior_size, owned_size, active_size;
    
    size_t memory () const CCTK_ATTRIBUTE_PURE;
    istream & input (istream & is);
    ostream & output (ostream & os) const;
  };
  
  struct full_dboxes {
    
    // Complete region description:
    
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
    
    bool operator== (full_dboxes const & b) const;
    bool operator!= (full_dboxes const & b) const
    {
      return not operator==(b);
    }
    
    size_t memory () const CCTK_ATTRIBUTE_PURE;
    istream & input (istream& is);
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
    
    bool do_init;               // the srpvects below are only defined
                                // if this is true
    srpvect fast_old2new_sync_sendrecv;
    srpvect fast_old2new_ref_prol_sendrecv;
    
    bool operator== (fast_dboxes const & b) const CCTK_ATTRIBUTE_PURE;
    bool operator!= (fast_dboxes const & b) const
    {
      return not operator==(b);
    }
    
    size_t memory () const CCTK_ATTRIBUTE_PURE;
    istream & input (istream & is);
    ostream & output (ostream & os) const;
  };
  
  typedef vector<dboxes> cboxes; // ... for each component
  typedef vector<cboxes> rboxes; // ... for each refinement level
  typedef vector<rboxes> mboxes; // ... for each multigrid level
  
  typedef vector<full_dboxes> full_cboxes; // ... for each component
  typedef vector<full_cboxes> full_rboxes; // ... for each refinement level
  typedef vector<full_rboxes> full_mboxes; // ... for each multigrid level
  
  typedef vector<fast_dboxes> fast_rboxes; // ... for each refinement level
  typedef vector<fast_rboxes> fast_mboxes; // ... for each multigrid level
  
private:
  
  void
  setup_bboxes ();
  
public:                         // should be readonly
  
  // Fields
  gh & h;                       // hierarchy
  gh::dh_handle gh_handle;
  
#if 0
  i2vect ghost_width;           // number of ghost zones
  i2vect buffer_width;          // number of buffer zones
  int prolongation_order_space; // order of spatial prolongation operator
#endif
  vector<i2vect> ghost_widths;  // number of ghost zones [rl]
  vector<i2vect> buffer_widths; // number of buffer zones [rl]
  vector<int> prolongation_orders_space; // order of spatial
                                         // prolongation operator [rl]
  
  mboxes boxes;                 // grid hierarchy
  fast_mboxes fast_boxes;       // grid hierarchy
  
  typedef list<ggf*>::iterator ggf_handle;
  list<ggf*> gfs;               // list of all grid functions
  
public:
  
  // Constructors
  dh (gh & h,
      vector<i2vect> const & ghost_widths, vector<i2vect> const & buffer_widths,
      vector<int> const & prolongation_orders_space);
  
  // Destructors
  ~dh ();
  
  // Helpers
  int prolongation_stencil_size (int rl) const CCTK_ATTRIBUTE_CONST;
  
  // Modifiers
  void regrid (bool do_init);
  void regrid_free (bool do_init);
  void recompose (int rl, bool do_prolongate);
  
private:
  int this_proc (int rl, int c) const CCTK_ATTRIBUTE_PURE;
  bool on_this_proc (int rl, int c) const CCTK_ATTRIBUTE_PURE;
  bool on_this_proc (int rl, int c, int cc) const CCTK_ATTRIBUTE_PURE;
  int this_oldproc (int rl, int c) const CCTK_ATTRIBUTE_PURE;
  bool on_this_oldproc (int rl, int c) const CCTK_ATTRIBUTE_PURE;
  
  static
  void
  broadcast_schedule (vector<fast_dboxes> & fast_level_otherprocs,
                      fast_dboxes & fast_level,
                      srpvect fast_dboxes::* const schedule_item);
  
public:
  // Grid function management
  ggf_handle add (ggf * f);
  void erase (ggf_handle fi);
  
  // Output
  size_t memory () const CCTK_ATTRIBUTE_PURE;
  static size_t allmemory () CCTK_ATTRIBUTE_PURE;
  ostream & output (ostream & os) const;
};



MPI_Datatype mpi_datatype (dh::dboxes const &) CCTK_ATTRIBUTE_CONST;
MPI_Datatype mpi_datatype (dh::fast_dboxes const &);
namespace dist {
  template<> inline MPI_Datatype mpi_datatype<dh::dboxes> ()
  CCTK_ATTRIBUTE_CONST;
  template<> inline MPI_Datatype mpi_datatype<dh::dboxes> ()
  { dh::dboxes dummy; return mpi_datatype(dummy); }
  template<> inline MPI_Datatype mpi_datatype<dh::fast_dboxes> ()
  CCTK_ATTRIBUTE_CONST;
  template<> inline MPI_Datatype mpi_datatype<dh::fast_dboxes> ()
  { dh::fast_dboxes dummy; return mpi_datatype(dummy); }
}

inline size_t memoryof (dh::dboxes const & b) CCTK_ATTRIBUTE_PURE;
inline size_t memoryof (dh::dboxes const & b)
{
  return b.memory ();
}

inline size_t memoryof (dh::full_dboxes const & b) CCTK_ATTRIBUTE_PURE;
inline size_t memoryof (dh::full_dboxes const & b)
{
  return b.memory ();
}

inline size_t memoryof (dh::fast_dboxes const & b) CCTK_ATTRIBUTE_PURE;
inline size_t memoryof (dh::fast_dboxes const & b)
{
  return b.memory ();
}

inline size_t memoryof (dh const & d) CCTK_ATTRIBUTE_PURE;
inline size_t memoryof (dh const & d)
{
  return d.memory ();
}

inline istream & operator>> (istream & is, dh::dboxes & b)
{
  return b.input (is);
}

inline istream & operator>> (istream & is, dh::full_dboxes & b)
{
  return b.input (is);
}

inline istream & operator>> (istream & is, dh::fast_dboxes & b)
{
  return b.input (is);
}

inline ostream & operator<< (ostream & os, dh::dboxes const & b)
{
  return b.output (os);
}

inline ostream & operator<< (ostream & os, dh::full_dboxes const & b)
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
