// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/ggf.hh,v 1.17 2003/11/05 16:18:39 schnetter Exp $

#ifndef GGF_HH
#define GGF_HH

#include <assert.h>

#include <iostream>
#include <string>
#include <vector>

#include "cctk.h"

#include "defs.hh"
#include "dh.hh"
#include "gdata.hh"
#include "gh.hh"
#include "th.hh"

using namespace std;



// Forward declaration
template<int D> class ggf;

// Output
template<int D>
ostream& operator<< (ostream& os, const ggf<D>& f);



// A generic grid function without type information
template<int D>
class ggf {

  // Types

  typedef vect<int,D>    ivect;
  typedef bbox<int,D>    ibbox;
  typedef bboxset<int,D> ibset;
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect;

  typedef gdata<D>*     tdata;  // data ...
  typedef vector<tdata> mdata;  // ... for each multigrid level
  typedef vector<mdata> cdata;  // ... for each component
  typedef vector<cdata> rdata;  // ... for each refinement level
  typedef vector<rdata> fdata;  // ... for each time level

public:				// should be readonly

  // Fields
  int varindex;                 // Cactus variable index
  
  th<D> &t;			// time hierarchy
  int tmin, tmax;		// timelevels
  int prolongation_order_time;	// order of temporal prolongation operator
  
  gh<D> &h;			// grid hierarchy
  dh<D> &d;			// data hierarchy

protected:
  fdata storage;		// storage

public:

  // Constructors
  ggf (const int varindex, th<D>& t, dh<D>& d,
       const int tmin, const int tmax,
       const int prolongation_order_time);

  // Destructors
  virtual ~ggf ();

  // Comparison
  bool operator== (const ggf<D>& f) const;



  // Modifiers
  // VGF
  void recompose (const int initialise_from, const bool do_prolongate);

  // Cycle the time levels by rotating the data sets
  void cycle (int rl, int c, int ml);
  
  // Flip the time levels by exchanging the data sets
  void flip (int rl, int c, int ml);
  
#if 0
  // Copy data from current time level to all previous time levels
  void copytoprevs (int rl, int c, int ml);
#endif
  
  
  
  // Helpers
  
protected:
  
  virtual gdata<D>* typed_data() = 0;
  
  
  
  // Operations
  
protected:
  
  // Copy a region
  void copycat (comm_state<D>& state,
                int tl1, int rl1, int c1, int ml1,
		const ibbox dh<D>::dboxes::* recv_list,
		int tl2, int rl2, int ml2,
		const ibbox dh<D>::dboxes::* send_list);

  // Copy regions
  void copycat (comm_state<D>& state,
                int tl1, int rl1, int c1, int ml1,
		const iblist dh<D>::dboxes::* recv_list,
		int tl2, int rl2, int ml2,
		const iblist dh<D>::dboxes::* send_list);

  // Copy regions
  void copycat (comm_state<D>& state,
                int tl1, int rl1, int c1, int ml1,
		const iblistvect dh<D>::dboxes::* recv_listvect,
		int tl2, int rl2, int ml2,
		const iblistvect dh<D>::dboxes::* send_listvect);
  
  // Interpolate a region
  void intercat (comm_state<D>& state,
                 int tl1, int rl1, int c1, int ml1,
		 const ibbox dh<D>::dboxes::* recv_list,
		 const vector<int> tl2s, int rl2, int ml2,
		 const ibbox dh<D>::dboxes::* send_list,
		 CCTK_REAL time);

  // Interpolate regions
  void intercat (comm_state<D>& state,
                 int tl1, int rl1, int c1, int ml1,
		 const iblist dh<D>::dboxes::* recv_list,
		 const vector<int> tl2s, int rl2, int ml2,
		 const iblist dh<D>::dboxes::* send_list,
		 CCTK_REAL time);

  // Interpolate regions
  void intercat (comm_state<D>& state,
                 int tl1, int rl1, int c1, int ml1,
		 const iblistvect dh<D>::dboxes::* recv_listvect,
		 const vector<int> tl2s, int rl2, int ml2,
		 const iblistvect dh<D>::dboxes::* send_listvect,
		 CCTK_REAL time);



public:

  // The grid boundaries have to be updated after calling mg_restrict,
  // mg_prolongate, ref_restrict, or ref_prolongate.

  // "Updating" means here that the boundaries have to be
  // synchronised.  They don't need to be prolongated.

  // Copy a component from the next time level
  void copy (comm_state<D>& state, int tl, int rl, int c, int ml);

  // Synchronise the boundaries of a component
  void sync (comm_state<D>& state, int tl, int rl, int c, int ml);

  // Prolongate the boundaries of a component
  void ref_bnd_prolongate (comm_state<D>& state, int tl, int rl, int c, int ml, CCTK_REAL time);

  // Restrict a multigrid level
  void mg_restrict (comm_state<D>& state, int tl, int rl, int c, int ml, CCTK_REAL time);

  // Prolongate a multigrid level
  void mg_prolongate (comm_state<D>& state, int tl, int rl, int c, int ml, CCTK_REAL time);

  // Restrict a refinement level
  void ref_restrict (comm_state<D>& state, int tl, int rl, int c, int ml, CCTK_REAL time);

  // Prolongate a refinement level
  void ref_prolongate (comm_state<D>& state, int tl, int rl, int c, int ml, CCTK_REAL time);
  
  
  
  // Access to the data
  virtual const gdata<D>* operator() (int tl, int rl, int c, int ml) const = 0;
  
  virtual gdata<D>* operator() (int tl, int rl, int c, int ml) = 0;
  
  
  
  // Output
  virtual ostream& output (ostream& os) const = 0;
};



template<int D>
inline ostream& operator<< (ostream& os, const ggf<D>& f) {
  return f.output(os);
}



#endif // GGF_HH
