#ifndef GGF_HH
#define GGF_HH

#include <cassert>
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
class ggf;

// Output
ostream& operator<< (ostream& os, const ggf& f);



// A generic grid function without type information
class ggf {

  // Types
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect;
  
  typedef gdata*        tdata;  // data ...
  typedef vector<tdata> fdata;  // ... for each time level
  typedef vector<fdata> cdata;  // ... for each component
  typedef vector<cdata> rdata;  // ... for each refinement level
  typedef vector<rdata> mdata;  // ... for each multigrid level
  
public:				// should be readonly
  
  // Fields
  const int varindex;           // Cactus variable index
  const operator_type transport_operator;
  
  const th &t;                  // time hierarchy
  const int prolongation_order_time; // order of temporal prolongation operator
  
  const gh &h;                  // grid hierarchy
  dh &d;			// data hierarchy

private:
  int timelevels_;              // time levels

protected:
  mdata storage;		// storage
  
public:
  const int vectorlength;       // vector length
  const int vectorindex;        // index of *this
  const ggf* vectorleader;      // first vector element
  
private:
  mdata oldstorage;            // temporary storage
  
public:

  // Constructors
  ggf (const int varindex, const operator_type transport_operator,
       th& t, dh& d,
       const int prolongation_order_time,
       const int vectorlength, const int vectorindex,
       ggf* const vectorleader);

  // Destructors
  virtual ~ggf ();

  // Comparison
  bool operator== (const ggf& f) const;
  
  // Querying
  int timelevels () const
  {
    return timelevels_;
  }
  
  
  
  // Modifiers
  void set_timelevels (int ml, int rl, int new_timelevels);

  void recompose_crop ();
  void recompose_allocate (int rl);
  void recompose_fill (comm_state& state, int rl, bool do_prolongate);
  void recompose_free (int rl);
  void recompose_bnd_prolongate (comm_state& state, int rl, bool do_prolongate);
  void recompose_sync (comm_state& state, int rl, bool do_prolongate);

  // Cycle the time levels by rotating the data sets
  void cycle (int rl, int c, int ml);
  
  // Flip the time levels by exchanging the data sets
  void flip (int rl, int c, int ml);
  
  
  
  // Helpers
  
protected:
  
  virtual gdata* typed_data (int tl, int rl, int c, int ml) = 0;
  
  
  
  // Operations
  
protected:
  
  // Copy a region
  void copycat (comm_state& state,
                int tl1, int rl1, int c1, int ml1,
		const ibbox dh::dboxes::* recv_list,
		int tl2, int rl2, int ml2,
		const ibbox dh::dboxes::* send_list);

  // Copy regions
  void copycat (comm_state& state,
                int tl1, int rl1, int c1, int ml1,
		const iblist dh::dboxes::* recv_list,
		int tl2, int rl2, int ml2,
		const iblist dh::dboxes::* send_list);

  // Copy regions
  void copycat (comm_state& state,
                int tl1, int rl1, int c1, int ml1,
		const iblistvect dh::dboxes::* recv_listvect,
		int tl2, int rl2, int ml2,
		const iblistvect dh::dboxes::* send_listvect);
  
  // Interpolate a region
  void intercat (comm_state& state,
                 int tl1, int rl1, int c1, int ml1,
		 const ibbox dh::dboxes::* recv_list,
		 const vector<int> tl2s, int rl2, int ml2,
		 const ibbox dh::dboxes::* send_list,
		 CCTK_REAL time);

  // Interpolate regions
  void intercat (comm_state& state,
                 int tl1, int rl1, int c1, int ml1,
		 const iblist dh::dboxes::* recv_list,
		 const vector<int> tl2s, int rl2, int ml2,
		 const iblist dh::dboxes::* send_list,
		 CCTK_REAL time);

  // Interpolate regions
  void intercat (comm_state& state,
                 int tl1, int rl1, int c1, int ml1,
		 const iblistvect dh::dboxes::* recv_listvect,
		 const vector<int> tl2s, int rl2, int ml2,
		 const iblistvect dh::dboxes::* send_listvect,
		 CCTK_REAL time);



public:

  // The grid boundaries have to be updated after calling mg_restrict,
  // mg_prolongate, ref_restrict, or ref_prolongate.

  // "Updating" means here that the boundaries have to be
  // synchronised.  They don't need to be prolongated.

  // Copy a component from the next time level
  void copy (comm_state& state, int tl, int rl, int c, int ml);

  // Synchronise the boundaries of a component
  void sync (comm_state& state, int tl, int rl, int c, int ml);

  // Prolongate the boundaries of a component
  void ref_bnd_prolongate (comm_state& state,
                           int tl, int rl, int c, int ml, CCTK_REAL time);

  // Restrict a multigrid level
  void mg_restrict (comm_state& state,
                    int tl, int rl, int c, int ml, CCTK_REAL time);

  // Prolongate a multigrid level
  void mg_prolongate (comm_state& state,
                      int tl, int rl, int c, int ml, CCTK_REAL time);

  // Restrict a refinement level
  void ref_restrict (comm_state& state,
                     int tl, int rl, int c, int ml, CCTK_REAL time);

  // Prolongate a refinement level
  void ref_prolongate (comm_state& state,
                       int tl, int rl, int c, int ml, CCTK_REAL time);
  
  
  
  // Access to the data
  virtual const gdata* operator() (int tl, int rl, int c, int ml) const = 0;
  
  virtual gdata* operator() (int tl, int rl, int c, int ml) = 0;
  
  
  
  // Output
  virtual ostream& output (ostream& os) const = 0;

private:
  ggf (); // canonical default construtor
  ggf ( const ggf & ); // canonical copy construtor
  ggf & operator =( const ggf & ); // canonical copy

};



template<int D>
inline ostream& operator<< (ostream& os, const ggf& f) {
  return f.output(os);
}



#endif // GGF_HH
