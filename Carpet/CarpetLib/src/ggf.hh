/***************************************************************************
                          ggf.hh  -  Generic Grid Function
                          grid function without type information
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/ggf.hh,v 1.2 2001/03/15 09:59:43 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GGF_HH
#define GGF_HH

#include <cassert>
#include <iostream>
#include <string>

#include "defs.hh"
#include "dh.hh"
#include "gdata.hh"
#include "gh.hh"
#include "th.hh"



// Forward declaration
template<int D> class generic_gf;

// Output
template<int D>
ostream& operator<< (ostream& os, const generic_gf<D>& f);



// A generic grid function without type information
template<int D>
class generic_gf {

  // Types
  typedef vect<int,D>    ivect;
  typedef bbox<int,D>    ibbox;
  typedef bboxset<int,D> ibset;
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect;

  typedef generic_data<D>* tdata; // data ...
  typedef vector<tdata>    mdata; // ... for each multigrid level
  typedef vector<mdata>    cdata; // ... for each component
  typedef vector<cdata>    rdata; // ... for each refinement level
  typedef vector<rdata>    fdata; // ... for each time level

public:				// should be readonly

  // Fields
  string name;

  gh<D> &h;			// grid hierarchy
  th<D> &t;			// time hierarchy
  dh<D> &d;			// data hierarchy
  int tmin, tmax;		// timelevels

protected:
  fdata storage;		// storage

public:

  // Constructors
  generic_gf (const string name, th<D>& t, dh<D>& d,
              const int tmin, const int tmax);

  // Destructors
  virtual ~generic_gf ();

  // Comparison
  virtual bool operator== (const generic_gf<D>& f) const;



  // Modifiers
  virtual void recompose ();
  
  
  
  // Helpers
  
protected:
  
  virtual generic_data<D>* typed_data() = 0;
  
  
  
  // Operations
  
protected:
  
  // Copy region for a component (between time levels)
  virtual void copycat (int tl1, int rl1, int c1, int ml1,
			const ibbox dh<D>::dboxes::* recv_list,
			int tl2, int rl2, int ml2,
			const ibbox dh<D>::dboxes::* send_list);

  // Copy regions for a component (between multigrid levels)
  virtual void copycat (int tl1, int rl1, int c1, int ml1,
			const iblist dh<D>::dboxes::* recv_list,
			int tl2, int rl2, int ml2,
			const iblist dh<D>::dboxes::* send_list);

  // Copy regions for a level (between refinement levels)
  virtual void copycat (int tl1, int rl1, int c1, int ml1,
			const iblistvect dh<D>::dboxes::* recv_listvect,
			int tl2, int rl2, int ml2,
			const iblistvect dh<D>::dboxes::* send_listvect);
  
  // Interpolate a component (between time levels)
  virtual void intercat (int tl1, int rl1, int c1, int ml1,
			 const ibbox dh<D>::dboxes::* recv_list,
			 int tl2, const double fact2,
			 int tl3, const double fact3,
			 int rl2, int ml2,
			 const ibbox dh<D>::dboxes::* send_list);

  // Interpolate a component (between multigrid levels)
  virtual void intercat (int tl1, int rl1, int c1, int ml1,
			 const iblist dh<D>::dboxes::* recv_list,
			 int tl2, const double fact2,
			 int tl3, const double fact3,
			 int rl2, int ml2,
			 const iblist dh<D>::dboxes::* send_list);

  // Interpolate a level (between refinement levels)
  virtual void intercat (int tl1, int rl1, int c1, int ml1,
			 const iblistvect dh<D>::dboxes::* recv_listvect,
			 int tl2, const double fact2,
			 int tl3, const double fact3,
			 int rl2, int ml2,
			 const iblistvect dh<D>::dboxes::* send_listvect);



public:

  // The grid boundaries have to be updated after calling mg_restrict,
  // mg_prolongate, ref_restrict, or ref_prolongate.

  // "Updating" means here that the boundaries have to be
  // synchronised.  They don't need to be prolongated.

  // Copy a component from the next time level
  virtual void copy (int tl, int rl, int c, int ml);

  // Synchronise the boundaries of a component
  virtual void sync (int tl, int rl, int c, int ml);

  // Prolongate the boundaries of a component
  virtual void ref_bnd_prolongate (int tl, int rl, int c, int ml);

  // Restrict a multigrid level
  virtual void mg_restrict (int tl, int rl, int c, int ml);

  // Prolongate a multigrid level
  virtual void mg_prolongate (int tl, int rl, int c, int ml);

  // Restrict a refinement level
  virtual void ref_restrict (int tl, int rl, int c, int ml);

  // Prolongate a refinement level
  virtual void ref_prolongate (int tl, int rl, int c, int ml);
  
  
  
  // Access to the data
  virtual const generic_data<D>* operator() (int tl, int rl, int c, int ml)
    const = 0;
  
  virtual generic_data<D>* operator() (int tl, int rl, int c, int ml) = 0;
  
  
  
  // Output
  friend ostream& operator<< <> (ostream& os, const generic_gf& f);
  
  virtual ostream& out (ostream& os) const = 0;
};



#if defined(TMPL_IMPLICIT)
#  include "ggf.cc"
#endif

#endif // GGF_HH
