/***************************************************************************
                          gf.hh  -  Grid Function
                          data for every element of a data hierarchy
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gf.hh,v 1.4 2001/06/12 14:56:59 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GF_HH
#define GF_HH

#include <assert.h>
#include <math.h>

#include <iostream>
#include <string>

#include "bbox.hh"
#include "bboxset.hh"
#include "data.hh"
#include "defs.hh"
#include "dh.hh"
#include "ggf.hh"
#include "th.hh"
#include "vect.hh"

using namespace std;



// A real grid function
template<class T,int D>
class gf: public generic_gf<D> {
  
  // Types
  typedef vect<int,D>    ivect;
  typedef bbox<int,D>    ibbox;
  typedef bboxset<int,D> ibset;
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect;
  
  typedef data<T,D>*    tdata;	        // data ...
  typedef vector<tdata> mdata;	        // ... for each multigrid level
  typedef vector<mdata> cdata;	        // ... for each component
  typedef vector<cdata> rdata;	        // ... for each refinement level
  typedef vector<rdata> fdata;          // ... for each time level

public:
  
  // Constructors
  gf (const string name, th& t, dh<D>& d, const int tmin, const int tmax);
  
  // Destructors
  virtual ~gf ();
  
  
  
  // Helpers
  
protected:
  
  virtual generic_data<D>* typed_data() { return new data<T,D>; }
  
  
  
  // Access to the data
  
public:
  
  virtual const data<T,D>* operator() (int tl, int rl, int c, int ml) const;
  
  virtual data<T,D>* operator() (int tl, int rl, int c, int ml);
  
  
  
  // Output
  virtual ostream& output (ostream& os) const;
};



#if defined(TMPL_IMPLICIT)
#  include "gf.cc"
#endif

#endif // GF_HH
