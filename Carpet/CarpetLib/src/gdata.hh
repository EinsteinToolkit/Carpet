/***************************************************************************
                          gdata.hh  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.hh,v 1.12 2002/01/08 12:03:55 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GDATA_HH
#define GDATA_HH

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <string>

#include "defs.hh"
#include "dgdata.hh"
#include "dist.hh"
#include "bbox.hh"
#include "vect.hh"

using namespace std;



// A generic data storage without type information
template<int D>
class generic_data: public dimgeneric_data {

  // Types
  typedef vect<int,D> ivect;
  typedef bbox<int,D> ibbox;

protected:                      // should be readonly

  // Fields
  ivect _shape, _stride;      	// shape and index order
  
  ibbox _extent;		// bbox for all data

public:

  // Constructors
  generic_data ();

  // Destructors
  virtual ~generic_data ();

  // Pseudo constructors
  virtual generic_data<D>* make_typed () const = 0;

  // Storage management
  virtual void transfer_from (generic_data<D>* src) = 0;
  
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const mem=0) = 0;
  virtual void free () = 0;
  
  // Accessors
  
  const ivect& shape () const {
    assert (_has_storage);
    return _shape;
  }

  const ivect& stride () const {
    assert (_has_storage);
    return _stride;
  }
  
  const ibbox& extent () const {
    assert (_has_storage);
    return _extent;
  }

  // Data accessors
  int offset (const ivect& index) const {
    assert (_has_storage);
    assert (all((index-extent().lower()) % extent().stride() == 0));
    ivect ind = (index-extent().lower()) / extent().stride();
    assert (all(ind>=0 && ind<=shape()));
    return dot(ind, stride());
  }

  // Data manipulators
  void copy_from (const generic_data* src, const ibbox& box);
  void interpolate_from (const vector<const generic_data*> srcs,
			 const vector<int> tls,
			 const ibbox& box, const int tl,
			 const int order_space,
			 const int order_time);
protected:
  virtual void
  copy_from_innerloop (const generic_data* src, const ibbox& box) = 0;
  virtual void
  interpolate_from_innerloop (const vector<const generic_data*> srcs,
			      const vector<int> tls,
			      const ibbox& box, const int tl,
			      const int order_space,
			      const int order_time) = 0;
  
public:
  
  // Output
  template<int DD>
  void write_ascii (ostream& os, const int time,
                    const vect<int,D>& org, const vect<int,DD>& dirs,
		    const int tl, const int rl,
                    const int c, const int ml,
		    const double ctime,
		    const vect<double,D>& coord_lower,
		    const vect<double,D>& coord_upper)
    const;
protected:
  virtual void write_ascii_output_element (ostream& os, const ivect& index)
    const = 0;
public:
//   void write_ieee (const string name, const int time,
// 		   const int tl, const int rl, const int c, const int ml)
//     const;
//   void write_hdf (const string name, const int time,
// 		  const int tl, const int rl, const int c, const int ml)
//     const;
//   void write_h5 (const string name, const int time,
// 		 const int tl, const int rl, const int c, const int ml)
//     const;
public:
};



#if defined(TMPL_IMPLICIT)
#  include "gdata.cc"
#endif

#endif // GDATA_HH
