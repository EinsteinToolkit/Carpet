/***************************************************************************
                          gdata.hh  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.hh,v 1.4 2001/03/14 11:00:26 eschnett Exp $

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

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "vect.hh"



// Forward declaration
template<int D> class generic_data;

// Output
template<int D>
ostream& operator<< (ostream& os, const generic_data<D>* f);



// A generic data storage without type information
template<int D>
class generic_data {

  // Types
  typedef vect<int,D> ivect;
  typedef bbox<int,D> ibbox;

protected:                      // should be readonly

  // Fields
  bool _has_storage;		// has storage associated (on some processor)
  bool _owns_storage;		// owns the storage
  ivect _shape, _stride;      	// shape and index order
  int _size;			// size

  int _proc;			// stored on processor

  ibbox _extent;		// bbox for all data

public:

  // Constructors
  generic_data ();

  // Destructors
  virtual ~generic_data ();

  // Pseudo constructors
  virtual generic_data* make_typed () const = 0;

  // Storage management
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const mem=0) = 0;
  virtual void free () = 0;
  virtual void transfer_from (generic_data* src) = 0;

  // Processor management
  virtual void change_processor (const int newproc, void* const mem=0) = 0;

  // Accessors
  bool has_storage () const {
    return _has_storage;
  }
  bool owns_storage () const {
    assert (_has_storage);
    return _owns_storage;
  }
  
  virtual const void* storage () const = 0;
  
  virtual void* storage () = 0;
  
  const ivect& shape () const {
    assert (_has_storage);
    return _shape;
  }

  const ivect& stride () const {
    assert (_has_storage);
    return _stride;
  }

  int size () const {
    assert (_has_storage);
    return _size;
  }

  int proc () const {
    assert (_has_storage);
    return _proc;
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
  virtual void copy_from (const generic_data* src,
                          const ibbox& b) = 0;
  virtual void interpolate_from (const generic_data* src,
				 const ibbox& box) = 0;
  virtual void interpolate_from (const generic_data* src, const double sfact,
				 const generic_data* trc, const double tfact,
				 const ibbox& box) = 0;
  
  // Output
  template<int DD>
  void write_ascii (const string name, const int time,
                    const vect<int,D>& org, const vect<int,DD>& dirs,
		    const int tl, const int rl,
                    const int c, const int ml)
    const;
protected:
  virtual void write_ascii_output_element (ofstream& file, const ivect& index)
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

  // Output
  friend ostream& operator<< <>(ostream& os, const generic_data* d);

  virtual ostream& out (ostream& os) const = 0;
};



#if defined(TMPL_IMPLICIT)
#  include "gdata.cc"
#endif

#endif // GDATA_HH
