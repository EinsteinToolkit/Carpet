/***************************************************************************
                          data.hh  -  Data storage
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/data.hh,v 1.4 2001/03/22 18:42:05 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DATA_HH
#define DATA_HH

#include <assert.h>

#include <string>

#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "gdata.hh"
#include "vect.hh"

using namespace std;



// Forward definition
template<class T,int D> class data;

// Output
template<class T,int D>
ostream& operator<< (ostream& os, const data<T,D>& d);



// A real data storage
template<class T,int D>
class data: public generic_data<D> {
  
  // Types
  typedef vect<int,D> ivect;
  typedef bbox<int,D> ibbox;

  // Fields
  T* _storage;			// the data (if located on this processor)

public:
  
  // Constructors
  data ();
  data (const ibbox& extent, const int proc);

  // Destructors
  virtual ~data ();

  // Pseudo constructors
  virtual data* make_typed () const;

  // Storage management
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const mem=0);
  virtual void free ();
  virtual void transfer_from (generic_data<D>* gsrc);

  // Processor management
  virtual void change_processor (const int newproc, void* const mem=0);

  // Accessors
  virtual const void* storage () const {
    assert (_has_storage);
    return _storage;
  }

  virtual void* storage () {
    assert (_has_storage);
    return _storage;
  }
  
  // Data accessors
  const T& operator[] (const ivect& index) const {
    assert (_storage);
    return _storage[offset(index)];
  }
  
  T& operator[] (const ivect& index) {
    assert (_storage);
    return _storage[offset(index)];
  }
  
  // Data manipulators
  void copy_from_innerloop (const generic_data<D>* src,
			    const ibbox& box);
  void interpolate_from_innerloop (const generic_data<D>* src,
				   const ibbox& box);
  void interpolate_from_innerloop (const generic_data<D>* src1, const int t1,
				   const generic_data<D>* src2, const int t2,
				   const ibbox& box, const int t);
  
  void write_ascii_output_element (ofstream& file, const ivect& index) const;
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
  ostream& out (ostream& os) const;
};



#if defined(TMPL_IMPLICIT)
#  include "data.cc"
#endif

#endif // DATA_HH
