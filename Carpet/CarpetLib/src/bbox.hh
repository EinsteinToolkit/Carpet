/***************************************************************************
                          bbox.hh  -  Bounding boxes
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/bbox.hh,v 1.6 2001/03/14 11:00:26 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BBOX_HH
#define BBOX_HH

#include <iostream>

#include "defs.hh"
#include "vect.hh"



// Forward definition
template<class T, int D> class bbox;

// Output
template<class T,int D>
ostream& operator<< (ostream& os, const bbox<T,D>& b);



// Bounding box class
template<class T, int D>
class bbox {
  
  // Fields
  vect<T,D> _lower, _upper, _stride; // bounds are inclusive
  
public:
  
  // Constructors
  bbox ();
  bbox (const bbox& b);
  bbox& operator= (const bbox& b);
  bbox (const vect<T,D>& lower, const vect<T,D>& upper,
	const vect<T,D>& stride);
  
  // Accessors
  // (Don't return references; *this might be a temporary)
  vect<T,D> lower () const { return _lower; }
  vect<T,D> upper () const { return _upper; }
  vect<T,D> stride () const { return _stride; }
  vect<T,D> shape () const { return _upper - _lower + _stride; }
  
  bool empty() const {
    return any(lower()>upper());
  }
  
  T size () const;
  T num_points () const;
  
  // Queries
  bool contains (const vect<T,D>& x) const;
  
  // Operators
  bool operator== (const bbox& b) const;
  bool operator!= (const bbox& b) const;
  bool operator< (const bbox& b) const;
  bool operator> (const bbox& b) const;
  bool operator<= (const bbox& b) const;
  bool operator>= (const bbox& b) const;
  
  // Intersection
  bbox operator& (const bbox& b) const;
  
  // Containment
  bool is_contained_in (const bbox& b) const;
  
  // Alignment check
  bool is_aligned_with (const bbox& b) const;
  
  // Expand the bbox a little by multiples of the stride
  bbox expand (const vect<T,D>& lo, const vect<T,D>& hi) const;
  
  // Find the smallest b-compatible box around *this
  bbox expanded_for (const bbox& b) const;
  
  // Find the largest b-compatible box inside *this
  bbox contracted_for (const bbox& b) const;
  
  // Smallest bbox containing both boxes
  bbox expanded_containing (const bbox<T,D>& b) const;
  
  // Iterators
  class iterator {
  protected:
    bbox box;
    vect<T,D> pos;
  public:
    iterator (const bbox& box, const vect<T,D>& pos);
    const vect<T,D>& operator* () const { return pos; }
    bool operator!= (const iterator& i) const;
    iterator& operator++ ();
  };
  
  iterator begin () const;
  iterator end () const;
  
  // Output
  friend ostream& operator<< <>(ostream& os, const bbox& b);
};



#if defined(TMPL_IMPLICIT)
#  include "bbox.cc"
#endif

#endif // BBOX_HH
