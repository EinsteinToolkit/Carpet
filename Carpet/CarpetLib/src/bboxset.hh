/***************************************************************************
                          bboxset.hh  -  Sets of bounding boxes
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/bboxset.hh,v 1.3 2001/03/12 16:54:25 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BBOXSET_HH
#define BBOXSET_HH

#include <cassert>
#include <iostream>
#include <set>

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"



// Forward definition
template<class T, int D> class bboxset;

template<class T,int D>
bboxset<T,D> operator- (const bbox<T,D>& b1, const bbox<T,D>& b2);
template<class T,int D>
bboxset<T,D> operator- (const bbox<T,D>& b, const bboxset<T,D>& s);

// Output
template<class T,int D>
ostream& operator<< (ostream& os, const bboxset<T,D>& s);



// Bounding box class
template<class T, int D>
class bboxset {
  
  // Types
  typedef bbox<T,D> box;
  typedef set<box> bset;
  
  // Fields
  bset bs;
  // Invariant:
  // All bboxes have the same stride.
  // No bbox is empty.
  // The bboxes don't overlap.
  
public:
  
  // Constructors
  bboxset ();
  bboxset (const box& b);
  bboxset (const bboxset& s);
  bboxset (const bset& bs);
  
  // Invariant
  bool invariant () const;
  
  // Normalisation
  void normalize ();
  
  // Accessors
  bool empty () const { return bs.empty(); }
  T size () const;
  
  // Add (bboxes that don't overlap)
  bboxset& operator+= (const box& b);
  bboxset& operator+= (const bboxset& s);
  bboxset operator+ (const box& b) const;
  bboxset operator+ (const bboxset& s) const;
  
  // Union
  bboxset& operator|= (const box& b);
  bboxset& operator|= (const bboxset& s);
  bboxset operator| (const box& b) const;
  bboxset operator| (const bboxset& s) const;
  
  // Intersection
  bboxset operator& (const box& b) const;
  bboxset operator& (const bboxset& s) const;
  bboxset& operator&= (const box& b);
  bboxset& operator&= (const bboxset& s);
  
  // Difference
  // friend bboxset operator- <T,D>(const box& b1, const box& b2);
  bboxset operator- (const box& b) const;
  bboxset& operator-= (const box& b);
  bboxset& operator-= (const bboxset& s);
  bboxset operator- (const bboxset& s) const;
  // friend bboxset operator- <T,D>(const box& b, const bboxset& s);
  
  // Iterators
  typedef bset::const_iterator const_iterator;
  typedef bset::iterator       iterator;
  
  iterator begin () const { return bs.begin(); }
  iterator end () const   { return bs.end(); }
  
  // Output
  friend ostream& operator<< <>(ostream& os, const bboxset& s);
};



#if defined(TMPL_IMPLICIT)
#  include "bboxset.cc"
#endif

#endif // BBOXSET_HH
