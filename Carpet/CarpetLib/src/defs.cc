/***************************************************************************
                          defs.cc  -  Commonly used definitions
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/defs.cc,v 1.4 2001/03/24 22:38:48 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <assert.h>

#include <iostream>
#include <list>
#include <set>
#include <vector>

#if !defined(TMPL_IMPLICIT) || !defined(DEFS_HH)
#  include "defs.hh"
#endif

using namespace std;



// List output
template<class T>
ostream& operator<< (ostream& os, const list<T>& l) {
  os << "[";
  for (list<T>::const_iterator ti=l.begin(); ti!=l.end(); ++ti) {
    if (ti!=l.begin()) os << ",";
    os << *ti;
  }
  os << "]";
  return os;
}

// Set output
template<class T>
ostream& operator<< (ostream& os, const set<T>& s) {
  os << "{";
  for (set<T>::const_iterator ti=s.begin(); ti!=s.end(); ++ti) {
    if (ti!=s.begin()) os << ",";
    os << *ti;
  }
  os << "}";
  return os;
}

// Vector output
template<class T>
ostream& operator<< (ostream& os, const vector<T>& v) {
  os << "[";
  int cnt=0;
  for (vector<T>::const_iterator ti=v.begin(); ti!=v.end(); ++ti) {
    if (ti!=v.begin()) os << ",";
    os << cnt++ << ":" << *ti;
  }
  os << "]";
  return os;
}



#if defined(TMPL_EXPLICIT)
#include "bbox.hh"
#include "bboxset.hh"

template ostream& operator<< (ostream& os, const list<bbox<int,3> >& l);
template ostream& operator<< (ostream& os, const set<bbox<int,3> >& s);
template ostream& operator<< (ostream& os, const set<bboxset<int,3> >& s);
template ostream& operator<< (ostream& os, const vector<int>& v);
template ostream& operator<< (ostream& os, const vector<list<bbox<int,3> > >& v);
template ostream& operator<< (ostream& os, const vector<vector<bbox<int,3> > >& v);
template ostream& operator<< (ostream& os, const vector<vector<vector<bbox<int,3> > > >& v);
#endif
