// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/defs.cc,v 1.14 2003/02/27 22:23:19 shawley Exp $

#include <assert.h>
#include <ctype.h>

#include <iostream>
#include <list>
#include <set>
#include <vector>

#include "defs.hh"

using namespace std;



void skipws (istream& is) {
  while (is.good() && isspace(is.peek())) {
    is.get();
  }
}



// Vector input
template<class T>
istream& input (istream& is, vector<T>& v) {
  v.clear();
  skipws (is);
  assert (is.peek() == '[');
  is.get();
  skipws (is);
  while (is.good() && is.peek() != ']') {
    T elem;
    is >> elem;
    v.push_back (elem);
    skipws (is);
    if (is.peek() != ',') break;
    is.get();
    skipws (is);
  }
  skipws (is);
  assert (is.peek() == ']');
  is.get();
  return is;
}



// List output
template<class T>
ostream& output (ostream& os, const list<T>& l) {
  os << "[";
  for (typename list<T>::const_iterator ti=l.begin(); ti!=l.end(); ++ti) {
    if (ti!=l.begin()) os << ",";
    os << *ti;
  }
  os << "]";
  return os;
}

// Set output
template<class T>
ostream& output (ostream& os, const set<T>& s) {
  os << "{";
  for (typename set<T>::const_iterator ti=s.begin(); ti!=s.end(); ++ti) {
    if (ti!=s.begin()) os << ",";
    os << *ti;
  }
  os << "}";
  return os;
}

// Vector output
template<class T>
ostream& output (ostream& os, const vector<T>& v) {
  os << "[";
  int cnt=0;
  for (typename vector<T>::const_iterator ti=v.begin(); ti!=v.end(); ++ti) {
    if (ti!=v.begin()) os << ",";
    os << cnt++ << ":" << *ti;
  }
  os << "]";
  return os;
}



#include "bbox.hh"
#include "bboxset.hh"

template istream& input (istream& os, vector<bbox<int,3> >& v);
template istream& input (istream& os, vector<vector<bbox<int,3> > >& v);
template istream& input (istream& os, vector<vector<vect<vect<bool,2>,3> > >& v);

template ostream& output (ostream& os, const list<bbox<int,3> >& l);
template ostream& output (ostream& os, const set<bbox<int,3> >& s);
template ostream& output (ostream& os, const set<bboxset<int,3> >& s);
template ostream& output (ostream& os, const vector<int>& v);
template ostream& output (ostream& os, const vector<bbox<int,3> >& v);
template ostream& output (ostream& os, const vector<list<bbox<int,3> > >& v);
template ostream& output (ostream& os, const vector<vector<int> >& v);
template ostream& output (ostream& os, const vector<vector<bbox<int,3> > >& v);
template ostream& output (ostream& os, const vector<vector<vect<vect<bool,2>,3> > >& v);
template ostream& output (ostream& os, const vector<vect<vect<bool,2>,3> >& v);
template ostream& output (ostream& os, const vector<vector<vector<bbox<int,3> > > >& v);
