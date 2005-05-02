#include <cassert>
#include <cctype>
#include <iostream>
#include <list>
#include <set>
#include <stack>
#include <vector>

#include "cctk.h"

#include "defs.hh"

using namespace std;



template <typename T>
inline T ipow_helper (T x, unsigned int y)
{
  T z = y&1 ? x : 1;
  while (y >>= 1)
  {
    x *= x;
    if (y & 1) z *= x;
  }
  return z;
}

template<class T>
T ipow (T x, int y)
{
  if (y < 0)
    return T(1) / ipow_helper(x, -y);
  else
    return ipow_helper(x, y);
}



void skipws (istream& is) {
  while (is.good() && isspace(is.peek())) {
    is.get();
  }
}



void expect (istream& is, const char c) {
  if (is.peek() == c) return;
  cout << "While reading characters from a stream:" << endl
       << "   Character '" << c << "' expected, but not found." << endl
       << "   The next up to 100 available characters are \"";
  for (int i=0; i<100; ++i) {
    const int uc = is.get();
    if (uc<0) break;
    cout << (unsigned char)uc;
  }
  cout << "\"." << endl;
  throw input_error();
}



void consume (istream& is, const char c) {
  expect (is, c);
  is.get();
}



// Vector input
template<class T>
istream& input (istream& is, vector<T>& v) {
  v.clear();
  try {
    skipws (is);
    consume (is, '[');
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
    consume (is, ']');
  } catch (input_error &err) {
    cout << "Input error while reading a vector<>" << endl
         << "   The following elements have been read so far: " << v << endl;
    throw err;
  }
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

// Stack output
template<class T>
ostream& output (ostream& os, const stack<T>& s) {
  stack<T> s2 (s);
  list<T> l;
  while (! s2.empty()) {
    l.insert (l.begin(), s2.top());
    s2.pop();
  }
  return output (os, l);
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
#include "vect.hh"

template int ipow (int x, int y);
template CCTK_REAL ipow (CCTK_REAL x, int y);
template vect<int,3> ipow (vect<int,3> x, int y);

template istream& input (istream& os, vector<int>& v);
template istream& input (istream& os, vector<bbox<int,3> >& v);
template istream& input (istream& os, vector<bbox<CCTK_REAL,3> >& v);
template istream& input (istream& os, vector<vector<bbox<int,3> > >& v);
template istream& input (istream& os, vector<vector<bbox<CCTK_REAL,3> > >& v);
template istream& input (istream& os, vector<vect<int,3> >& v);
template istream& input (istream& os, vector<vect<vect<bool,2>,3> >& v);
template istream& input (istream& os, vector<vector<vect<vect<bool,2>,3> > >& v);

template ostream& output (ostream& os, const list<bbox<int,3> >& l);
template ostream& output (ostream& os, const set<bbox<int,3> >& s);
template ostream& output (ostream& os, const set<bboxset<int,3> >& s);
template ostream& output (ostream& os, const stack<bbox<int,3> >& s);
template ostream& output (ostream& os, const vector<bool>& v);
template ostream& output (ostream& os, const vector<int>& v);
template ostream& output (ostream& os, const vector<CCTK_REAL>& v);
template ostream& output (ostream& os, const vector<bbox<int,3> >& v);
template ostream& output (ostream& os, const vector<bbox<CCTK_REAL,3> >& v);
template ostream& output (ostream& os, const vector<list<bbox<int,3> > >& v);
template ostream& output (ostream& os, const vector<vector<int> >& v);
template ostream& output (ostream& os, const vector<vector<CCTK_REAL> >& v);
template ostream& output (ostream& os, const vector<vector<bbox<int,3> > >& v);
template ostream& output (ostream& os, const vector<vector<bbox<CCTK_REAL,3> > >& v);
template ostream& output (ostream& os, const vector<vector<vect<vect<bool,2>,3> > >& v);
template ostream& output (ostream& os, const vector<vect<int,3> >& v);
template ostream& output (ostream& os, const vector<vect<vect<bool,2>,3> >& v);
template ostream& output (ostream& os, const vector<vector<vector<bbox<int,3> > > >& v);
