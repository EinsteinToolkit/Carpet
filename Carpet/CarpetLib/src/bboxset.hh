#ifndef BBOXSET_HH
#define BBOXSET_HH

#include <cassert>
#include <iostream>
#include <list>
#include <set>
#include <vector>

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"

using namespace std;



// Forward declaration
template<class T, int D> class bboxset;

// template<class T,int D>
// bboxset<T,D> operator+ (const bbox<T,D>& b1, const bbox<T,D>& b2);
// template<class T,int D>
// bboxset<T,D> operator+ (const bbox<T,D>& b, const bboxset<T,D>& s);

// template<class T,int D>
// bboxset<T,D> operator- (const bbox<T,D>& b1, const bbox<T,D>& b2);
// template<class T,int D>
// bboxset<T,D> operator- (const bbox<T,D>& b, const bboxset<T,D>& s);

// Output
template<class T,int D>
ostream& operator<< (ostream& os, const bboxset<T,D>& s);



// Bounding box class
template<class T, int D>
class bboxset {
  
  // Types
  typedef bbox<T,D> box;
  //S typedef set<box> bset;
  typedef list<box> bset;
  
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
  
  bboxset (const list<box>& lb);
  bboxset (const vector<list<box> >& vlb);
  
  // Invariant
  bool invariant () const;
  
  // Normalisation
  void normalize ();
  
  // Accessors
  bool empty () const { return bs.empty(); }
  // T size () const;
  typedef typename box::size_type size_type;
  size_type size () const;
  int setsize () const { return bs.size(); }
  
  // Find out whether this bboxset intersects the bbox b
  bool intersects (const box& b) const;
  
  // Add (bboxes that don't overlap)
  bboxset& operator+= (const box& b)
  {
    if (b.empty()) return *this;
    // This is very slow when there are many bboxes
#if 0 && ! defined(CARPET_OPTIMISE)
    // check for overlap
    for (const_iterator bi=begin(); bi!=end(); ++bi) {
      assert (not (*bi).intersects(b));
    }
#endif
    //S bs.insert(b);
    bs.push_back(b);
    assert (invariant());
    return *this;
  }
  
  bboxset& operator+= (const bboxset& s);
  bboxset& add_transfer (bboxset& s);
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
  static bboxset minus (const box& b1, const box& b2);
  bboxset operator- (const box& b) const;
  bboxset& operator-= (const box& b)
  {
    *this = *this - b;
    assert (invariant());
    return *this;
  }
  bboxset& operator-= (const bboxset& s);
  bboxset operator- (const bboxset& s) const;
  // friend bboxset operator- <T,D>(const box& b, const bboxset& s);
  static bboxset minus (const box& b, const bboxset& s);
  
  // Equality
  bool operator== (const bboxset& s) const;
  bool operator!= (const bboxset& s) const;
  bool operator< (const bboxset& s) const;
  bool operator<= (const bboxset& s) const;
  bool operator> (const bboxset& s) const;
  bool operator>= (const bboxset& s) const;
  
  // Iterators
  typedef typename bset::const_iterator const_iterator;
  typedef typename bset::iterator       iterator;
  
  const_iterator begin () const { return bs.begin(); }
  const_iterator end () const   { return bs.end(); }
//   iterator begin () const { return bs.begin(); }
//   iterator end () const   { return bs.end(); }
  
  // Memory usage
  size_t memory () const { return memoryof (bs); }
  
  // Output
  void output (ostream& os) const;
};



template<class T,int D>
inline bboxset<T,D> operator+ (const bbox<T,D>& b1, const bbox<T,D>& b2) {
  return bboxset<T,D>(b1) + bboxset<T,D>(b2);
}

template<class T,int D>
inline bboxset<T,D> operator+ (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) + s;
}

template<class T,int D>
inline bboxset<T,D> operator- (const bbox<T,D>& b1, const bbox<T,D>& b2) {
  return bboxset<T,D>::minus(b1,b2);
}

template<class T,int D>
inline bboxset<T,D> operator- (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>::minus(b,s);
}



template<class T,int D>
inline bboxset<T,D> operator| (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return s | b;
}

template<class T,int D>
inline bboxset<T,D> operator& (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return s & b;
}



template<class T,int D>
inline bool operator== (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) == s;
}

template<class T,int D>
inline bool operator!= (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) != s;
}

template<class T,int D>
inline bool operator< (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) < s;
}

template<class T,int D>
inline bool operator<= (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) <= s;
}

template<class T,int D>
inline bool operator> (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) > s;
}

template<class T,int D>
inline bool operator>= (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) >= s;
}



template<class T,int D>
inline bool operator== (const bboxset<T,D>& s, const bbox<T,D>& b) {
  return s == bboxset<T,D>(b);
}

template<class T,int D>
inline bool operator!= (const bboxset<T,D>& s, const bbox<T,D>& b) {
  return s != bboxset<T,D>(b);
}

template<class T,int D>
inline bool operator< (const bboxset<T,D>& s, const bbox<T,D>& b) {
  return s < bboxset<T,D>(b);
}

template<class T,int D>
inline bool operator<= (const bboxset<T,D>& s, const bbox<T,D>& b) {
  return s <= bboxset<T,D>(b);
}

template<class T,int D>
inline bool operator> (const bboxset<T,D>& s, const bbox<T,D>& b) {
  return s > bboxset<T,D>(b);
}

template<class T,int D>
inline bool operator>= (const bboxset<T,D>& s, const bbox<T,D>& b) {
  return s >= bboxset<T,D>(b);
}



// Memory usage
template<class T, int D>
inline size_t memoryof (bboxset<T,D> const & s) { return s.memory(); }



// Output
template<class T,int D>
inline ostream& operator<< (ostream& os, const bboxset<T,D>& s) {
  s.output(os);
  return os;
}



#endif // BBOXSET_HH
