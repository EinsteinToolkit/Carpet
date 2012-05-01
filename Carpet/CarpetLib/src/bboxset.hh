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
template<typename T, int D> class bboxset;

// template<typename T,int D>
// bboxset<T,D> operator+ (const bbox<T,D>& b1, const bbox<T,D>& b2);
// template<typename T,int D>
// bboxset<T,D> operator+ (const bbox<T,D>& b, const bboxset<T,D>& s);

// template<typename T,int D>
// bboxset<T,D> operator- (const bbox<T,D>& b1, const bbox<T,D>& b2);
// template<typename T,int D>
// bboxset<T,D> operator- (const bbox<T,D>& b, const bboxset<T,D>& s);

// Input
template<typename T,int D>
istream& operator>> (istream& is, bboxset<T,D>& s);

// Output
template<typename T,int D>
ostream& operator<< (ostream& os, const bboxset<T,D>& s);



// Bounding box set class
template<typename T, int D>
class bboxset {
  
  // Cost annotations depend on the number of bset elements n, and
  // assume that normalization is skipped.
  
  struct skip_normalize_t {
    bboxset<T,D>& s;
    bool const saved_skip_normalize;
    skip_normalize_t(bboxset<T,D>& s_)
      : s(s_),
        saved_skip_normalize(s.skip_normalize)
    {
      s.skip_normalize = true;
    }
    ~skip_normalize_t()
    {
      s.skip_normalize = saved_skip_normalize;
      s.normalize();
    }
  };
  
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
  
  bool skip_normalize;
  
public:
  
  // Constructors
  bboxset ();                   // cost: O(1)
  bboxset (const box& b);       // cost: O(1)
  bboxset (const bboxset& s);   // cost: O(n)
  
  bboxset (const list<box>& lb);
  bboxset (const vector<box>& lb);
  bboxset (const vector<list<box> >& vlb);
  template<typename U>
  bboxset (const vector<U>& vb, const bbox<T,D> U::* const v);
  template<typename U>
  bboxset (const vector<U>& vb, const bboxset U::* const v);
  
  static bboxset poison ();
  
  // Invariant
  bool invariant () const CCTK_MEMBER_ATTRIBUTE_PURE;
  
private:
  // Normalisation
  void normalize ();
  
public:
  // Accessors
  bool empty () const { return bs.empty(); } // cost: O(1)
  // T size () const;
  typedef typename box::size_type size_type;
  size_type size () const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  int setsize () const { return bs.size(); } // cost: O(1)
  
  // Find out whether this bboxset intersects the bbox b
  bool intersects (const box& b) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  
  // Add (bboxes that don't overlap)
  bboxset& operator+= (const box& b);     // cost: O(1)
  bboxset& operator+= (const bboxset& s); // cost: O(n)
  bboxset& add_transfer (bboxset& s);     // cost: O(1)
  bboxset operator+ (const box& b)
    const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  bboxset operator+ (const bboxset& s)
    const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  
  // Union
  bboxset& operator|= (const box& b);     // cost: O(n)
  bboxset& operator|= (const bboxset& s); // cost: O(n^2)
  bboxset operator| (const box& b)
    const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  bboxset operator| (const bboxset& s)
    const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n^2)
  
  // Intersection
  bboxset operator& (const box& b)
    const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  bboxset operator& (const bboxset& s)
    const CCTK_MEMBER_ATTRIBUTE_PURE;     // cost: O(n)
  bboxset& operator&= (const box& b);     // cost: O(n)
  bboxset& operator&= (const bboxset& s); // cost: O(n)
  
  // Difference
  // friend bboxset operator- <T,D>(const box& b1, const box& b2);
  static bboxset minus (const box& b1, const box& b2)
    CCTK_ATTRIBUTE_PURE;        // cost: O(1)
  bboxset operator- (const box& b)
    const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  // cost: O(n)
  bboxset& operator-= (const box& b)
  {
    *this = *this - b;
    assert (invariant());
    return *this;
  }
  bboxset& operator-= (const bboxset& s);     // cost: O(n^2)
  bboxset operator- (const bboxset& s)
    const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n^2)
  // friend bboxset operator- <T,D>(const box& b, const bboxset& s);
  static bboxset minus (const box& b, const bboxset& s)
    CCTK_ATTRIBUTE_PURE;        // cost: O(n^2)
  
  /** Find a bbox containing the whole set.  */
  box container () const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  /** Find the pseudo-inverse.  */
  bboxset pseudo_inverse (const int n) const CCTK_MEMBER_ATTRIBUTE_PURE;
  
  /** Expand (enlarge) the bboxset by multiples of the stride.  */
  bboxset expand (const vect<T,D>& lo, const vect<T,D>& hi)
    const CCTK_MEMBER_ATTRIBUTE_PURE;
  bboxset expand (const vect<vect<T,D>,2>& lohi) const
  { return expand (lohi[0], lohi[1]); }
  
  /** Shift the bboxset by multiples of the stride.  */
  bboxset shift (const vect<T,D>& v) const
  { return expand (-v, v); }
  
  /** Expand (enlarge) the bboxset by multiples of a fraction of the
      stride.  */
  // cost: O(n^2) in general, but only O(n) for shifting
  bboxset expand (const vect<T,D>& lo, const vect<T,D>& hi,
                  const vect<T,D>& denom) const CCTK_MEMBER_ATTRIBUTE_PURE;
  // cost: O(n^2) in general, but only O(n) for shifting
  bboxset expand (const vect<vect<T,D>,2>& lohi, const vect<T,D>& denom) const
  { return expand (lohi[0], lohi[1], denom); }
  
  /** Shift the bboxset by multiples of a fraction of the stride.  */
// cost: O(n)
  bboxset shift (const vect<T,D>& v, const vect<T,D>& denom) const
  { return expand (-v, v, denom); }
  
  /** Find the smallest b-compatible box around this bbox.
      ("compatible" means having the same stride.)  */
  bboxset expanded_for (const box& b) const CCTK_MEMBER_ATTRIBUTE_PURE;
  
#warning "TODO: this is incorrect"
#if 1
  /** Find the largest b-compatible box inside this bbox.  */
  bboxset contracted_for (const box& b) const CCTK_MEMBER_ATTRIBUTE_PURE;
#endif
  
  // Equality
  bool operator== (const bboxset& s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator!= (const bboxset& s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator< (const bboxset& s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator<= (const bboxset& s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator> (const bboxset& s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator>= (const bboxset& s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  
  // Iterators
  typedef typename bset::const_iterator const_iterator;
  typedef typename bset::iterator       iterator;
  
  const_iterator begin () const { return bs.begin(); }
  const_iterator end () const   { return bs.end(); }
//   iterator begin () const { return bs.begin(); }
//   iterator end () const   { return bs.end(); }
  
  // Memory usage
  size_t memory () const { return memoryof (bs); }
  
  // Input
  istream& input (istream& is);
  
  // Output
  ostream& output (ostream& os) const;
};



template<typename T,int D>
inline bboxset<T,D> operator+ (const bbox<T,D>& b1, const bbox<T,D>& b2) {
  return bboxset<T,D>(b1) + bboxset<T,D>(b2);
}

template<typename T,int D>
inline bboxset<T,D> operator+ (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>(b) + s;
}

// cost: O(1)
template<typename T,int D>
inline bboxset<T,D> operator- (const bbox<T,D>& b1, const bbox<T,D>& b2) {
  return bboxset<T,D>::minus(b1,b2);
}

// cost: O(n^2)
template<typename T,int D>
inline bboxset<T,D> operator- (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return bboxset<T,D>::minus(b,s);
}



template<typename T,int D>
inline bboxset<T,D> operator| (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return s | b;
}

template<typename T,int D>
inline bboxset<T,D> operator& (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return s & b;
}



template<typename T,int D>
inline bool operator== (const bbox<T,D>& b, const bboxset<T,D>& s)
{
  return bboxset<T,D>(b) == s;
}

template<typename T,int D>
inline bool operator!= (const bbox<T,D>& b, const bboxset<T,D>& s)
{
  return bboxset<T,D>(b) != s;
}

template<typename T,int D>
inline bool operator< (const bbox<T,D>& b, const bboxset<T,D>& s)
{
  return bboxset<T,D>(b) < s;
}

template<typename T,int D>
inline bool operator<= (const bbox<T,D>& b, const bboxset<T,D>& s)
{
  return bboxset<T,D>(b) <= s;
}

template<typename T,int D>
inline bool operator> (const bbox<T,D>& b, const bboxset<T,D>& s)
{
  return bboxset<T,D>(b) > s;
}

template<typename T,int D>
inline bool operator>= (const bbox<T,D>& b, const bboxset<T,D>& s)
{
  return bboxset<T,D>(b) >= s;
}



template<typename T,int D>
inline bool operator== (const bboxset<T,D>& s, const bbox<T,D>& b)
{
  return s == bboxset<T,D>(b);
}

template<typename T,int D>
inline bool operator!= (const bboxset<T,D>& s, const bbox<T,D>& b)
{
  return s != bboxset<T,D>(b);
}

template<typename T,int D>
inline bool operator< (const bboxset<T,D>& s, const bbox<T,D>& b)
{
  return s < bboxset<T,D>(b);
}

template<typename T,int D>
inline bool operator<= (const bboxset<T,D>& s, const bbox<T,D>& b)
{
  return s <= bboxset<T,D>(b);
}

template<typename T,int D>
inline bool operator> (const bboxset<T,D>& s, const bbox<T,D>& b)
{
  return s > bboxset<T,D>(b);
}

template<typename T,int D>
inline bool operator>= (const bboxset<T,D>& s, const bbox<T,D>& b)
{
  return s >= bboxset<T,D>(b);
}



// Memory usage
template<typename T, int D>
inline size_t memoryof (bboxset<T,D> const & s)
{ return s.memory(); }



// Input
template<typename T,int D>
inline istream& operator>> (istream& is, bboxset<T,D>& s) {
  return s.input(is);
}



// Output
template<typename T,int D>
inline ostream& operator<< (ostream& os, const bboxset<T,D>& s) {
  return s.output(os);
}



#endif // BBOXSET_HH
