// $Header:$

#ifndef BBOX_HH
#define BBOX_HH

#include <iostream>

#include "defs.hh"
#include "vect.hh"

using namespace std;



// Forward declaration
template<class T, int D> class bbox;

// Input/Output
template<class T, int D>
istream& operator>> (istream& is, bbox<T,D>& b);
template<class T, int D>
ostream& operator<< (ostream& os, const bbox<T,D>& b);



/**
 * A bounding box, i.e. a rectangle with lower and upper bound and a
 * stride.
 */
template<class T, int D>
class bbox {
  
  // Fields
  
  /** Bounding box bounds and stride.  The bounds are inclusive.  */
  vect<T,D> _lower, _upper, _stride;
  
public:
  
  // Constructors
  
  /** Construct an empty bbox.  */
  bbox ();
  
  /** Copy constructor.  */
  bbox (const bbox& b);
  
  /** Assignment operator.  */
  bbox& operator= (const bbox& b);
  
  /** Create a bbox from bounds and stride.  */
  bbox (const vect<T,D>& lower, const vect<T,D>& upper,
	const vect<T,D>& stride);
  
  // Accessors
  // (Don't return references; *this might be a temporary)
  
  /** Get lower bound.  */
  vect<T,D> lower () const { return _lower; }
  
  /** Get upper bound.  */
  vect<T,D> upper () const { return _upper; }
  
  /** Get stride.  */
  vect<T,D> stride () const { return _stride; }
  
  /** Get the shape (or extent).  */
  vect<T,D> shape () const { return _upper - _lower + _stride; }
  
  /** Determine whether the bbox is empty.  */
  bool empty() const {
    return any(lower()>upper());
  }
  
  /** Return the size, which is the product of the shape.  */
  T size () const;
  
  // Queries
  
  /** Find out whether the bbox contains the point x.  */
  bool contains (const vect<T,D>& x) const;
  
  // Operators
  bool operator== (const bbox& b) const;
  bool operator!= (const bbox& b) const;
  bool operator< (const bbox& b) const;
  bool operator> (const bbox& b) const;
  bool operator<= (const bbox& b) const;
  bool operator>= (const bbox& b) const;
  
  /** Calculate the intersection (the set of common points) with the
      bbox b.  */
  bbox operator& (const bbox& b) const;
  
  /** Find out whether this bbox is contained in the bbox b.  */
  bool is_contained_in (const bbox& b) const;
  
  /** Find out whether this bbox is aligned with the bbox b.
      ("aligned" means that both bboxes have the same stride and that
      their boundaries are commesurate.)  */
  bool is_aligned_with (const bbox& b) const;
  
  /** Expand (enlarge) the bbox by multiples of the stride.  */
  bbox expand (const vect<T,D>& lo, const vect<T,D>& hi) const;
  
  /** Find the smallest b-compatible box around this bbox.
      ("compatible" means having the same stride.)  */
  bbox expanded_for (const bbox& b) const;
  
  /** Find the largest b-compatible box inside this bbox.  */
  bbox contracted_for (const bbox& b) const;
  
  /** Find the smallest bbox containing both boxes.  */
  bbox expanded_containing (const bbox<T,D>& b) const;
  
  // Iterators
  
  /** An iterator over all points in a bbox.  */
  class iterator {
  protected:
    /** The bbox over which we iterate.  */
    const bbox& box;
    /** Current position.  */
    vect<T,D> pos;
  public:
    /** Constructor.  */
    iterator (const bbox& box, const vect<T,D>& pos);
    /** Accessor.  */
    const vect<T,D>& operator* () const { return pos; }
    /** Check whether the position is the same.  */
    bool operator!= (const iterator& i) const;
    /** Advance.  */
    iterator& operator++ ();
  };
  
  /** Create an iterator that points to the first point in a bbox.  */
  iterator begin () const;
  /** Create an iterator that points "after the last point" in a bbox,
      which means that it also points to the first point.  */
  iterator end () const;
  
  // Input/Output helpers
  void input (istream& is);
  void output (ostream& os) const;
};



// Input

/** Read a formatted bbox from a stream.  */
template<class T,int D>
inline istream& operator>> (istream& is, bbox<T,D>& b) {
  b.input(is);
  return is;
}



// Output

/** Write a bbox formatted to a stream.  */
template<class T,int D>
inline ostream& operator<< (ostream& os, const bbox<T,D>& b) {
  b.output(os);
  return os;
}



#endif // BBOX_HH
