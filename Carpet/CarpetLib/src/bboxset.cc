#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <set>
#include <stack>

#include "defs.hh"

#include "bboxset.hh"

using namespace std;



// Constructors
template<class T, int D>
bboxset<T,D>::bboxset () {
  assert (invariant());
}

template<class T, int D>
bboxset<T,D>::bboxset (const box& b) {
  //S if (!b.empty()) bs.insert(b);
  if (!b.empty()) bs.push_back(b);
  assert (invariant());
}

template<class T, int D>
bboxset<T,D>::bboxset (const bboxset& s): bs(s.bs) {
  assert (invariant());
}

template<class T, int D>
bboxset<T,D>::bboxset (const list<box>& lb) {
  for (typename list<box>::const_iterator
         li = lb.begin(); li != lb.end(); ++ li)
  {
    *this |= *li;
  }
  normalize();
}

template<class T, int D>
bboxset<T,D>::bboxset (const vector<list<box> >& vlb) {
  for (typename vector<list<box> >::const_iterator
         vli = vlb.begin(); vli != vlb.end(); ++ vli)
  {
    *this |= bboxset (*vli);
  }
  normalize();
}



// Invariant
template<class T, int D>
bool bboxset<T,D>::invariant () const {
// This is very slow when there are many bboxes
#if 0 && defined(CARPET_DEBUG)
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    if ((*bi).empty()) return false;
    if (not (*bi).is_aligned_with(*bs.begin())) return false;
    // check for overlap (quadratic -- expensive)
    for (const_iterator bi2=begin(); bi2!=bi; ++bi2) {
      if ((*bi2).intersects(*bi)) return false;
    }
  }
#endif
  return true;
}



// Normalisation
template<class T, int D>
void bboxset<T,D>::normalize ()
{
  assert (invariant());
  
  bboxset const oldbs = * this;
  size_type const oldsize = this->size();
  
  // Split all bboxes into small pieces which have all their
  // boundaries aligned.
  for (int d=0; d<D; ++d) {
    // Find all boundaries
    //S typedef set<T> buf;
    typedef vector<T> buf;
    buf sbnds;
    sbnds.reserve (2 * bs.size());
    for (typename bset::const_iterator si = bs.begin(); si != bs.end(); ++ si) {
      box const & b = * si;
      int const bstr = b.stride()[d];
      int const blo = b.lower()[d];
      int const bhi = b.upper()[d] + bstr;
      //S sbnds.insert (blo);
      //S sbnds.insert (bhi);
      sbnds.push_back (blo);
      sbnds.push_back (bhi);
    }
    sort (sbnds.begin(), sbnds.end());
    typename buf::iterator const last = unique (sbnds.begin(), sbnds.end());
    sbnds.resize (last - sbnds.begin());
    // Split bboxes
    bset nbs;
    for (typename bset::const_iterator si = bs.begin(); si != bs.end(); ++ si) {
      box const & b = * si;
      int const bstr = b.stride()[d];
      int const blo = b.lower()[d];
      int const bhi = b.upper()[d] + bstr;
      typename buf::const_iterator const ilo
        = lower_bound (sbnds.begin(), sbnds.end(), blo);
      typename buf::const_iterator const ihi
        = lower_bound (sbnds.begin(), sbnds.end(), bhi);
      assert (ilo != sbnds.end());
      assert (ihi != sbnds.end());
      assert (* ilo == blo);
      assert (* ihi == bhi);
      // Split one bbox
      for (typename buf::const_iterator curr = ilo; curr != ihi; ++ curr) {
        typename buf::const_iterator next = curr;
        advance (next, 1);
        int const nblo = * curr;
        int const nbhi = * next;
        assert (nbhi > nblo);   // ensure that the set remains sorted
        box nb (b.lower().replace(d, nblo),
                b.upper().replace(d, nbhi - bstr),
                b.stride());
        //S nbs.insert (nb);
        nbs.push_back (nb);
      }
    }
    // Replace old set
    bs.swap (nbs);
    assert (invariant());
  }
  
  // Combine bboxes if possible
  for (int d=0; d<D; ++d) {
    bset nbs;
    while (not bs.empty()) {
      typename bset::iterator si = bs.begin();
      assert (si != bs.end());
      
      box const b = * si;
      int const bstr = b.stride()[d];
      int const blo = b.lower()[d];
      int const bhi = b.upper()[d] + bstr;
      
      for (typename bset::iterator nsi = nbs.begin(); nsi != nbs.end(); ++ nsi)
      {
        box const nb = * nsi;
        int const nblo = nb.lower()[d];
        int const nbhi = nb.upper()[d] + bstr;
        
        box const mb (nb.lower().replace(d, blo),
                      nb.upper().replace(d, bhi - bstr),
                      nb.stride());
        
        // Check whether the other dimensions match
        if (b == mb) {
          // Check whether the bboxes are adjacent in this dimension
          if (nbhi == blo) {
            // Combine boxes, nb < b
            box const cb (b.lower().replace(d, nblo),
                          b.upper(),
                          b.stride());
            bs.erase (si);
            nbs.erase (nsi);
            //S bs.insert (cb);
            bs.push_back (cb);
            goto done;
          } else if (bhi == nblo) {
            // Combine boxes, b < nb
            box const cb (b.lower(),
                          b.upper().replace(d, nbhi - bstr),
                          b.stride());
            bs.erase (si);
            nbs.erase (nsi);
            //S bs.insert (cb);
            bs.push_back (cb);
            goto done;
          }
        }
      }
      bs.erase (si);
      //S nbs.insert (b);
      nbs.push_back (b);
    done:;
    }
    bs.swap (nbs);
    assert (invariant());
  }
  
  size_type const newsize = this->size();
  
  assert (*this == oldbs);
  assert (newsize <= oldsize);
}



// Accessors
template<class T, int D>
typename bboxset<T,D>::size_type bboxset<T,D>::size () const {
  size_type s=0;
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    const size_type bsz = (*bi).size();
    assert (numeric_limits<size_type>::max() - bsz >= s);
    s += bsz;
  }
  return s;
}




// Queries

// Intersection
template<class T, int D>
bool bboxset<T,D>::intersects (const box& b) const {
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    if ((*bi).intersects(b)) return true;
  }
  return false;
}



// Add (bboxes that don't overlap)
template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator+= (const bboxset& s) {
  for (const_iterator bi=s.begin(); bi!=s.end(); ++bi) {
    *this += *bi;
  }
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::add_transfer (bboxset& s) {
  bs.splice (bs.end(), s.bs);
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator+ (const box& b) const {
  bboxset r(*this);
  r += b;
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator+ (const bboxset& s) const {
  bboxset r(*this);
  r += s;
  assert (r.invariant());
  return r;
}



// Union
template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator|= (const box& b) {
  bboxset tmp = b - *this;
  add_transfer (tmp);
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator|= (const bboxset& s) {
  bboxset tmp = s - *this;
  add_transfer (tmp);
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator| (const box& b) const {
  bboxset r(*this);
  r |= b;
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator| (const bboxset& s) const {
  bboxset r(*this);
  r |= s;
  assert (r.invariant());
  return r;
}



// Intersection
template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator& (const box& b) const {
  // start with an empty set
  bboxset r;
  // walk all my elements
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    // insert the intersection with the bbox
    r += *bi & b;
  }
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator& (const bboxset& s) const {
  // start with an empty set
  bboxset r;
  // walk all the bboxes
  for (const_iterator bi=s.begin(); bi!=s.end(); ++bi) {
    // insert the intersection with this bbox
    bboxset tmp = *this & *bi;
    r.add_transfer (tmp);
  }
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator&= (const box& b) {
  *this = *this & b;
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator&= (const bboxset& s) {
  *this = *this & s;
  assert (invariant());
  return *this;
}



// Difference
template<class T, int D>
bboxset<T,D> bboxset<T,D>::minus (const bbox<T,D>& b1, const bbox<T,D>& b2) {
  assert (b1.is_aligned_with(b2));
  if (b1.empty()) return bboxset<T,D>();
  if (not b1.intersects(b2)) return bboxset<T,D>(b1);
  vect<T,D> const b2lo = max (b1.lower(), b2.lower());
  vect<T,D> const b2up = min (b1.upper(), b2.upper());
  vect<T,D> const & b1lo = b1.lower();
  vect<T,D> const & b1up = b1.upper();
  vect<T,D> const & str = b1.stride();
  bboxset<T,D> r;
  for (int d=0; d<D; ++d) {
    // make resulting bboxes as large as possible in x-direction (for
    // better consumption by Fortranly ordered arrays)
    vect<T,D> lb, ub;
    bbox<T,D> b;
    for (int dd=0; dd<D; ++dd) {
      if (dd<d) {
	lb[dd] = b2lo[dd];
	ub[dd] = b2up[dd];
      } else if (dd>d) {
	lb[dd] = b1lo[dd];
	ub[dd] = b1up[dd];
      }
    }
    lb[d] = b1lo[d];
    ub[d] = b2lo[d] - str[d];
    b = bbox<T,D>(lb,ub,str);
    r += b;
    lb[d] = b2up[d] + str[d];
    ub[d] = b1up[d];
    b = bbox<T,D>(lb,ub,str);
    r += b;
  }
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator- (const box& b) const {
  // start with an empty set
  bboxset r;
  // walk all my elements
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    // insert the difference with the bbox
    bboxset tmp = *bi - b;
    r.add_transfer (tmp);
  }
  assert (r.invariant());
  return r;
}
  
template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator-= (const bboxset& s) {
  for (const_iterator bi=s.begin(); bi!=s.end(); ++bi) {
    *this -= *bi;
  }
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator- (const bboxset& s) const {
  bboxset r(*this);
  r -= s;
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::minus (const bbox<T,D>& b, const bboxset<T,D>& s) {
  bboxset<T,D> r = bboxset<T,D>(b) - s;
  assert (r.invariant());
  return r;
}



// Equality
template<class T, int D>
bool bboxset<T,D>::operator<= (const bboxset<T,D>& s) const {
  return (*this - s).empty();
}

template<class T, int D>
bool bboxset<T,D>::operator< (const bboxset<T,D>& s) const {
  return (*this - s).empty() && ! (s - *this).empty();
}

template<class T, int D>
bool bboxset<T,D>::operator>= (const bboxset<T,D>& s) const {
  return s <= *this;
}

template<class T, int D>
bool bboxset<T,D>::operator> (const bboxset<T,D>& s) const {
  return s < *this;
}

template<class T, int D>
bool bboxset<T,D>::operator== (const bboxset<T,D>& s) const {
  return (*this <= s) && (*this >= s);
}

template<class T, int D>
bool bboxset<T,D>::operator!= (const bboxset<T,D>& s) const {
  return ! (*this == s);
}



// Output
template<class T,int D>
void bboxset<T,D>::output (ostream& os) const {
  T Tdummy;
  os << "bboxset<" << typestring(Tdummy) << "," << D << ">:"
     << "size=" << size() << ","
     << "setsize=" << setsize() << ","
     << "set=" << bs;
}



template class bboxset<int,3>;
