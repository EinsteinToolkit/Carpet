// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gf.cc,v 1.14 2003/11/05 16:18:39 schnetter Exp $

#include <assert.h>

#include "cctk.h"

#include "defs.hh"

#include "gf.hh"

using namespace std;



// Constructors
// VGF
template<class T,int D>
gf<T,D>::gf (const int varindex, th<D>& t, dh<D>& d,
	     const int tmin, const int tmax, const int prolongation_order_time)
  : ggf<D>(varindex, t, d, tmin, tmax, prolongation_order_time)
{
  // VGF
  this->recompose (0, true);
}

// Destructors
template<class T,int D>
gf<T,D>::~gf () { }



// Access to the data
template<class T,int D>
const data<T,D>* gf<T,D>::operator() (int tl, int rl, int c, int ml) const {
  assert (tl>=this->tmin && tl<=this->tmax);
  assert (rl>=0 && rl<this->h.reflevels());
  assert (c>=0 && c<this->h.components(rl));
  assert (ml>=0 && ml<this->h.mglevels(rl,c));
  return (const data<T,D>*)this->storage[tl-this->tmin][rl][c][ml];
}

template<class T,int D>
data<T,D>* gf<T,D>::operator() (int tl, int rl, int c, int ml) {
  assert (tl>=this->tmin && tl<=this->tmax);
  assert (rl>=0 && rl<this->h.reflevels());
  assert (c>=0 && c<this->h.components(rl));
  assert (ml>=0 && ml<this->h.mglevels(rl,c));
  return (data<T,D>*)this->storage[tl-this->tmin][rl][c][ml];
}



// Output
template<class T,int D>
ostream& gf<T,D>::output (ostream& os) const {
  T Tdummy;
  os << "gf<" << typestring(Tdummy) << "," << D << ">:"
     << this->varindex << "[" << CCTK_VarName(this->varindex) << "],"
     << "dt=[" << this->tmin << ":" << this->tmax<< "]";
  return os;
}



#define INSTANTIATE(T)				\
template class gf<T,3>;

#include "instantiate"

#undef INSTANTIATE
