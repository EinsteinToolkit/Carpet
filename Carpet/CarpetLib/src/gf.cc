#include <cassert>

#include "cctk.h"

#include "defs.hh"

#include "gf.hh"

using namespace std;



// Constructors
template<typename T>
gf<T>::gf (const int varindex, const operator_type transport_operator,
           th& t, dh& d,
           const int tmin, const int tmax, const int prolongation_order_time,
           const int vectorlength, const int vectorindex,
           gf* const vectorleader)
  : ggf(varindex, transport_operator,
        t, d, tmin, tmax, prolongation_order_time,
        vectorlength, vectorindex, vectorleader)
{
  // this->recompose ();
  this->recompose_crop ();
  for (int rl=0; rl<this->h.reflevels(); ++rl) {
    this->recompose_allocate (rl);
#if 0
    for (comm_state state; !state.done(); state.step()) {
      this->recompose_fill (state, rl, false);
    }
#endif
    this->recompose_free (rl);
#if 0
    for (comm_state state; !state.done(); state.step()) {
      this->recompose_bnd_prolongate (state, rl, false);
    }
    for (comm_state state; !state.done(); state.step()) {
      this->recompose_sync (state, rl, false);
    }
#endif
  } // for rl
}

// Destructors
template<typename T>
gf<T>::~gf () { }



// Access to the data
template<typename T>
const data<T>* gf<T>::operator() (int tl, int rl, int c, int ml) const {
  assert (tl>=this->tmin && tl<=this->tmax);
  assert (rl>=0 && rl<this->h.reflevels());
  assert (c>=0 && c<this->h.components(rl));
  assert (ml>=0 && ml<this->h.mglevels(rl,c));
  return (const data<T>*)this->storage.at(tl-this->tmin).at(rl).at(c).at(ml);
}

template<typename T>
data<T>* gf<T>::operator() (int tl, int rl, int c, int ml) {
  assert (tl>=this->tmin && tl<=this->tmax);
  assert (rl>=0 && rl<this->h.reflevels());
  assert (c>=0 && c<this->h.components(rl));
  assert (ml>=0 && ml<this->h.mglevels(rl,c));
  return (data<T>*)this->storage.at(tl-this->tmin).at(rl).at(c).at(ml);
}



// Output
template<typename T>
ostream& gf<T>::output (ostream& os) const {
  T Tdummy;
  os << "gf<" << typestring(Tdummy) << ">:"
     << this->varindex << "[" << CCTK_VarName(this->varindex) << "],"
     << "dt=[" << this->tmin << ":" << this->tmax<< "]";
  return os;
}



#define INSTANTIATE(T)				\
template class gf<T>;

#include "instantiate"

#undef INSTANTIATE
