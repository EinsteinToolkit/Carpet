#include <cassert>

#include "cctk.h"

#include "defs.hh"

#include "gf.hh"

using namespace std;



// Constructors
template<typename T>
gf<T>::gf (const int varindex_, const operator_type transport_operator_,
           th& t_, dh& d_,
           const int tmin_, const int tmax_,
           const int prolongation_order_time_,
           const int vectorlength_, const int vectorindex_,
           gf* const vectorleader_)
  : ggf(varindex_, transport_operator_,
        t_, d_, tmin_, tmax_, prolongation_order_time_,
        vectorlength_, vectorindex_, vectorleader_)
{
  // recompose ();
  recompose_crop ();
  for (int rl=0; rl<h.reflevels(); ++rl) {
    recompose_allocate (rl);
#if 0
    for (comm_state state; !state.done(); state.step()) {
      recompose_fill (state, rl, false);
    }
#endif
    recompose_free (rl);
#if 0
    for (comm_state state; !state.done(); state.step()) {
      recompose_bnd_prolongate (state, rl, false);
    }
    for (comm_state state; !state.done(); state.step()) {
      recompose_sync (state, rl, false);
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
  assert (tl>=tmin && tl<=tmax);
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  return (const data<T>*)storage.at(tl-tmin).at(rl).at(c).at(ml);
}

template<typename T>
data<T>* gf<T>::operator() (int tl, int rl, int c, int ml) {
  assert (tl>=tmin && tl<=tmax);
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  return (data<T>*)storage.at(tl-tmin).at(rl).at(c).at(ml);
}



// Output
template<typename T>
ostream& gf<T>::output (ostream& os) const {
  T Tdummy;
  os << "gf<" << typestring(Tdummy) << ">:"
     << varindex << "[" << CCTK_VarName(varindex) << "],"
     << "dt=[" << tmin << ":" << tmax<< "]";
  return os;
}



#define INSTANTIATE(T)				\
template class gf<T>;

#include "instantiate"

#undef INSTANTIATE
