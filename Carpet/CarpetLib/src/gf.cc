#include <cassert>

#include "cctk.h"

#include "defs.hh"

#include "gf.hh"

using namespace std;



// Constructors
template<typename T>
gf<T>::gf (const int varindex_, const operator_type transport_operator_,
           th& t_, dh& d_,
           const int prolongation_order_time_,
           const int vectorlength_, const int vectorindex_,
           gf* const vectorleader_)
  : ggf(varindex_, transport_operator_,
        t_, d_, prolongation_order_time_,
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
gf<T>::~gf ()
{
}



// Access to the data
template<typename T>
const data<T>* gf<T>::operator() (int tl, int rl, int c, int ml) const
{
  assert (rl>=0 and rl<h.reflevels());
  assert (c>=0 and c<h.components(rl));
  assert (ml>=0 and ml<h.mglevels());
  assert (tl>=0 and tl<timelevels(ml, rl));
  return (const data<T>*)storage.AT(ml).AT(rl).AT(c).AT(tl);
}

template<typename T>
data<T>* gf<T>::operator() (int tl, int rl, int c, int ml)
{
  assert (rl>=0 and rl<h.reflevels());
  assert (c>=0 and c<h.components(rl));
  assert (ml>=0 and ml<h.mglevels());
  assert (tl>=0 and tl<timelevels(ml, rl));
  return (data<T>*)storage.AT(ml).AT(rl).AT(c).AT(tl);
}



// Output
template<typename T>
ostream& gf<T>::output (ostream& os) const
{
  T Tdummy;
  os << "gf<" << typestring(Tdummy) << ">:"
     << varindex << "[" << CCTK_VarName(varindex) << "],"
     << "tls=" << timelevels_;
  return os;
}



#define INSTANTIATE(T)				\
template class gf<T>;

#include "instantiate"

#undef INSTANTIATE
