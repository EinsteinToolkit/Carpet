#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "commstate.hh"



// Communication state control
comm_state::comm_state ()
  : thestate(state_recv)
{
}

void comm_state::step ()
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (thestate!=state_done);
  if (combine_recv_send) {
    switch (thestate) {
    case state_recv:
      assert (tmps1.empty());
      thestate = state_wait;
      break;
    case state_send:
      assert (0);
    case state_wait:
      assert (tmps1.empty());
      assert (tmps2.empty());
      thestate = state_done;
      break;
    case state_done:
      assert (0);
    default:
      assert (0);
    }
  } else {
    switch (thestate) {
    case state_recv:
      assert (tmps2.empty());
      thestate = state_send;
      break;
    case state_send:
      assert (tmps1.empty());
      thestate = state_wait;
      break;
    case state_wait:
      assert (tmps1.empty());
      assert (tmps2.empty());
      thestate = state_done;
      break;
    case state_done:
      assert (0);
    default:
      assert (0);
    }
  }
}

bool comm_state::done ()
{
  return thestate==state_done;
}

comm_state::~comm_state ()
{
  assert (thestate==state_recv || thestate==state_done);
  assert (tmps1.empty());
  assert (tmps2.empty());
  assert (requests.empty());
}



comm_state::gcommbuf::gcommbuf ()
{
}

comm_state::gcommbuf::~gcommbuf ()
{
}




template<typename T>
comm_state::commbuf<T>::commbuf (ibbox const & box)
{
  data.resize (prod (box.shape() / box.stride()));
}

template<typename T>
comm_state::commbuf<T>::~commbuf ()
{
}

template<typename T>
void const *
comm_state::commbuf<T>::pointer ()
  const
{
  return &data.front();
}

template<typename T>
void *
comm_state::commbuf<T>::pointer ()
{
  return &data.front();
}

template<typename T>
int
comm_state::commbuf<T>::size ()
  const
{
  return data.size();
}

template<typename T>
MPI_Datatype
comm_state::commbuf<T>::datatype ()
  const
{
  T dummy;
  return dist::datatype (dummy);
}



#define INSTANTIATE(T)                          \
  template class comm_state::commbuf<T>;

#include "instantiate"

#undef INSTANTIATE
