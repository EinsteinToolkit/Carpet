#include "cctk.h"
#include "cctk_Parameters.h"

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
