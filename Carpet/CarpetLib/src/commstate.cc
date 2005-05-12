#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "commstate.hh"


// Communication state control
comm_state::comm_state (int vartype_) : vartype(vartype_)
{
  // If this comm state is created with a valid (ie. non-negative) CCTK vartype
  // then it is assumed to be used in collective communications for that type
  // later on, ie. it will step through
  //   state_get_buffer_sizes
  //   state_fill_send_buffers
  //   state_empty_recv_buffers
  // This path needs the size of the vartype, the corresponding MPI datatype
  // and the collective communication buffers initialized also.
  //
  // If this comm state is created with a negative vartype then it will be
  // for single-component communications and set up to step through
  //   state_recv
  //   state_send
  //   state_wait

  DECLARE_CCTK_PARAMETERS;

  uses_collective_communication_buffers =
    use_collective_communication_buffers && vartype >= 0;

  if (uses_collective_communication_buffers) {
    vartypesize = CCTK_VarTypeSize (vartype);
    assert (vartypesize > 0);

    thestate = state_get_buffer_sizes;
    switch (vartype) {
#define TYPECASE(N,T)                                                     \
      case  N: { T dummy; datatype = dist::datatype(dummy); } break;
#include "carpet_typecase.hh"
#undef TYPECASE
      default: assert (0);
    }
    collbufs.resize (dist::size());
    rrequests.resize (dist::size(), MPI_REQUEST_NULL);
    srequests.resize (dist::size(), MPI_REQUEST_NULL);
    recvbuffers_ready.resize (dist::size());

  } else {
    thestate = state_recv;
  }
}

void comm_state::step ()
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (thestate!=state_done);
  if (! uses_collective_communication_buffers && combine_recv_send) {
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

    case state_get_buffer_sizes:
      // The sizes of the collective communication buffers are known
      // so now allocate them.
      // The receive operations are also posted here already
      // (a clever MPI layer may take advantage of such early posting).
      num_posted_recvs = num_completed_recvs = 0;

      for (size_t i = 0; i < collbufs.size(); i++) {
        collbufs[i].sendbufbase = new char[collbufs[i].sendbufsize*vartypesize];
        collbufs[i].recvbufbase = new char[collbufs[i].recvbufsize*vartypesize];
        collbufs[i].sendbuf = collbufs[i].sendbufbase;
        collbufs[i].recvbuf = collbufs[i].recvbufbase;

        if (collbufs[i].recvbufsize > 0) {
          MPI_Irecv (collbufs[i].recvbufbase, collbufs[i].recvbufsize,
                     datatype, i, MPI_ANY_TAG, dist::comm, &rrequests[i]);
          num_posted_recvs++;
        }
      }

      // Now go and get the send buffers filled with data.
      // Once a buffer is full it will get posted right away
      // (see gdata::copy_into_sendbuffer() and
      //  gdata::interpolate_into_sendbuffer()).
      thestate = state_fill_send_buffers;
      break;

    case state_fill_send_buffers:
      // Now fall through to the next state in which the recv buffers
      // are emptied as soon as data has arrived.
      thestate = state_empty_recv_buffers;

    case state_empty_recv_buffers:
      // Finish (at least one of) the posted communications
      if (! AllPostedCommunicationsFinished (use_waitall)) {
        // No state change if there are still outstanding communications;
        // do another comm_state loop iteration.
      } else {
        // Everything is done so release the collective communication buffers.
        for (size_t i = 0; i < collbufs.size(); i++) {
          delete[] collbufs[i].sendbufbase;
          delete[] collbufs[i].recvbufbase;
        }
        thestate = state_done;
      }
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
  assert (thestate == state_done ||
          thestate == uses_collective_communication_buffers ?
                      state_get_buffer_sizes : state_recv);
  assert (tmps1.empty());
  assert (tmps2.empty());
  assert (requests.empty());
}


// wait for completion of posted collective buffer sends/receives
//
// Depending on use_waitall, this function will wait for all at once (true)
// or at least one (false) of the posted receive operations to finish.
//
// It returns true if all posted communications have been completed.
bool comm_state::AllPostedCommunicationsFinished (bool use_waitall)
{
  // check if all outstanding receives have been completed already
  if (num_posted_recvs == num_completed_recvs) {
    // finalize the outstanding sends in one go
    MPI_Waitall (srequests.size(), &srequests.front(), MPI_STATUSES_IGNORE);

    return true;
  }

  // reset completion flag for all receive buffers
  for (size_t i = 0; i < recvbuffers_ready.size(); i++) {
    recvbuffers_ready[i] = false;
  }

  if (use_waitall) {
    // mark all posted recveive buffers as ready
    for (size_t i = 0; i < rrequests.size(); i++) {
      recvbuffers_ready[i] = rrequests[i] != MPI_REQUEST_NULL;
    }

    // wait for completion of all posted receive operations
    MPI_Waitall (rrequests.size(), &rrequests.front(), MPI_STATUSES_IGNORE);
    num_completed_recvs = num_posted_recvs;
  } else {
    int num_completed_recvs_ = 0;
    vector<int> completed_recvs(rrequests.size(), -1);

    // wait for completion of at least one posted receive operation
    MPI_Waitsome (rrequests.size(), &rrequests.front(), &num_completed_recvs_,
                  &completed_recvs.front(), MPI_STATUSES_IGNORE);
    assert (0 < num_completed_recvs_);
    num_completed_recvs += num_completed_recvs_;

    // mark the recveive buffers of completed communications as ready
    for (int i = 0; i < num_completed_recvs_; i++) {
      assert (rrequests.at(completed_recvs.at(i)) == MPI_REQUEST_NULL);
      recvbuffers_ready.at(completed_recvs.at(i)) = true;
    }
  }

  return (false);
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
