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
      // so now allocate the send buffers.
      for (int i = 0; i < collbufs.size(); i++) {
//fprintf (stderr, "proc %d: allocating buf[%d] = %d elems\n", dist::rank(), i, collbufs[i].sendbufsize);
        collbufs[i].sendbufbase = new char[collbufs[i].sendbufsize*vartypesize];
        collbufs[i].sendbuf = collbufs[i].sendbufbase;
      }

      // Now go and get the send buffers filled with data.
      thestate = state_fill_send_buffers;
      break;

    case state_fill_send_buffers:
      // The send buffers contain all the data so now allocate the recv buffers.
      for (int i = 0; i < collbufs.size(); i++) {
        collbufs[i].recvbufbase = new char[collbufs[i].recvbufsize*vartypesize];
        collbufs[i].recvbuf = collbufs[i].recvbufbase;
      }

      // Exchange pairs of send/recv buffers between all processors.
      ExchangeBuffers ();

      // The send buffers are not needed anymore.
      for (int i = 0; i < collbufs.size(); i++) {
        delete[] collbufs[i].sendbufbase;
      }

      // All data has been received so now go and empty the recv buffers.
      thestate = state_empty_recv_buffers;
      break;
    case state_empty_recv_buffers:
      // Finally release the recv buffers.
      for (int i = 0; i < collbufs.size(); i++) {
        delete[] collbufs[i].recvbufbase;
      }

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
  assert (thestate == state_done ||
          thestate == uses_collective_communication_buffers ?
                      state_get_buffer_sizes : state_recv);
  assert (tmps1.empty());
  assert (tmps2.empty());
  assert (requests.empty());
}


// exchange the collective communication buffers between processors
void comm_state::ExchangeBuffers ()
{
  vector<MPI_Request> rrequests(collbufs.size(), MPI_REQUEST_NULL);
  vector<MPI_Request> srequests(collbufs.size(), MPI_REQUEST_NULL);

  // loop over all buffers (alias processors) and initiate the communications
  for (int proc = 0; proc < collbufs.size(); proc++) {

    // don't send zero-sized messages (this covers 'proc == rank')
    if (collbufs[proc].recvbufsize > 0) {
      MPI_Irecv (collbufs[proc].recvbufbase, collbufs[proc].recvbufsize,
                 datatype, proc, MPI_ANY_TAG, dist::comm, &rrequests[proc]);
    }
    if (collbufs[proc].sendbufsize > 0) {
      MPI_Isend (collbufs[proc].sendbufbase, collbufs[proc].sendbufsize,
                 datatype, proc, 0,           dist::comm, &srequests[proc]);
    }
  }

  // finalize the outstanding communications
  // FIXME: use MPI_Waitsome() on rrequests
  MPI_Waitall (rrequests.size(), &rrequests.front(), MPI_STATUSES_IGNORE);
  MPI_Waitall (srequests.size(), &srequests.front(), MPI_STATUSES_IGNORE);
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
