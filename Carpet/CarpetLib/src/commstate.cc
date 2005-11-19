#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "commstate.hh"


// Communication state control
comm_state::comm_state ()
{
  // If CarpetLib::use_collective_communication_buffers is set to true,
  // this comm_state object will use collective communications,
  // ie. it will step through
  //   state_get_buffer_sizes
  //   state_fill_send_buffers
  //   state_empty_recv_buffers
  //
  // If CarpetLib::use_collective_communication_buffers is false
  // then individual communications on single components are used
  // by stepping through
  //   state_post
  //   state_wait

  DECLARE_CCTK_PARAMETERS;

  thestate = use_collective_communication_buffers ?
             state_get_buffer_sizes : state_post;

  typebufs.resize (dist::c_ndatatypes());
  for (size_t type = 0; type < typebufs.size(); type++) {
    typebufs[type].procbufs.resize(dist::size());
  }

#define INSTANTIATE(T)                                                 \
  {                                                                    \
    T dummy;                                                           \
    int type = dist::c_datatype (dummy);                               \
    typebufs.at(type).datatypesize = sizeof (dummy);                   \
    typebufs.at(type).mpi_datatype = dist::datatype (dummy);           \
  }
#include "instantiate"
#undef INSTANTIATE

  recvbuffers_ready.resize (dist::c_ndatatypes() * dist::size());
  srequests.resize (dist::c_ndatatypes() * dist::size(), MPI_REQUEST_NULL);
  rrequests.resize (dist::c_ndatatypes() * dist::size(), MPI_REQUEST_NULL);
}


void comm_state::step ()
{
  assert (thestate != state_done);
  switch (thestate) {
    case state_post:
      // After all sends/recvs have been posted in 'state_post'
      // now wait for their completion in 'state_wait'.
      thestate = state_wait;
      break;

    case state_wait:
      thestate = state_done;
      break;

    case state_get_buffer_sizes:
      // The sizes of the collective communication buffers are known
      // so now allocate them.
      // The receive operations are also posted here already
      // (a clever MPI layer may take advantage of such early posting).
      num_posted_recvs = num_completed_recvs = 0;

      for (size_t type = 0; type < typebufs.size(); type++) {

        // skip unused datatype buffers
        if (not typebufs[type].in_use) {
          continue;
        }

        int& datatypesize = typebufs[type].datatypesize;
        for (size_t proc = 0; proc < typebufs[type].procbufs.size(); proc++) {
          procbufdesc& procbuf = typebufs[type].procbufs[proc];

          procbuf.sendbufbase = new char[procbuf.sendbufsize*datatypesize];
          procbuf.recvbufbase = new char[procbuf.recvbufsize*datatypesize];
          procbuf.sendbuf = procbuf.sendbufbase;
          procbuf.recvbuf = procbuf.recvbufbase;

          if (procbuf.recvbufsize > 0) {
            MPI_Irecv (procbuf.recvbufbase, procbuf.recvbufsize,
                       typebufs[type].mpi_datatype, proc, type,
                       dist::comm, &rrequests[dist::size()*type + proc]);
            num_posted_recvs++;
          }
        }
      }

      // Now go and get the send buffers filled with data.
      // Once a buffer is full it will be posted right away
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
      if (not AllPostedCommunicationsFinished ()) {
        // No state change if there are still outstanding communications;
        // do another comm_state loop iteration.
      } else {
        // Everything is done so release the collective communication buffers.
        for (size_t type = 0; type < typebufs.size(); type++) {
          for (size_t proc = 0; proc < typebufs[type].procbufs.size(); proc++) {
            delete[] typebufs[type].procbufs[proc].sendbufbase;
            delete[] typebufs[type].procbufs[proc].recvbufbase;
          }
        }
        thestate = state_done;
      }
      break;

    default:
      assert (0 && "invalid state");
  }
}


bool comm_state::done ()
{
  return thestate == state_done;
}


comm_state::~comm_state ()
{
  DECLARE_CCTK_PARAMETERS;

  assert (thestate == state_done or
          thestate == (use_collective_communication_buffers ?
                       state_get_buffer_sizes : state_post));
  assert (requests.empty());
}


// wait for completion of posted collective buffer sends/receives
//
// Depending on the parameter CarpetLib::use_waitall, this function will wait
// for all (true) or at least one (false) of the posted receive operations
// to finish.
//
// It returns true if all posted communications have been completed.
bool comm_state::AllPostedCommunicationsFinished ()
{
  DECLARE_CCTK_PARAMETERS;
 
  // check if all outstanding receives have been completed already
  if (num_posted_recvs == num_completed_recvs) {
    // finalize the outstanding sends in one go
    MPI_Waitall (srequests.size(), &srequests.front(), MPI_STATUSES_IGNORE);

    return true;
  }

  // reset completion flag for all receive buffers
  recvbuffers_ready.assign (recvbuffers_ready.size(), false);

  if (use_waitall) {
    // mark all posted recveive buffers as ready
    for (size_t i = 0; i < recvbuffers_ready.size(); i++) {
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
