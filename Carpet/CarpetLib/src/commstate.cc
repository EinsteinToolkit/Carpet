#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "commstate.hh"
#include "timestat.hh"



using namespace std;
using namespace CarpetLib;



// Communication state control
comm_state::comm_state ()
{
  // A comm_state object will step through
  //    state_get_buffer_sizes
  //    state_fill_send_buffers
  //    state_empty_recv_buffers

  DECLARE_CCTK_PARAMETERS;

  static Timer timer ("commstate::create");
  timer.start ();
  thestate = state_get_buffer_sizes;

  typebufs.resize (dist::c_ndatatypes());
#define INSTANTIATE(T)                                          \
  {                                                             \
    T dummy;                                                    \
    int const type = dist::c_datatype (dummy);                  \
    assert (typebufs.AT(type).datatypesize == 0);               \
    typebufs.AT(type).datatypesize = sizeof dummy;              \
    typebufs.AT(type).mpi_datatype = dist::datatype (dummy);    \
    typebufs.AT(type).procbufs.resize (dist::size());           \
  }
#include "instantiate"
#undef INSTANTIATE

  srequests.resize (dist::c_ndatatypes() * dist::size(), MPI_REQUEST_NULL);
  rrequests.resize (dist::c_ndatatypes() * dist::size(), MPI_REQUEST_NULL);
  
  timer.stop (0);
}


void comm_state::step ()
{
  DECLARE_CCTK_PARAMETERS;
  static Timer total ("commstate::step");
  total.start ();
  assert (thestate != state_done);
  switch (thestate) {
    
  case state_get_buffer_sizes:
    // The sizes of the collective communication buffers are known so
    // now allocate them.
    // The receive operations are also posted here already (a clever
    // MPI layer may take advantage of such early posting).
    num_posted_recvs = num_completed_recvs = 0;
    
    for (int proc1 = 0; proc1 < dist::size(); ++ proc1) {
      size_t const proc =
        interleave_communications
        ? (proc1 + dist::rank()) % dist::size()
        : proc1;
      
      for (size_t type = 0; type < typebufs.size(); type++) {
        
        // skip unused datatype buffers
        if (not typebufs.AT(type).in_use) continue;
        
        int datatypesize = typebufs.AT(type).datatypesize;
        procbufdesc& procbuf = typebufs.AT(type).procbufs.AT(proc);
        
        assert (procbuf.sendbufbase.empty());
        assert (procbuf.recvbufbase.empty());
        procbuf.sendbufbase.resize (procbuf.sendbufsize*datatypesize);
        procbuf.recvbufbase.resize (procbuf.recvbufsize*datatypesize);
        // TODO: this may be a bit extreme, and it is only for
        // internal consistency checking
        if (poison_new_memory) {
          memset (&procbuf.sendbufbase.front(), poison_value,
                  procbuf.sendbufsize*datatypesize);
          memset (&procbuf.recvbufbase.front(), poison_value,
                  procbuf.recvbufsize*datatypesize);
        }
        procbuf.sendbuf = &procbuf.sendbufbase.front();
        procbuf.recvbuf = &procbuf.recvbufbase.front();
        
        if (procbuf.recvbufsize > 0) {
          static Timer timer ("commstate_sizes_irecv");
          timer.start ();
          int const tag =
            vary_tags
            ? (dist::rank() + dist::size() * (proc + dist::size() * type)) % 32768
            : type;
          if (commstate_verbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "About to MPI_Irecv from %d", (int)proc);
          }
          MPI_Irecv (&procbuf.recvbufbase.front(), procbuf.recvbufsize,
                     typebufs.AT(type).mpi_datatype, proc, tag,
                     dist::comm(), &rrequests.AT(dist::size()*type + proc));
          if (commstate_verbose) {
            CCTK_INFO ("Finished MPI_Irecv");
          }
          timer.stop (procbuf.recvbufsize * datatypesize);
          num_posted_recvs++;
        }
      }
    }
    
    if (barrier_between_stages) {
      // Add a barrier, to try to ensure that all Irecvs are posted
      // before the first Isends are made
      // (Alternative: Use MPI_Alltoallv instead)
      MPI_Barrier (dist::comm());
    }
    
    // Now go and get the send buffers filled with data.
    // Once a buffer is full it will be posted right away
    // (see gdata::copy_into_sendbuffer() and
    //  gdata::interpolate_into_sendbuffer()).
    thestate = state_fill_send_buffers;
    break;
    
  case state_fill_send_buffers:
    if (combine_sends) {
      // Send the data.  Do not send them sequentially, but try to
      // intersperse the communications
      for (int proc1 = 0; proc1 < dist::size(); ++ proc1) {
        int const proc =
          interleave_communications
          ? (proc1 + dist::size() - dist::rank()) % dist::size()
          : proc1;
        
        for (size_t type = 0; type < typebufs.size(); type++) {
          // skip unused datatype buffers
          if (not typebufs.AT(type).in_use) continue;
          
          int const datatypesize = typebufs.AT(type).datatypesize;
          procbufdesc const & procbuf = typebufs.AT(type).procbufs.AT(proc);
          
          size_t const fillstate =
            procbuf.sendbuf - &procbuf.sendbufbase.front();
          assert (fillstate == procbuf.sendbufsize * datatypesize);
          
          if (procbuf.sendbufsize > 0) {
            int const tag =
              vary_tags
              ? (proc + dist::size() * (dist::rank() + dist::size() * type)) % 32768
              : type;
            if (use_mpi_send) {
              // use MPI_Send
              static Timer timer ("commstate_send");
              timer.start ();
              if (commstate_verbose) {
                CCTK_VInfo (CCTK_THORNSTRING,
                            "About to MPI_Send to %d", (int)proc);
              }
              MPI_Send (const_cast<char*>(&procbuf.sendbufbase.front()),
                        procbuf.sendbufsize,
                        typebufs.AT(type).mpi_datatype, proc, tag,
                        dist::comm());
              if (commstate_verbose) {
                CCTK_INFO ("Finished MPI_Send");
              }
              srequests.AT(dist::size()*type + proc) = MPI_REQUEST_NULL;
              timer.stop (procbuf.sendbufsize * datatypesize);
            } else if (use_mpi_ssend) {
              // use MPI_Ssend
              static Timer timer ("commstate_ssend");
              timer.start ();
              if (commstate_verbose) {
                CCTK_VInfo (CCTK_THORNSTRING,
                            "About to MPI_Ssend to %d", (int)proc);
              }
              MPI_Ssend (const_cast<char*>(&procbuf.sendbufbase.front()),
                         procbuf.sendbufsize,
                         typebufs.AT(type).mpi_datatype, proc, tag,
                         dist::comm());
              if (commstate_verbose) {
                CCTK_INFO ("Finished MPI_Ssend");
              }
              srequests.AT(dist::size()*type + proc) = MPI_REQUEST_NULL;
              timer.stop (procbuf.sendbufsize * datatypesize);
            } else {
              // use MPI_Isend
              static Timer timer ("commstate_isend");
              timer.start ();
              if (commstate_verbose) {
                CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                            "About to MPI_Isend to %d", (int)proc);
              }
              MPI_Isend (const_cast<char*>(&procbuf.sendbufbase.front()),
                         procbuf.sendbufsize,
                         typebufs.AT(type).mpi_datatype, proc, tag,
                         dist::comm(), &srequests.AT(dist::size()*type + proc));
              if (commstate_verbose) {
                CCTK_INFO ("Finished MPI_Isend");
              }
              timer.stop (procbuf.sendbufsize * datatypesize);
            }
          }
          
        } // for type
        
      } // for proc
    }
    
    // Now fall through to the next state in which the recv buffers
    // are emptied as soon as data has arrived.
    thestate = state_do_some_work;
    break;
    
  case state_do_some_work:
    // Now fall through to the next state in which the recv buffers
    // are emptied as soon as data has arrived.
    thestate = state_empty_recv_buffers;
    
  case state_empty_recv_buffers:
    // Finish (at least one of) the posted communications
    if (not AllPostedCommunicationsFinished ()) {
      // No state change if there are still outstanding
      // communications; do another comm_state loop iteration.
    } else {
      // Everything is done so release the collective communication buffers.
      for (size_t type = 0; type < typebufs.size(); type++) {
        for (size_t proc = 0; proc < typebufs.AT(type).procbufs.size(); proc++) {
          typebufs.AT(type).procbufs.AT(proc).sendbufbase.clear();
          typebufs.AT(type).procbufs.AT(proc).recvbufbase.clear();
        }
      }
      thestate = state_done;
    }
    break;
    
  default:
    assert (0 && "invalid state");
  }
  total.stop (0);
}



bool comm_state::done ()
{
  return thestate == state_done;
}


comm_state::~comm_state ()
{
  DECLARE_CCTK_PARAMETERS;

  assert (thestate == state_done or
          thestate == state_get_buffer_sizes);
}


// wait for completion of posted collective buffer sends/receives
//
// This function will wait for all of the posted receive operations to
// finish.
//
// It returns true if all posted communications have been completed.
bool comm_state::AllPostedCommunicationsFinished ()
{
  DECLARE_CCTK_PARAMETERS;
 
  // check if all outstanding receives have been completed already
  if (num_posted_recvs == num_completed_recvs) {
    // finalize the outstanding sends in one go
    if (reduce_mpi_waitall) {
      size_t nreqs = 0;
      for (size_t i=0; i<srequests.size(); ++i) {
        if (srequests.AT(i) != MPI_REQUEST_NULL) {
          ++nreqs;
        }
      }
      vector<MPI_Request> reqs(nreqs);
      nreqs = 0;
      for (size_t i=0; i<srequests.size(); ++i) {
        if (srequests.AT(i) != MPI_REQUEST_NULL) {
          reqs.AT(nreqs) = srequests.AT(i);
          ++nreqs;
        }
      }
      assert (nreqs == reqs.size());
      static Timer timer ("commstate_waitall_final");
      timer.start ();
      if (commstate_verbose) {
        CCTK_INFO ("About to MPI_Waitall");
      }
      MPI_Waitall (reqs.size(), &reqs.front(), MPI_STATUSES_IGNORE);
      if (commstate_verbose) {
        CCTK_INFO ("Finished MPI_Waitall");
      }
      timer.stop (0);
    } else {
      static Timer timer ("commstate_waitall_final");
      timer.start ();
      if (commstate_verbose) {
        CCTK_INFO ("About to MPI_Waitall");
      }
      MPI_Waitall (srequests.size(), &srequests.front(), MPI_STATUSES_IGNORE);
      if (commstate_verbose) {
        CCTK_INFO ("Finished MPI_Waitall");
      }
      timer.stop (0);
    }

    return true;
  }

  // wait for completion of all posted receive operations
  if (reduce_mpi_waitall) {
    size_t nreqs = 0;
    for (size_t i=0; i<rrequests.size(); ++i) {
      if (rrequests.AT(i) != MPI_REQUEST_NULL) {
        ++nreqs;
      }
    }
    vector<MPI_Request> reqs(nreqs);
    nreqs = 0;
    for (size_t i=0; i<rrequests.size(); ++i) {
      if (rrequests.AT(i) != MPI_REQUEST_NULL) {
        reqs.AT(nreqs) = rrequests.AT(i);
        ++nreqs;
      }
    }
    assert (nreqs == reqs.size());
    static Timer timer ("commstate_waitall");
    timer.start ();
    if (commstate_verbose) {
      CCTK_INFO ("About to MPI_Waitall");
    }
    MPI_Waitall (reqs.size(), &reqs.front(), MPI_STATUSES_IGNORE);
    if (commstate_verbose) {
      CCTK_INFO ("Finished MPI_Waitall");
    }
    timer.stop (0);
  } else {
    static Timer timer ("commstate_waitall");
    timer.start ();
    if (commstate_verbose) {
      CCTK_INFO ("About to MPI_Waitall");
    }
    MPI_Waitall (rrequests.size(), &rrequests.front(), MPI_STATUSES_IGNORE);
    if (commstate_verbose) {
      CCTK_INFO ("Finished MPI_Waitall");
    }
    timer.stop (0);
  }
  num_completed_recvs = num_posted_recvs;
  
  return false;
}



void
comm_state::
reserve_send_space (unsigned int const type,
                    int const proc,
                    int const npoints)
{
  assert (type < dist::c_ndatatypes());
  assert (proc >= 0 and proc < dist::size());
  assert (npoints >= 0);
  typebufdesc & typebuf = typebufs.AT(type);
  procbufdesc & procbuf = typebuf.procbufs.AT(proc);
  procbuf.sendbufsize += npoints;
  typebuf.in_use = true;
}

void
comm_state::
reserve_recv_space (unsigned int const type,
                    int const proc,
                    int const npoints)
{
  assert (type < dist::c_ndatatypes());
  assert (proc >= 0 and proc < dist::size());
  assert (npoints >= 0);
  typebufdesc & typebuf = typebufs.AT(type);
  procbufdesc & procbuf = typebuf.procbufs.AT(proc);
  procbuf.recvbufsize += npoints;
  typebuf.in_use = true;
}

void *
comm_state::
send_buffer (unsigned int const type,
             int const proc,
             int const npoints)
{
  assert (type < dist::c_ndatatypes());
  assert (proc >= 0 and proc < dist::size());
  typebufdesc const & typebuf = typebufs.AT(type);
  procbufdesc const & procbuf = typebuf.procbufs.AT(proc);
  
  assert (procbuf.sendbuf + npoints * typebuf.datatypesize <=
          &procbuf.sendbufbase.front() +
          procbuf.sendbufsize * typebuf.datatypesize);
  
  return procbuf.sendbuf;
}

void *
comm_state::
recv_buffer (unsigned int const type,
             int const proc,
             int const npoints)
{
  assert (type < dist::c_ndatatypes());
  assert (proc >= 0 and proc < dist::size());
  typebufdesc const & typebuf = typebufs.AT(type);
  procbufdesc const & procbuf = typebuf.procbufs.AT(proc);
  
  assert (procbuf.recvbuf + npoints * typebuf.datatypesize <=
          &procbuf.recvbufbase.front() +
          procbuf.recvbufsize * typebuf.datatypesize);
  
  return procbuf.recvbuf;
}

void
comm_state::
commit_send_space (unsigned int const type,
                   int const proc,
                   int const npoints)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (type < dist::c_ndatatypes());
  assert (proc >= 0 and proc < dist::size());
  assert (npoints >= 0);
  typebufdesc & typebuf = typebufs.AT(type);
  procbufdesc & procbuf = typebuf.procbufs.AT(proc);
  procbuf.sendbuf += npoints * typebuf.datatypesize;
  assert (procbuf.sendbuf <=
          &procbuf.sendbufbase.front() +
          procbuf.sendbufsize * typebuf.datatypesize);
  
  if (not combine_sends) {
    // post the send if the buffer is full
    if (procbuf.sendbuf ==
        &procbuf.sendbufbase.front() +
        procbuf.sendbufsize * typebuf.datatypesize)
    {
      static Timer timer ("commit_send_space::isend");
      timer.start ();
      if (commstate_verbose) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "About to MPI_Isend to %d", (int)proc);
      }
      MPI_Isend (&procbuf.sendbufbase.front(),
                 procbuf.sendbufsize, typebuf.mpi_datatype,
                 proc, type, dist::comm(),
                 & srequests.AT(type * dist::size() + proc));
      if (commstate_verbose) {
        CCTK_INFO ("Finished MPI_Isend");
      }
      timer.stop (procbuf.sendbufsize * typebuf.datatypesize);
    }
  }
}

void
comm_state::
commit_recv_space (unsigned int const type,
                   int const proc,
                   int const npoints)
{
  assert (type < dist::c_ndatatypes());
  assert (proc >= 0 and proc < dist::size());
  assert (npoints >= 0);
  typebufdesc & typebuf = typebufs.AT(type);
  procbufdesc & procbuf = typebuf.procbufs.AT(proc);
  procbuf.recvbuf += npoints * typebuf.datatypesize;
  assert (procbuf.recvbuf <=
          &procbuf.recvbufbase.front() +
          procbuf.recvbufsize * typebuf.datatypesize);
}
