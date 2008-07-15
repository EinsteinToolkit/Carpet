#ifndef COMMSTATE_HH
#define COMMSTATE_HH

#include <cstdlib>
#include <queue>
#include <vector>

#include <mpi.h>

#include "cctk_Parameters.h"

#include "dist.hh"
#include "timestat.hh"

using namespace std;
using namespace CarpetLib;



// State information for communications

// A comm state object will step through the state transitions in the
// given order:
enum astate {
  state_get_buffer_sizes,
  state_fill_send_buffers,
  state_do_some_work,
  state_empty_recv_buffers,
  state_done
};

struct comm_state {
  astate thestate;

  comm_state ();
  void step ();
  bool done ();
  ~comm_state ();

private:
  // Forbid copying and passing by value
  comm_state (comm_state const &);
  comm_state& operator= (comm_state const &);

public:

  //////////////////////////////////////////////////////////////////////////
  // the following members are used for collective communications
  //////////////////////////////////////////////////////////////////////////

public:
  // structure describing a per-processor buffer for collective communications
  struct procbufdesc {
    // the allocated communication buffers
    vector<char> sendbufbase;
    vector<char> recvbufbase;

    // the sizes of communication buffers (in elements of type <datatype>)
    size_t sendbufsize;
    size_t recvbufsize;

    // pointers to step through the communication buffers
    // (these get advanced by the routines which fill/empty the buffers)
    char* sendbuf;
    char* recvbuf;

    // constructor for an instance of this structure
    procbufdesc() : sendbufsize(0), recvbufsize(0),
                    sendbuf(NULL), recvbuf(NULL)
    {
    }
  };

  // structure describing a collective communications buffer for a C datatype
  struct typebufdesc {
    // flag indicating whether this buffer is in use
    bool in_use;

    // the size of this datatype (in bytes)
    int datatypesize;

    // the corresponding MPI datatype
    MPI_Datatype mpi_datatype;

    // per-processor buffers
    vector<procbufdesc> procbufs;      // [dist::size()]

    // constructor for an instance of this structure
    typebufdesc() : in_use(false), datatypesize(0),
                    mpi_datatype(MPI_DATATYPE_NULL)
    {
    }
  };

  // list of datatype buffers
  vector<typebufdesc> typebufs;        // [dist::c_ndatatypes()]

  void
  reserve_send_space (unsigned int type,
                      int proc,
                      int npoints);

  void
  reserve_recv_space (unsigned int type,
                      int proc,
                      int npoints);

  void *
  send_buffer (unsigned int type,
               int proc,
               int npoints);

  void *
  recv_buffer (unsigned int type,
               int proc,
               int npoints);

  void
  commit_send_space (unsigned int type,
                     int proc,
                     int npoints);

  void
  commit_recv_space (unsigned int type,
                     int proc,
                     int npoints);

private:
  // lists of outstanding requests for posted send/recv communications
  vector<MPI_Request> srequests;       // [dist::size() * dist::c_ndatatypes()]
  vector<MPI_Request> rrequests;       // [dist::size() * dist::c_ndatatypes()]

  // number of posted and already completed receive communications
  int num_posted_recvs;
  int num_completed_recvs;

  // wait for completion of posted collective buffer sends/receives
  bool AllPostedCommunicationsFinished();
};



#endif  // COMMSTATE_HH
