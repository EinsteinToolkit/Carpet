#ifndef COMMSTATE_HH
#define COMMSTATE_HH

#include <queue>
#include <vector>

#include <mpi.h>

#include "dist.hh"

using namespace std;


// State information for communications
//
// Depending on how a comm state object was created,
// it will step through one of two state transitions (in the given order):
enum astate {
  // these are used for collective communications
  state_get_buffer_sizes, state_fill_send_buffers, state_empty_recv_buffers,

  // these are used for communications on individual components
  state_post, state_wait,

  // all transition graphs must end with here
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
  // the following members are used for single-component communications
  //////////////////////////////////////////////////////////////////////////

  // List of MPI requests for use_waitall
  vector<MPI_Request> requests;

  // Lists of communication buffers for use_lightweight_buffers
  struct gcommbuf {
    gcommbuf () {};
    virtual ~gcommbuf () {};
    MPI_Request request;
    virtual void const * pointer () const = 0;
    virtual void * pointer () = 0;
    virtual int size () const = 0;
    virtual MPI_Datatype datatype () const = 0;
  };

  template<typename T>
  struct commbuf : gcommbuf {
    commbuf (ibbox const & box);
    virtual ~commbuf ();
    virtual void const * pointer () const;
    virtual void * pointer ();
    virtual int size () const;
    virtual MPI_Datatype datatype () const;
    vector<T> data;
  };
  queue<gcommbuf*> recvbufs, sendbufs;


  //////////////////////////////////////////////////////////////////////////
  // the following members are used for collective communications
  //////////////////////////////////////////////////////////////////////////

public:
  // structure describing a per-processor buffer for collective communications
  struct procbufdesc {
    // the allocated communication buffers
    char* sendbufbase;
    char* recvbufbase;

    // the sizes of communication buffers (in elements of type <datatype>)
    size_t sendbufsize;
    size_t recvbufsize;

    // pointers to step through the communication buffers
    // (these get advanced by the routines which fill/empty the buffers)
    char* sendbuf;
    char* recvbuf;

    // constructor for an instance of this structure
    procbufdesc() : sendbufbase(NULL), recvbufbase(NULL),
                    sendbufsize(0), recvbufsize(0),
                    sendbuf(NULL), recvbuf(NULL) { }
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
      procbufs.resize (dist::size());
    }
  };

  // list of datatype buffers
  vector<typebufdesc> typebufs;        // [dist::c_ndatatypes()]

  // flags indicating which receive buffers are ready to be emptied
  vector<bool> recvbuffers_ready;      // [dist::size() * dist::c_ndatatypes()]

  // lists of outstanding requests for posted send/recv communications
  vector<MPI_Request> srequests;       // [dist::size() * dist::c_ndatatypes()]
private:
  vector<MPI_Request> rrequests;       // [dist::size() * dist::c_ndatatypes()]

  // number of posted and already completed receive communications
  int num_posted_recvs;
  int num_completed_recvs;

  // wait for completion of posted collective buffer sends/receives
  bool AllPostedCommunicationsFinished();
};



#endif  // COMMSTATE_HH
