#ifndef COMMSTATE_HH
#define COMMSTATE_HH

#include <queue>
#include <vector>

#include <mpi.h>



using namespace std;



class gdata;



// State information for communications
//
// Depending on how a comm state object was created,
// it will step through one of two state transitions (in the given order):
enum astate {
  // "recv -> send -> wait" are used for single-component communications
  state_recv, state_send, state_wait,

  // these are used for collective communications
  state_get_buffer_sizes, state_fill_send_buffers, state_empty_recv_buffers,

  // all transition graphs must end with here
  state_done
};

struct comm_state {
  astate thestate;

  // If no vartype is given the comm state will be set up
  // for single-component communications.
  comm_state (int vartype = -1);
  void step ();
  bool done ();
  ~comm_state ();
  
private:
  // Forbid copying and passing by value
  comm_state (comm_state const &);
  comm_state& operator= (comm_state const &);

  // flag to indicate whether this comm state object is used for
  // collective (true) or single-component communications (false)
  bool use_collective_communication_buffers;

public:
  
  //////////////////////////////////////////////////////////////////////////
  // the following members are used for single-component communications
  //////////////////////////////////////////////////////////////////////////

  // Lists of temporary data objects
  queue<gdata*> tmps1, tmps2;
  
  // List of MPI requests for use_waitall
  vector<MPI_Request> requests;
  
  // Lists of communication buffers for use_lightweight_buffers
  struct gcommbuf {
    gcommbuf ();
    virtual ~gcommbuf ();
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

  // CCTK vartype used for this comm_state object
  int vartype;

  // size of CCTK vartype
  // (used as stride for advancing the char-based buffer pointers)
  int vartypesize;

  // MPI datatype corresponding to CCTK vartype
  MPI_Datatype datatype;

  // buffers for collective communications
  struct collbufdesc {
    // the sizes of communication buffers (in elements of type <vartype>)
    size_t sendbufsize;
    size_t recvbufsize;

    // pointers to step through the communication buffers
    // (these get advanced by the routines which fill/empty the buffers)
    char* sendbuf;
    char* recvbuf;

    collbufdesc() : sendbufsize(0), recvbufsize(0),
                    sendbuf(NULL), recvbuf(NULL),
                    sendbufbase(NULL), recvbufbase(NULL) {}

// FIXME: why can't these be made private ??
//private:
    // the allocated communication buffers
    char* sendbufbase;
    char* recvbufbase;
  };
  vector<collbufdesc> collbufs;          // [nprocs]

private:
  // Exchange pairs of send/recv buffers between all processors.
  void ExchangeBuffers();
};



#endif  // COMMSTATE_HH
