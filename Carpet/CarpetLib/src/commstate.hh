#ifndef COMMSTATE_HH
#define COMMSTATE_HH

#include <queue>
#include <vector>

#include <mpi.h>



using namespace std;



class gdata;



// State information for communications
enum astate { state_recv, state_send, state_wait, state_done };

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
};



#endif  // COMMSTATE_HH
