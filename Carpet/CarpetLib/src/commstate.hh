#ifndef COMMSTATE_HH
#define COMMSTATE_HH

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef CCTK_MPI
#  include <mpi.h>
#else
#  include "nompi.h"
#endif

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

char const * tostring (astate const & thestate);

inline ostream& operator<< (ostream& os, astate const & thestate)
{
  return os << tostring(thestate);
}



struct comm_state {
  astate thestate;
  
  comm_state ();
  void step ();
  bool done () const;
  ~comm_state ();
  
  void init();

private:
  // Forbid copying and passing by value
  comm_state (comm_state const &);
  comm_state& operator= (comm_state const &);
  
  
  
  // structure describing a per-processor buffer
  struct procbufdesc {
    // allocated communication buffers
    vector<char> sendbufbase;
    vector<char> recvbufbase;
    
    // sizes of the communication buffers (in elements of type <datatype>)
    size_t sendbufsize;
    size_t recvbufsize;
    
    // pointers to step through the communication buffers
    // (these get advanced by the routines which fill/empty the buffers)
    char* sendbuf;
    char* recvbuf;
    
    bool did_post_send;
    bool did_post_recv;
    
    // constructor for an instance of this structure
    procbufdesc()
    {
        init();
    }

    void init(){
        did_post_send = false;
        did_post_recv = false;
        sendbuf = NULL;
        recvbuf = NULL;
        sendbufsize = 0;
        recvbufsize = 0;
    }

  };
  
  
  
  // structure describing a collective communications buffer for a C datatype
  struct typebufdesc {
    // flag indicating whether this buffer is in use
    bool in_use;
    
    // the MPI datatype
    MPI_Datatype mpi_datatype;
    
    // the size of this datatype (in bytes)
    int datatypesize;
    
    // per-processor buffers
    vector<procbufdesc> procbufs; // [dist::size()]
    
    // constructor for an instance of this structure
    typebufdesc()
    {
        init();
    }

    void init(){
        in_use = false;
        mpi_datatype = MPI_DATATYPE_NULL;
        datatypesize = 0;
    }
  };
  
  
  
  // datatype buffers
  vector<typebufdesc> typebufs; // [type]
  
  
  
  // outstanding requests for posted send/recv communications
  vector<MPI_Request> srequests;
  vector<MPI_Request> rrequests;
  
  static inline
  MPI_Request & push_back (vector<MPI_Request> & reqs)
  {
    reqs.push_back (MPI_REQUEST_NULL);
    return reqs.back();
  }
  
  
  
public:
  
  void
  reserve_send_space (unsigned type,
                      int proc,
                      int npoints);
  
  void
  reserve_recv_space (unsigned type,
                      int proc,
                      int npoints);
  
  void *
  send_buffer (unsigned type,
               int proc,
               int npoints);
  
  void *
  recv_buffer (unsigned type,
               int proc,
               int npoints);
  
  void
  commit_send_space (unsigned type,
                     int proc,
                     int npoints);
  
  void
  commit_recv_space (unsigned type,
                     int proc,
                     int npoints);
};

#endif  // COMMSTATE_HH
