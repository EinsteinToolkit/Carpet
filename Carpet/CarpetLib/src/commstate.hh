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
  
  queue<gdata*> tmps1, tmps2;
  vector<MPI_Request> requests; // for use_waitall
};



#endif  // COMMSTATE_HH
