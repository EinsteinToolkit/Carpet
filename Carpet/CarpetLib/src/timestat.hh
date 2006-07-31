#ifndef TIMESTAT_HH
#define TIMESTAT_HH

#include <iostream>

using namespace std;



// Time (in seconds) spend during various operations
class timestat {
  
public:
  double wtime;
  double wtime2;
  double wmin;
  double wmax;
  double count;
  
public:
  timestat ();
  
private:
  void addstat (double const t);
  
private:
  bool running;
  double starttime;
  
public:
  void start();
  void stop();
};

ostream& operator<< (ostream& os, const timestat& wt);



extern timestat wtime_copyfrom_recv;
extern timestat wtime_copyfrom_send;
extern timestat wtime_copyfrom_wait;

extern timestat wtime_copyfrom_recv_maketyped;
extern timestat wtime_copyfrom_recv_allocate;
extern timestat wtime_copyfrom_recv_changeproc_recv;
extern timestat wtime_copyfrom_send_copyfrom_nocomm1;
extern timestat wtime_copyfrom_send_copyfrom_nocomm2;
extern timestat wtime_copyfrom_send_changeproc_send;
extern timestat wtime_copyfrom_wait_changeproc_wait;
extern timestat wtime_copyfrom_wait_copyfrom_nocomm;
extern timestat wtime_copyfrom_wait_delete;

extern timestat wtime_copyfrom_recvinner_allocate;
extern timestat wtime_copyfrom_recvinner_recv;
extern timestat wtime_copyfrom_sendinner_allocate;
extern timestat wtime_copyfrom_sendinner_copy;
extern timestat wtime_copyfrom_sendinner_send;
extern timestat wtime_copyfrom_recvwaitinner_wait;
extern timestat wtime_copyfrom_recvwaitinner_copy;
extern timestat wtime_copyfrom_recvwaitinner_delete;
extern timestat wtime_copyfrom_sendwaitinner_wait;
extern timestat wtime_copyfrom_sendwaitinner_delete;

extern timestat wtime_changeproc_recv;
extern timestat wtime_changeproc_send;
extern timestat wtime_changeproc_wait;

extern timestat wtime_irecv;
extern timestat wtime_isend;
extern timestat wtime_irecvwait;
extern timestat wtime_isendwait;

extern timestat wtime_commstate_sizes_irecv;
extern timestat wtime_commstate_waitall_final;
extern timestat wtime_commstate_waitall;
extern timestat wtime_commstate_waitsome;
extern timestat wtime_commstate_isend;
extern timestat wtime_commstate_memcpy;
extern timestat wtime_commstate_interpolate_irecv;
extern timestat wtime_commstate_interpolate_from_isend;
extern timestat wtime_commstate_interpolate_to_isend;

#endif  // TIMESTAT_HH
