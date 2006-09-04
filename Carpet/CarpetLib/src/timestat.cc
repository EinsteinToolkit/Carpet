#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "timestat.hh"



using namespace std;



timestat::timestat ()
  : wtime(0.0), wtime2(0.0), wmin(0.0), wmax(0.0),
    bytes(0.0), bytes2(0.0), bmin(0.0), bmax(0.0),
    count(0.0),
    running(false)
{
}

void timestat::addstat (double const t, double const b)
{
  wtime += t;
  wtime2 += t*t;
  wmin = min (wmin, t);
  wmax = max (wmax, t);
  
  bytes += b;
  bytes2 += b*b;
  bmin = min (bmin, b);
  bmax = max (bmax, b);
  
  ++count;
}

void timestat::start ()
{
  assert (! running);
  running = true;
  starttime = MPI_Wtime();
}

void timestat::stop (double const b)
{
  assert (running);
  running = false;
  double const endtime = MPI_Wtime();
  addstat (endtime - starttime, b);
}



ostream& operator<< (ostream& os, const timestat& wt)
{
  double const avg = wt.wtime / wt.count;
  double const stddev = sqrt (max (0.0, wt.wtime2 / wt.count - pow (avg, 2)));
  double const bavg = wt.bytes / wt.count;
  double const bstddev = sqrt (max (0.0, wt.bytes2 / wt.count - pow (bavg, 2)));
  os << "timestat[seconds]:"
     << " cnt: " << wt.count
     << "   time: sum: " << wt.wtime
     << " avg: " << avg
     << " stddev: " << stddev
     << " min: " << wt.wmin
     << " max: " << wt.wmax
     << "   bytes: sum: " << wt.bytes
     << " avg: " << bavg
     << " stddev: " << bstddev
     << " min: " << wt.bmin
     << " max: " << wt.bmax;
  return os;
}



timestat wtime_copyfrom_recv;
timestat wtime_copyfrom_send;
timestat wtime_copyfrom_wait;

timestat wtime_copyfrom_recv_maketyped;
timestat wtime_copyfrom_recv_allocate;
timestat wtime_copyfrom_recv_changeproc_recv;
timestat wtime_copyfrom_send_copyfrom_nocomm1;
timestat wtime_copyfrom_send_copyfrom_nocomm2;
timestat wtime_copyfrom_send_changeproc_send;
timestat wtime_copyfrom_wait_changeproc_wait;
timestat wtime_copyfrom_wait_copyfrom_nocomm;
timestat wtime_copyfrom_wait_delete;

timestat wtime_copyfrom_recvinner_allocate;
timestat wtime_copyfrom_recvinner_recv;
timestat wtime_copyfrom_sendinner_allocate;
timestat wtime_copyfrom_sendinner_copy;
timestat wtime_copyfrom_sendinner_send;
timestat wtime_copyfrom_recvwaitinner_wait;
timestat wtime_copyfrom_recvwaitinner_copy;
timestat wtime_copyfrom_recvwaitinner_delete;
timestat wtime_copyfrom_sendwaitinner_wait;
timestat wtime_copyfrom_sendwaitinner_delete;

timestat wtime_changeproc_recv;
timestat wtime_changeproc_send;
timestat wtime_changeproc_wait;

timestat wtime_irecv;
timestat wtime_isend;
timestat wtime_isendwait;
timestat wtime_irecvwait;

timestat wtime_commstate_sizes_irecv;
timestat wtime_commstate_waitall_final;
timestat wtime_commstate_waitall;
timestat wtime_commstate_waitsome;
timestat wtime_commstate_isend;
timestat wtime_commstate_memcpy;
timestat wtime_commstate_interpolate_irecv;
timestat wtime_commstate_interpolate_from_isend;
timestat wtime_commstate_interpolate_to_isend;



extern "C" void CarpetLib_printtimestats (CCTK_ARGUMENTS);

void CarpetLib_printtimestats (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (print_timestats or
      (print_timestats_every and
       cctk_iteration % print_timestats_every == 0))
  {
    cout << "Timing statistics from CarpetLib:" << endl
         << "   wtime_copyfrom_recv:                    " << wtime_copyfrom_recv                  << endl
         << "   wtime_copyfrom_send:                    " << wtime_copyfrom_send                  << endl
         << "   wtime_copyfrom_wait:                    " << wtime_copyfrom_wait                  << endl
         << endl
         << "   wtime_copyfrom_recv_maketyped:          " << wtime_copyfrom_recv_maketyped        << endl   
         << "   wtime_copyfrom_recv_allocate:           " << wtime_copyfrom_recv_allocate         << endl
         << "   wtime_copyfrom_recv_changeproc_recv:    " << wtime_copyfrom_recv_changeproc_recv  << endl
         << "   wtime_copyfrom_send_copyfrom_nocomm1:   " << wtime_copyfrom_send_copyfrom_nocomm1 << endl
         << "   wtime_copyfrom_send_copyfrom_nocomm2:   " << wtime_copyfrom_send_copyfrom_nocomm2 << endl
         << "   wtime_copyfrom_send_changeproc_send:    " << wtime_copyfrom_send_changeproc_send  << endl
         << "   wtime_copyfrom_wait_changeproc_wait:    " << wtime_copyfrom_wait_changeproc_wait  << endl
         << "   wtime_copyfrom_wait_copyfrom_nocomm2:   " << wtime_copyfrom_wait_copyfrom_nocomm  << endl
         << "   wtime_copyfrom_wait_delete:             " << wtime_copyfrom_wait_delete           << endl
         << endl
         << "   wtime_copyfrom_recvinner_allocate:      " << wtime_copyfrom_recvinner_allocate    << endl
         << "   wtime_copyfrom_recvinner_recv:          " << wtime_copyfrom_recvinner_recv        << endl
         << "   wtime_copyfrom_sendinner_allocate:      " << wtime_copyfrom_sendinner_allocate    << endl
         << "   wtime_copyfrom_sendinner_copy:          " << wtime_copyfrom_sendinner_copy        << endl
         << "   wtime_copyfrom_sendinner_send:          " << wtime_copyfrom_sendinner_send        << endl
         << "   wtime_copyfrom_recvwaitinner_wait:      " << wtime_copyfrom_recvwaitinner_wait    << endl
         << "   wtime_copyfrom_recvwaitinner_copy:      " << wtime_copyfrom_recvwaitinner_copy    << endl
         << "   wtime_copyfrom_recvwaitinner_delete:    " << wtime_copyfrom_recvwaitinner_delete  << endl
         << "   wtime_copyfrom_sendwaitinner_wait:      " << wtime_copyfrom_sendwaitinner_wait    << endl
         << "   wtime_copyfrom_sendwaitinner_delete:    " << wtime_copyfrom_sendwaitinner_delete  << endl
         << endl
         << "   wtime_changeproc_recv:                  " << wtime_changeproc_recv                << endl
         << "   wtime_changeproc_send:                  " << wtime_changeproc_send                << endl
         << "   wtime_changeproc_wait:                  " << wtime_changeproc_wait                << endl
         << endl
         << "   wtime_irecv:                            " << wtime_irecv                          << endl
         << "   wtime_isend:                            " << wtime_isend                          << endl
         << "   wtime_isendwait:                        " << wtime_isendwait                      << endl
         << "   wtime_irecvwait:                        " << wtime_irecvwait                      << endl
         << endl
         << "   wtime_commstate_sizes_irecv:            " << wtime_commstate_sizes_irecv            << endl
         << "   wtime_commstate_waitall_final:          " << wtime_commstate_waitall_final          << endl
         << "   wtime_commstate_waitall:                " << wtime_commstate_waitall                << endl
         << "   wtime_commstate_waitsome:               " << wtime_commstate_waitsome               << endl
         << "   wtime_commstate_isend:                  " << wtime_commstate_isend                  << endl
         << "   wtime_commstate_memcpy:                 " << wtime_commstate_memcpy                 << endl
         << "   wtime_commstate_interpolate_irecv:      " << wtime_commstate_interpolate_irecv      << endl
         << "   wtime_commstate_interpolate_from_isend: " << wtime_commstate_interpolate_from_isend << endl
         << "   wtime_commstate_interpolate_to_isend:   " << wtime_commstate_interpolate_to_isend   << endl
         << endl;
  }
}
