#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "gdata.hh"

using namespace std;



// Communication state control
template<int D>
comm_state<D>::comm_state ()
  : thestate(state_recv)
{
}

template<int D>
void comm_state<D>::step ()
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (thestate!=state_done);
  if (combine_recv_send) {
    switch (thestate) {
    case state_recv:
      assert (tmps1.empty());
      thestate = state_wait;
      break;
    case state_send:
      assert (0);
    case state_wait:
      assert (tmps1.empty());
      assert (tmps2.empty());
      thestate = state_done;
      break;
    case state_done:
      assert (0);
    default:
      assert (0);
    }
  } else {
    switch (thestate) {
    case state_recv:
      assert (tmps2.empty());
      thestate = state_send;
      break;
    case state_send:
      assert (tmps1.empty());
      thestate = state_wait;
      break;
    case state_wait:
      assert (tmps1.empty());
      assert (tmps2.empty());
      thestate = state_done;
      break;
    case state_done:
      assert (0);
    default:
      assert (0);
    }
  }
}

template<int D>
bool comm_state<D>::done ()
{
  return thestate==state_done;
}

template<int D>
comm_state<D>::~comm_state ()
{
  assert (thestate==state_recv || thestate==state_done);
  assert (tmps1.empty());
  assert (tmps2.empty());
  assert (requests.empty());
}



// Hand out the next MPI tag
static int nexttag ()
{
  DECLARE_CCTK_PARAMETERS;
  
  int const min_tag = 100;
  static int last = 0;
  ++last;
  if (last >= max_mpi_tags) last = 0;
  return min_tag + last;
}



timestat::timestat () : wtime(0.0), wtime2(0.0), count(0.0), running(false)
{
}

void timestat::addstat (double const t)
{
  wtime += t;
  wtime2 += t*t;
  ++count;
}

void timestat::start ()
{
  assert (! running);
  running = true;
  starttime = MPI_Wtime();
}

void timestat::stop ()
{
  assert (running);
  running = false;
  double const endtime = MPI_Wtime();
  addstat (endtime - starttime);
}

ostream& operator<< (ostream& os, const timestat& wt)
{
  double const avg = wt.wtime / wt.count;
  double const stddev = sqrt(max(0.0, wt.wtime2 / wt.count - avg * avg));
  os << "timestat[seconds]:"
     << " cnt: " << wt.count
     << " sum: " << wt.wtime
     << " avg: " << avg
     << " stddev: " << stddev;
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

timestat wtime_changeproc_recv;
timestat wtime_changeproc_send;
timestat wtime_changeproc_wait;

timestat wtime_irecv;
timestat wtime_isend;
timestat wtime_isendwait;
timestat wtime_irecvwait;

extern "C" void CarpetLib_printtimestats (CCTK_ARGUMENTS);
void CarpetLib_printtimestats (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (print_timestats) {
    cout << "Timing statistics from CarpetLib:" << endl
         << "   wtime_copyfrom_recv:   " << wtime_copyfrom_recv   << endl
         << "   wtime_copyfrom_send:   " << wtime_copyfrom_send   << endl
         << "   wtime_copyfrom_wait:   " << wtime_copyfrom_wait   << endl
         << endl
         << "   wtime_copyfrom_recv_maketyped:        " << wtime_copyfrom_recv_maketyped        << endl   
         << "   wtime_copyfrom_recv_allocate:         " << wtime_copyfrom_recv_allocate         << endl
         << "   wtime_copyfrom_recv_changeproc_recv:  " << wtime_copyfrom_recv_changeproc_recv  << endl
         << "   wtime_copyfrom_send_copyfrom_nocomm1: " << wtime_copyfrom_send_copyfrom_nocomm1 << endl
         << "   wtime_copyfrom_send_copyfrom_nocomm2: " << wtime_copyfrom_send_copyfrom_nocomm2 << endl
         << "   wtime_copyfrom_send_changeproc_send:  " << wtime_copyfrom_send_changeproc_send  << endl
         << "   wtime_copyfrom_wait_changeproc_wait:  " << wtime_copyfrom_wait_changeproc_wait  << endl
         << "   wtime_copyfrom_wait_copyfrom_nocomm2: " << wtime_copyfrom_wait_copyfrom_nocomm  << endl
         << "   wtime_copyfrom_wait_delete:           " << wtime_copyfrom_wait_delete           << endl
         << endl
         << "   wtime_changeproc_recv: " << wtime_changeproc_recv << endl
         << "   wtime_changeproc_send: " << wtime_changeproc_send << endl
         << "   wtime_changeproc_wait: " << wtime_changeproc_wait << endl
         << endl
         << "   wtime_irecv:           " << wtime_irecv           << endl
         << "   wtime_isend:           " << wtime_isend           << endl
         << "   wtime_isendwait:       " << wtime_isendwait       << endl
         << "   wtime_irecvwait:       " << wtime_irecvwait       << endl
         << endl;
  }
}



// Constructors
template<int D>
gdata<D>::gdata (const int varindex_, const operator_type transport_operator_)
  : varindex(varindex_), transport_operator(transport_operator_),
    _has_storage(false),
    comm_active(false),
    tag(nexttag())
{
  DECLARE_CCTK_PARAMETERS;
  if (barriers) {
    MPI_Barrier (dist::comm);
  }
}

// Destructors
template<int D>
gdata<D>::~gdata ()
{
  DECLARE_CCTK_PARAMETERS;
  if (barriers) {
    MPI_Barrier (dist::comm);
  }
}



// Processor management
template<int D>
void gdata<D>::change_processor (comm_state<D>& state,
                                 const int newproc,
                                 void* const mem)
{
  DECLARE_CCTK_PARAMETERS;
  
  switch (state.thestate) {
  case state_recv:
    if (combine_recv_send) {
      change_processor_recv (state, newproc, mem);
      change_processor_send (state, newproc, mem);
    } else {
      change_processor_recv (state, newproc, mem);
    }
    break;
  case state_send:
    if (combine_recv_send) {
      // do nothing
    } else {
      change_processor_send (state, newproc, mem);
    }
    break;
  case state_wait:
    change_processor_wait (state, newproc, mem);
    break;
  default:
    assert(0);
  }
}



// Data manipulators
template<int D>
void gdata<D>::copy_from (comm_state<D>& state,
                          const gdata* src, const ibbox& box)
{
  DECLARE_CCTK_PARAMETERS;
  
  switch (state.thestate) {
  case state_recv:
    if (combine_recv_send) {
      copy_from_recv (state, src, box);
      copy_from_send (state, src, box);
    } else {
      copy_from_recv (state, src, box);
    }
    break;
  case state_send:
    if (combine_recv_send) {
      // do nothing
    } else {
      copy_from_send (state, src, box);
    }
    break;
  case state_wait:
    copy_from_wait (state, src, box);
    break;
  default:
    assert(0);
  }
}



template<int D>
void gdata<D>::copy_from_nocomm (const gdata* src, const ibbox& box)
{
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()
	      && box.stride()==src->extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      && (box.lower()-src->extent().lower())%box.stride() == 0));
  
  if (box.empty()) return;
  
  assert (has_storage() && src->has_storage());
  assert (proc() == src->proc());
  
  // copy on same processor
  if (lives_on_this_processor()) {
    copy_from_innerloop (src, box);
  }
}



template<int D>
void gdata<D>::copy_from_recv (comm_state<D>& state,
                               const gdata* src, const ibbox& box)
{
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()
	      && box.stride()==src->extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      && (box.lower()-src->extent().lower())%box.stride() == 0));
  
  if (box.empty()) return;
  
  wtime_copyfrom_recv.start();
  
  if (proc() == src->proc()) {
    // copy on same processor
    
  } else {
    
    // copy to different processor
    wtime_copyfrom_recv_maketyped.start();
    gdata<D>* const tmp = make_typed(varindex, transport_operator);
    wtime_copyfrom_recv_maketyped.stop();
    state.tmps1.push (tmp);
    wtime_copyfrom_recv_allocate.start();
    tmp->allocate (box, src->proc());
    wtime_copyfrom_recv_allocate.stop();
    wtime_copyfrom_recv_changeproc_recv.start();
    tmp->change_processor_recv (state, proc());
    wtime_copyfrom_recv_changeproc_recv.stop();
    
  }
  
  wtime_copyfrom_recv.stop();
}



template<int D>
void gdata<D>::copy_from_send (comm_state<D>& state,
                               const gdata* src, const ibbox& box)
{
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()
	      && box.stride()==src->extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      && (box.lower()-src->extent().lower())%box.stride() == 0));
  
  if (box.empty()) return;
  
  wtime_copyfrom_send.start();
  
  if (proc() == src->proc()) {
    // copy on same processor
    
    wtime_copyfrom_send_copyfrom_nocomm1.start();
    copy_from_nocomm (src, box);
    wtime_copyfrom_send_copyfrom_nocomm1.stop();
    
  } else {
    
    // copy to different processor
    gdata<D>* const tmp = state.tmps1.front();
    state.tmps1.pop();
    state.tmps2.push (tmp);
    assert (tmp);
    wtime_copyfrom_send_copyfrom_nocomm2.start();
    tmp->copy_from_nocomm (src, box);
    wtime_copyfrom_send_copyfrom_nocomm2.stop();
    wtime_copyfrom_send_changeproc_send.start();
    tmp->change_processor_send (state, proc());
    wtime_copyfrom_send_changeproc_send.stop();
    
  }
  
  wtime_copyfrom_send.stop();
}



template<int D>
void gdata<D>::copy_from_wait (comm_state<D>& state,
                               const gdata* src, const ibbox& box)
{
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()
	      && box.stride()==src->extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      && (box.lower()-src->extent().lower())%box.stride() == 0));
  
  if (box.empty()) return;
  
  wtime_copyfrom_wait.start();
  
  if (proc() == src->proc()) {
    // copy on same processor
    
  } else {
    
    // copy to different processor
    gdata<D>* const tmp = state.tmps2.front();
    state.tmps2.pop();
    assert (tmp);
    wtime_copyfrom_wait_changeproc_wait.start();
    tmp->change_processor_wait (state, proc());
    wtime_copyfrom_wait_changeproc_wait.stop();
    wtime_copyfrom_wait_copyfrom_nocomm.start();
    copy_from_nocomm (tmp, box);
    wtime_copyfrom_wait_copyfrom_nocomm.stop();
    wtime_copyfrom_wait_delete.start();
    delete tmp;
    wtime_copyfrom_wait_delete.stop();
    
  }
  
  wtime_copyfrom_wait.stop();
}



template<int D>
void gdata<D>
::interpolate_from (comm_state<D>& state,
                    const vector<const gdata*> srcs,
                    const vector<CCTK_REAL> times,
                    const ibbox& box, const CCTK_REAL time,
                    const int order_space,
                    const int order_time)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (transport_operator != op_error);
  if (transport_operator == op_none) return;
  switch (state.thestate) {
  case state_recv:
    if (combine_recv_send) {
      interpolate_from_recv (state, srcs, times, box, time, order_space, order_time);
      interpolate_from_send (state, srcs, times, box, time, order_space, order_time);
    } else {
      interpolate_from_recv (state, srcs, times, box, time, order_space, order_time);
    }
    break;
  case state_send:
    if (combine_recv_send) {
      // do nothing
    } else {
      interpolate_from_send (state, srcs, times, box, time, order_space, order_time);
    }
    break;
  case state_wait:
    interpolate_from_wait (state, srcs, times, box, time, order_space, order_time);
    break;
  default:
    assert(0);
  }
}



template<int D>
void gdata<D>
::interpolate_from_nocomm (const vector<const gdata*> srcs,
                           const vector<CCTK_REAL> times,
                           const ibbox& box, const CCTK_REAL time,
                           const int order_space,
                           const int order_time)
{
  assert (has_storage());
  assert (all(box.lower()>=extent().lower()));
  assert (all(box.upper()<=extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  assert (srcs.size() == times.size() && srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs.at(t)->has_storage());
    assert (all(box.lower()>=srcs.at(t)->extent().lower()));
    assert (all(box.upper()<=srcs.at(t)->extent().upper()));
  }
  
  assert (! box.empty());
  if (box.empty()) return;
  
  assert (proc() == srcs.at(0)->proc());
  
  assert (transport_operator != op_error);
  assert (transport_operator != op_none);
  
  // interpolate on same processor
  if (lives_on_this_processor()) {
    interpolate_from_innerloop
      (srcs, times, box, time, order_space, order_time);
  }
}



template<int D>
void gdata<D>
::interpolate_from_recv (comm_state<D>& state,
                         const vector<const gdata*> srcs,
                         const vector<CCTK_REAL> times,
                         const ibbox& box, const CCTK_REAL time,
                         const int order_space,
                         const int order_time)
{
  assert (has_storage());
  assert (all(box.lower()>=extent().lower()));
  assert (all(box.upper()<=extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  assert (srcs.size() == times.size() && srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs.at(t)->has_storage());
    assert (all(box.lower()>=srcs.at(t)->extent().lower()));
    assert (all(box.upper()<=srcs.at(t)->extent().upper()));
  }
  
  assert (! box.empty());
  if (box.empty()) return;
  
  if (proc() == srcs.at(0)->proc()) {
    // interpolate on same processor
    
  } else {
    // interpolate from other processor
    
    gdata<D>* const tmp = make_typed(varindex, transport_operator);
    state.tmps1.push (tmp);
    tmp->allocate (box, srcs.at(0)->proc());
    tmp->change_processor_recv (state, proc());
    
  }
}



template<int D>
void gdata<D>
::interpolate_from_send (comm_state<D>& state,
                         const vector<const gdata*> srcs,
                         const vector<CCTK_REAL> times,
                         const ibbox& box, const CCTK_REAL time,
                         const int order_space,
                         const int order_time)
{
  assert (has_storage());
  assert (all(box.lower()>=extent().lower()));
  assert (all(box.upper()<=extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  assert (srcs.size() == times.size() && srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs.at(t)->has_storage());
    assert (all(box.lower()>=srcs.at(t)->extent().lower()));
    assert (all(box.upper()<=srcs.at(t)->extent().upper()));
  }
  
  assert (! box.empty());
  if (box.empty()) return;
  
  if (proc() == srcs.at(0)->proc()) {
    // interpolate on same processor
    
    interpolate_from_nocomm (srcs, times, box, time, order_space, order_time);
    
  } else {
    // interpolate from other processor
    
    gdata<D>* const tmp = state.tmps1.front();
    state.tmps1.pop();
    state.tmps2.push (tmp);
    assert (tmp);
    tmp->interpolate_from_nocomm (srcs, times, box, time, order_space, order_time);
    tmp->change_processor_send (state, proc());
    
  }
}



template<int D>
void gdata<D>
::interpolate_from_wait (comm_state<D>& state,
                         const vector<const gdata*> srcs,
                         const vector<CCTK_REAL> times,
                         const ibbox& box, const CCTK_REAL time,
                         const int order_space,
                         const int order_time)
{
  assert (has_storage());
  assert (all(box.lower()>=extent().lower()));
  assert (all(box.upper()<=extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  assert (srcs.size() == times.size() && srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs.at(t)->has_storage());
    assert (all(box.lower()>=srcs.at(t)->extent().lower()));
    assert (all(box.upper()<=srcs.at(t)->extent().upper()));
  }
  
  assert (! box.empty());
  if (box.empty()) return;
  
  if (proc() == srcs.at(0)->proc()) {
    // interpolate on same processor
    
  } else {
    // interpolate from other processor
    
    gdata<D>* const tmp = state.tmps2.front();
    state.tmps2.pop();
    assert (tmp);
    tmp->change_processor_wait (state, proc());
    copy_from_nocomm (tmp, box);
    delete tmp;
    
  }
}


template<int D>
bool gdata<D>
::this_processor_is (int procno)
{
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  return rank == procno;
}

template<int D>
bool gdata<D>
::lives_on_this_processor ()
{
  return this_processor_is( proc() );
}

template<int D>
gdata<D>& gdata<D>
::operator = ( const gdata & ) // canonical copy
{
	return *this; // does nothing
}

template class comm_state<3>;
template class gdata<3>;
