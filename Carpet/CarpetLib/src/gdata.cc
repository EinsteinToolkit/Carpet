// $Header:$

#include <assert.h>
#include <stdlib.h>

#include <iostream>

#include "cctk.h"
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



// Constructors
template<int D>
gdata<D>::gdata (const int varindex_, const operator_type transport_operator_)
  : varindex(varindex_), transport_operator(transport_operator_),
    wtime_isend(0.0), wtime_isendwait(0.0),
    wtime_irecv(0.0), wtime_irecvwait(0.0),
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
  
  if (proc() == src->proc()) {
    // copy on same processor
    
  } else {
    
    // copy to different processor
    gdata<D>* const tmp = make_typed(varindex, transport_operator);
    state.tmps1.push (tmp);
    tmp->allocate (box, src->proc());
    tmp->change_processor_recv (state, proc());
    
  }
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
  
  if (proc() == src->proc()) {
    // copy on same processor
    
    copy_from_nocomm (src, box);
    
  } else {
    
    // copy to different processor
    gdata<D>* const tmp = state.tmps1.front();
    state.tmps1.pop();
    state.tmps2.push (tmp);
    assert (tmp);
    tmp->copy_from_nocomm (src, box);
    tmp->change_processor_send (state, proc());
    
  }
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
  
  if (proc() == src->proc()) {
    // copy on same processor
    
  } else {
    
    // copy to different processor
    gdata<D>* const tmp = state.tmps2.front();
    state.tmps2.pop();
    assert (tmp);
    tmp->change_processor_wait (state, proc());
    copy_from_nocomm (tmp, box);
    delete tmp;
    
  }
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
