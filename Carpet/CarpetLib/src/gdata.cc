#include <cassert>
#include <cstdlib>
#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "bbox.hh"
#include "commstate.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "gdata.hh"

using namespace std;



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
gdata::gdata (const int varindex_, const operator_type transport_operator_)
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
gdata::~gdata ()
{
  DECLARE_CCTK_PARAMETERS;
  if (barriers) {
    MPI_Barrier (dist::comm);
  }
}

// Assignment
gdata & gdata::operator= (gdata const & from)
{
  return * this = from;
}



// Processor management
void gdata::change_processor (comm_state& state,
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
void gdata::copy_from (comm_state& state,
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



void gdata::copy_from_nocomm (const gdata* src, const ibbox& box)
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
  if (dist::rank() == proc()) {
    copy_from_innerloop (src, box);
  }
}



void gdata::copy_from_recv (comm_state& state,
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
    gdata* const tmp = make_typed(varindex, transport_operator);
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



void gdata::copy_from_send (comm_state& state,
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
    gdata* const tmp = state.tmps1.front();
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



void gdata::copy_from_wait (comm_state& state,
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
    gdata* const tmp = state.tmps2.front();
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



void gdata
::interpolate_from (comm_state& state,
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



void gdata
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
  if (dist::rank() == proc()) {
    interpolate_from_innerloop
      (srcs, times, box, time, order_space, order_time);
  }
}



void gdata
::interpolate_from_recv (comm_state& state,
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
    
    gdata* const tmp = make_typed(varindex, transport_operator);
    state.tmps1.push (tmp);
    tmp->allocate (box, srcs.at(0)->proc());
    tmp->change_processor_recv (state, proc());
    
  }
}



void gdata
::interpolate_from_send (comm_state& state,
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
    
    gdata* const tmp = state.tmps1.front();
    state.tmps1.pop();
    state.tmps2.push (tmp);
    assert (tmp);
    tmp->interpolate_from_nocomm (srcs, times, box, time, order_space, order_time);
    tmp->change_processor_send (state, proc());
    
  }
}



void gdata
::interpolate_from_wait (comm_state& state,
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
    
    gdata* const tmp = state.tmps2.front();
    state.tmps2.pop();
    assert (tmp);
    tmp->change_processor_wait (state, proc());
    copy_from_nocomm (tmp, box);
    delete tmp;
    
  }
}
