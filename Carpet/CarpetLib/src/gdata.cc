#include <cassert>
#include <cstdlib>
#include <iostream>

#include <mpi.h>

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
gdata::gdata (const int varindex_,
              const operator_type transport_operator_,
              const int tag_)
  : varindex(varindex_),
    transport_operator(transport_operator_),
    _has_storage(false),
    comm_active(false),
    tag(tag_ >= 0 ? tag_ : nexttag())
{
  DECLARE_CCTK_PARAMETERS;
  if (barriers) {
    MPI_Barrier (dist::comm());
  }
}

// Destructors
gdata::~gdata ()
{
  DECLARE_CCTK_PARAMETERS;
  if (barriers) {
    MPI_Barrier (dist::comm());
  }
}



// Processor management
void gdata::change_processor (comm_state& state,
                              const int newproc,
                              void* const mem)
{
  // if this function is being called with collective commbuffers turned on,
  // mimic the old state transitions here
  switch (state.thestate) {
  case state_post:
  case state_get_buffer_sizes:
    change_processor_recv (state, newproc, mem);
    change_processor_send (state, newproc, mem);
    break;
  case state_wait:
  case state_fill_send_buffers:
    change_processor_wait (state, newproc, mem);
    break;
  case state_empty_recv_buffers:
    break;
  default:
    assert(0 and "invalid state");
  }
}



// Data manipulators
void gdata::copy_from (comm_state& state,
                       const gdata* src, const ibbox& box)
{
  assert (has_storage() and src->has_storage());
  assert (all(box.lower()>=extent().lower()
          and box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
          and box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()
          and box.stride()==src->extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
          and (box.lower()-src->extent().lower())%box.stride() == 0));

  if (box.empty()) return;
  if (dist::rank() != proc() and dist::rank() != src->proc()) return;

  switch (state.thestate) {
  case state_post:
    if (proc() == src->proc()) {
      copy_from_innerloop (src, box);
    } else {
      copy_from_post (state, src, box);
    }
    break;

  case state_wait:
    if (proc() != src->proc()) {
      copy_from_wait (state, src, box);
    }
    break;

  case state_get_buffer_sizes:
    // don't count processor-local copies
    if (proc() != src->proc()) {
      // if this is a destination processor: advance its recv buffer size
      vector<comm_state::procbufdesc>& procbufs =
        state.typebufs.at(c_datatype()).procbufs;
      if (proc() == dist::rank()) {
        procbufs.at(src->proc()).recvbufsize += box.size();
        state.typebufs.at(c_datatype()).in_use = true;
      }
      // if this is a source processor: increment its send buffer size
      if (src->proc() == dist::rank()) {
        procbufs.at(proc()).sendbufsize += box.size();
        state.typebufs.at(c_datatype()).in_use = true;
      }
    }
    break;

  case state_fill_send_buffers:
    // if this is a source processor: copy its data into the send buffer
    // (the processor-local case is also handled here)
    if (src->proc() == dist::rank()) {
      if (proc() == src->proc()) {
        copy_from_innerloop (src, box);
      } else {
        copy_into_sendbuffer (state, src, box);
      }
    }
    break;

  case state_empty_recv_buffers:
    // if this is a destination processor and data has already been received
    // from the source processor: copy it from the recv buffer
    if (proc() == dist::rank() and
        state.recvbuffers_ready.at(dist::size()*c_datatype() + src->proc())) {
      copy_from_recvbuffer (state, src, box);
    }
    break;

  default:
    assert(0 and "invalid state");
  }
}


void gdata::copy_from_post (comm_state& state,
                            const gdata* src, const ibbox& box)
{
  wtime_copyfrom_recv.start();

  if (dist::rank() == proc()) {

    // this processor receives data

    wtime_copyfrom_recvinner_allocate.start();
    comm_state::gcommbuf * b = make_typed_commbuf (box);
    wtime_copyfrom_recvinner_allocate.stop();

    wtime_copyfrom_recvinner_recv.start();
    MPI_Irecv (b->pointer(), b->size(), b->datatype(), src->proc(),
               tag, dist::comm(), &b->request);
    wtime_copyfrom_recvinner_recv.stop();
    state.requests.push_back (b->request);
    state.recvbufs.push (b);

  } else {
    // this processor sends data

    wtime_copyfrom_sendinner_allocate.start();
    comm_state::gcommbuf * b = src->make_typed_commbuf (box);
    wtime_copyfrom_sendinner_allocate.stop();

    // copy data into send buffer
    wtime_copyfrom_sendinner_copy.start();
    const ibbox& ext = src->extent();
    ivect myshape = ext.shape() / ext.stride();
    ivect items = (box.upper() - box.lower()) / box.stride() + 1;
    ivect offs  = (box.lower() - ext.lower()) / ext.stride();
    char* send_buffer = (char*) b->pointer();
    int& datatypesize = state.typebufs.at(c_datatype()).datatypesize;

    for (int k = 0; k < items[2]; k++) {
      for (int j = 0; j < items[1]; j++) {
        int i = offs[0] + myshape[0]*((j+offs[1]) + myshape[1]*(k+offs[2]));
        memcpy (send_buffer, ((char*) src->storage()) + datatypesize*i,
                datatypesize * items[0]);
        send_buffer += datatypesize * items[0];
      }
    }
    wtime_copyfrom_sendinner_copy.stop();

    wtime_copyfrom_sendinner_send.start();
    MPI_Isend (b->pointer(), b->size(), b->datatype(), proc(),
               tag, dist::comm(), &b->request);
    wtime_copyfrom_sendinner_send.stop();
    state.requests.push_back (b->request);
    state.sendbufs.push (b);
  }

  wtime_copyfrom_recv.stop();
}


void gdata::copy_from_wait (comm_state& state,
                            const gdata* src, const ibbox& box)
{
  wtime_copyfrom_wait.start();

  wtime_copyfrom_recvwaitinner_wait.start();
  if (not state.requests.empty()) {
    // wait for all requests at once
    MPI_Waitall (state.requests.size(), &state.requests.front(),
                 MPI_STATUSES_IGNORE);
    state.requests.clear();
  }
  wtime_copyfrom_recvwaitinner_wait.stop();

  queue<comm_state::gcommbuf*>* const bufs =
    dist::rank() == proc() ? &state.recvbufs : &state.sendbufs;
  comm_state::gcommbuf* b = bufs->front();

  // copy data out of receive buffer
  if (bufs == &state.recvbufs) {
    wtime_copyfrom_recvwaitinner_copy.start();
    const ibbox& ext = extent();
    ivect myshape = ext.shape() / ext.stride();
    ivect items = (box.upper() - box.lower()) / box.stride() + 1;
    ivect offs  = (box.lower() - ext.lower()) / ext.stride();
    const char* recv_buffer = (const char*) b->pointer();
    int& datatypesize = state.typebufs.at(c_datatype()).datatypesize;

    for (int k = 0; k < items[2]; k++) {
      for (int j = 0; j < items[1]; j++) {
        int i = offs[0] + myshape[0]*((j+offs[1]) + myshape[1]*(k+offs[2]));
        memcpy (((char*) storage()) + datatypesize*i, recv_buffer,
                datatypesize * items[0]);
        recv_buffer += datatypesize * items[0];
      }
    }
    wtime_copyfrom_recvwaitinner_copy.stop();
  }

  wtime_copyfrom_recvwaitinner_delete.start();
  bufs->pop();
  delete b;
  wtime_copyfrom_recvwaitinner_delete.stop();

  wtime_copyfrom_wait.stop();
}


// Copy processor-local source data into communication send buffer
// of the corresponding destination processor
void gdata::copy_into_sendbuffer (comm_state& state,
                                  const gdata* src, const ibbox& box)
{
  if (proc() == src->proc()) {
    // copy on same processor
    copy_from_innerloop (src, box);
  } else {
    // copy to remote processor
    assert (src->_has_storage);
    int datatypesize = state.typebufs.at(c_datatype()).datatypesize;
    comm_state::procbufdesc& procbuf =
      state.typebufs.at(c_datatype()).procbufs.at(proc());
    assert (procbuf.sendbuf - procbuf.sendbufbase <=
            ((int)procbuf.sendbufsize - box.size()) * datatypesize);
    int const fillstate = procbuf.sendbuf + (int)box.size()*datatypesize -
                          procbuf.sendbufbase;
    assert (fillstate <= (int)procbuf.sendbufsize * datatypesize);

    // copy this processor's data into the send buffer
    const ibbox& ext = src->extent();
    ivect myshape = ext.shape() / ext.stride();
    ivect items = (box.upper() - box.lower()) / box.stride() + 1;
    ivect offs  = (box.lower() - ext.lower()) / ext.stride();
  
    assert (dim == 3);
    for (int k = 0; k < items[2]; k++) {
      for (int j = 0; j < items[1]; j++) {
        int i = offs[0] + myshape[0]*((j+offs[1]) + myshape[1]*(k+offs[2]));
        memcpy (procbuf.sendbuf,
                ((const char*) src->storage()) + datatypesize*i,
                datatypesize * items[0]);
        procbuf.sendbuf += datatypesize * items[0];
      }
    }
  
    // post the send if the buffer is full
    if (fillstate == (int)procbuf.sendbufsize * datatypesize) {
      wtime_commstate_isend.start();
      MPI_Isend (procbuf.sendbufbase, procbuf.sendbufsize,
                 state.typebufs.at(c_datatype()).mpi_datatype,
                 proc(), c_datatype(), dist::comm(),
                 &state.srequests.at(dist::size()*c_datatype() + proc()));
    }
  }
}


// Copy processor-local destination data from communication recv buffer
// of the corresponding source processor
void gdata::copy_from_recvbuffer (comm_state& state,
                                  const gdata* src, const ibbox& box)
{
  int& datatypesize = state.typebufs.at(c_datatype()).datatypesize;
  comm_state::procbufdesc& procbuf =
    state.typebufs.at(c_datatype()).procbufs.at(src->proc());
  assert (procbuf.recvbuf - procbuf.recvbufbase <=
          ((int)procbuf.recvbufsize-box.size()) * datatypesize);

  // copy this processor's data from the recv buffer
  const ibbox& ext = extent();
  ivect myshape = ext.shape() / ext.stride();
  ivect items = (box.upper() - box.lower()) / box.stride() + 1;
  ivect offs  = (box.lower() - ext.lower()) / ext.stride();

  assert (dim == 3);
  for (int k = 0; k < items[2]; k++) {
    for (int j = 0; j < items[1]; j++) {
      int i = offs[0] + myshape[0]*((j+offs[1]) + myshape[1]*(k+offs[2]));
      memcpy (((char*) storage()) + datatypesize*i,
              procbuf.recvbuf, datatypesize * items[0]);
      procbuf.recvbuf += datatypesize * items[0];
    }
  }
}


void gdata
::interpolate_from (comm_state& state,
                    const vector<const gdata*> srcs,
                    const vector<CCTK_REAL> times,
                    const ibbox& box, const CCTK_REAL time,
                    const int order_space,
                    const int order_time)
{
  assert (transport_operator != op_error);
  if (transport_operator == op_none) return;

  assert (has_storage());
  assert (all(box.lower()>=extent().lower()));
  assert (all(box.upper()<=extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  assert (srcs.size() == times.size() and srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs.at(t)->has_storage());
    assert (all(box.lower()>=srcs.at(t)->extent().lower()));
    assert (all(box.upper()<=srcs.at(t)->extent().upper()));
  }
  assert (not box.empty());
  const gdata* src = srcs.at(0);
  if (dist::rank() != proc() and dist::rank() != src->proc()) return;

  switch (state.thestate) {
  case state_post:
    if (proc() == src->proc()) {
      interpolate_from_innerloop(srcs, times, box, time,
                                 order_space, order_time);
    } else {
      interpolate_from_post(state, srcs, times, box, time,
                            order_space, order_time);
    }
    break;
  case state_wait:
    if (proc() != src->proc()) {
      copy_from_wait (state, src, box);
    }
    break;
  case state_get_buffer_sizes:
    // don't count processor-local interpolations
    if (proc() != src->proc()) {
      // if this is a destination processor: increment its recv buffer size
      vector<comm_state::procbufdesc>& procbufs =
        state.typebufs.at(c_datatype()).procbufs;
      if (proc() == dist::rank()) {
        procbufs.at(src->proc()).recvbufsize += box.size();
        state.typebufs.at(c_datatype()).in_use = true;
      }
      // if this is a source processor: increment its send buffer size
      if (src->proc() == dist::rank()) {
        procbufs.at(proc()).sendbufsize += box.size();
        state.typebufs.at(c_datatype()).in_use = true;
      }
    }
    break;
  case state_fill_send_buffers:
    // if this is a source processor: interpolate its data into the send buffer
    // (the processor-local case is also handled here)
    if (src->proc() == dist::rank()) {
      if (proc() == src->proc()) {
        interpolate_from_innerloop(srcs, times, box, time,
                                   order_space, order_time);
      } else {
        interpolate_into_sendbuffer(state, srcs, times, box,
                                    time, order_space, order_time);
      }
    }
    break;
  case state_empty_recv_buffers:
    // if this is a destination processor and data has already been received
    // from the source processor: copy it from the recv buffer
    // (the processor-local case is not handled here)
    if (proc() == dist::rank() and
        state.recvbuffers_ready.at(dist::size()*c_datatype() + src->proc())) {
      copy_from_recvbuffer(state, src, box);
    }
    break;
  default:
    assert(0 and "invalid state");
  }
}


void gdata
::interpolate_from_post (comm_state& state,
                         const vector<const gdata*> srcs,
                         const vector<CCTK_REAL> times,
                         const ibbox& box,
                         const CCTK_REAL time,
                         const int order_space,
                         const int order_time)
{
  const gdata* src = srcs.at(0);
  if (dist::rank() == proc()) {
    // interpolate from other processor

    // this processor receives data

    comm_state::gcommbuf * b = make_typed_commbuf (box);

    MPI_Irecv (b->pointer(), b->size(), b->datatype(), src->proc(),
               tag, dist::comm(), &b->request);
    state.requests.push_back (b->request);
    state.recvbufs.push (b);
  } else {
    // this processor sends data

    comm_state::gcommbuf * b = src->make_typed_commbuf (box);

    gdata * tmp = src->make_typed (varindex, transport_operator, tag);
    tmp->allocate (box, src->proc(), b->pointer());
    tmp->interpolate_from_innerloop (srcs, times, box, time,
                                     order_space, order_time);
    delete tmp;

    MPI_Isend (b->pointer(), b->size(), b->datatype(), proc(),
               tag, dist::comm(), &b->request);
    state.requests.push_back (b->request);
    state.sendbufs.push (b);
  }
}


// Interpolate processor-local source data into communication send buffer
// of the corresponding destination processor
void gdata
::interpolate_into_sendbuffer (comm_state& state,
                               const vector<const gdata*> srcs,
                               const vector<CCTK_REAL> times,
                               const ibbox& box,
                               const CCTK_REAL time,
                               const int order_space,
                               const int order_time)
{
  if (proc() == srcs.at(0)->proc()) {
    // interpolate on same processor
    interpolate_from_innerloop (srcs, times, box, time,
                                order_space, order_time);
  } else {
    // interpolate to remote processor
    const gdata* src = srcs.at(0);
    assert (src->_has_storage);
    int& datatypesize = state.typebufs.at(c_datatype()).datatypesize;
    comm_state::procbufdesc& procbuf =
      state.typebufs.at(c_datatype()).procbufs.at(proc());
    assert (procbuf.sendbuf - procbuf.sendbufbase <=
            ((int)procbuf.sendbufsize - box.size()) * datatypesize);
    assert (src->has_storage());
    int const fillstate = (procbuf.sendbuf + box.size()*datatypesize) -
                          procbuf.sendbufbase;
    assert (fillstate <= (int)procbuf.sendbufsize * datatypesize);

    // interpolate this processor's data into the send buffer
    gdata* tmp = src->make_typed (varindex, transport_operator, tag);
    tmp->allocate (box, src->proc(), procbuf.sendbuf);
    tmp->interpolate_from_innerloop (srcs, times, box, time,
                                     order_space, order_time);
    delete tmp;
  
    // advance send buffer to point to the next ibbox slot
    procbuf.sendbuf += datatypesize * box.size();
  
    // post the send if the buffer is full
    if (fillstate == procbuf.sendbufsize*datatypesize) {
      MPI_Isend (procbuf.sendbufbase, procbuf.sendbufsize,
                 state.typebufs.at(c_datatype()).mpi_datatype,
                 proc(), c_datatype(), dist::comm(),
                 &state.srequests.at(dist::size()*c_datatype() + proc()));
    }
  }
}
