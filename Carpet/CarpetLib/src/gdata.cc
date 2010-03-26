#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "bbox.hh"
#include "commstate.hh"
#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"
#include "vect.hh"

#include "gdata.hh"

using namespace std;
using namespace CarpetLib;



// Hand out the next MPI tag
static int nexttag ()
{
  DECLARE_CCTK_PARAMETERS;

  int const min_tag = 100;
  static int last = 0;
  ++last;
  if (last >= 30000) last = 0;
  return min_tag + last;
}



// Constructors
gdata::gdata (const int varindex_,
              const centering cent_,
              const operator_type transport_operator_,
              const int tag_)
  : _storage(NULL),
    varindex(varindex_),
    cent(cent_),
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



// Data manipulators

void
gdata::
copy_from (comm_state & state,
           gdata const * const src,
           ibbox const & box)
{
  vector <gdata const *> srcs (1, src);
  CCTK_REAL const time = 0.0;
  vector <CCTK_REAL> times (1, time);
  int const order_space = cent == vertex_centered ? 1 : 0;
  int const order_time = 0;
  transfer_from (state,
                 srcs, times,
                 box, box,
                 time, order_space, order_time);
}



void
gdata::
transfer_from (comm_state & state,
               vector<gdata const *> const & srcs,
               vector<CCTK_REAL>     const & times,
               ibbox const & dstbox,
               ibbox const & srcbox,
               CCTK_REAL const time,
               int const order_space,
               int const order_time)
{
  assert (has_storage());
  assert (not dstbox.empty());
  assert (all(dstbox.lower() >= extent().lower()));
  assert (all(dstbox.upper() <= extent().upper()));
  assert (all(dstbox.stride() == extent().stride()));
  assert (all((dstbox.lower() - extent().lower()) % dstbox.stride() == 0));
  
  assert (not srcbox.empty());
  assert (srcs.size() == times.size() and srcs.size() > 0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs.AT(t)->has_storage());
    assert (all(srcbox.lower() >= srcs.AT(t)->extent().lower()));
    assert (all(srcbox.upper() <= srcs.AT(t)->extent().upper()));
  }
  gdata const * const src = srcs.AT(0);
  
  assert (transport_operator != op_error);
  if (transport_operator == op_none) return;
  
  // Return early if this communication does not concern us
  if (dist::rank() != proc() and dist::rank() != src->proc()) return;
  
  // Interpolate either on the source or on the destination processor,
  // depending on whether this increases or reduces the amount of data
  int timelevel0, ntimelevels;
  find_source_timelevel (times, time, order_time, timelevel0, ntimelevels);
  assert (int (srcs.size()) >= ntimelevels);
  int const dstpoints = dstbox.size();
  int const srcpoints = srcbox.size() * ntimelevels;
  bool const interp_on_src = dstpoints <= srcpoints;
  int const npoints = interp_on_src ? dstpoints : srcpoints;
  
  switch (state.thestate) {
    
  case state_get_buffer_sizes:
    // don't count processor-local copies
    if (proc() != src->proc()) {
      // if this is a destination processor: advance its recv buffer
      // size
      if (proc() == dist::rank()) {
        state.reserve_recv_space (c_datatype(), src->proc(), npoints);
      }
      // if this is a source processor: increment its send buffer size
      if (src->proc() == dist::rank()) {
        state.reserve_send_space (c_datatype(), proc(), npoints);
      }
    }
    break;
    
  case state_fill_send_buffers:
    // if this is a source processor: copy its data into the send
    // buffer
    if (proc() != src->proc()) {
      if (src->proc() == dist::rank()) {
        if (interp_on_src) {
          size_t const sendbufsize = c_datatype_size() * dstbox.size();
          void * const sendbuf =
            state.send_buffer (c_datatype(), proc(), dstbox.size());
          gdata * const buf =
            make_typed (varindex, cent, transport_operator, tag);
          buf->allocate (dstbox, src->proc(), sendbuf, sendbufsize);
          buf->transfer_from_innerloop
            (srcs, times, dstbox, time, order_space, order_time);
          delete buf;
          state.commit_send_space (c_datatype(), proc(), dstbox.size());
        } else {
          for (int tl = timelevel0; tl < timelevel0 + ntimelevels; ++ tl) {
            size_t const sendbufsize = c_datatype_size() * srcbox.size();
            void * const sendbuf =
              state.send_buffer (c_datatype(), proc(), srcbox.size());
            gdata * const buf =
              make_typed (varindex, cent, transport_operator, tag);
            buf->allocate (srcbox, src->proc(), sendbuf, sendbufsize);
            buf->copy_from_innerloop (srcs.AT(tl), srcbox);
            delete buf;
            state.commit_send_space (c_datatype(), proc(), srcbox.size());
          }
        }
      }
    }
    break;
    
  case state_do_some_work:
    // handle the processor-local case
    if (proc() == src->proc()) {
      if (proc() == dist::rank()) {
        transfer_from_innerloop
          (srcs, times, dstbox, time, order_space, order_time);
      }
    }
    break;
    
  case state_empty_recv_buffers:
    // if this is a destination processor: copy it from the recv
    // buffer
    if (proc() != src->proc()) {
      if (proc() == dist::rank()) {
        if (interp_on_src) {
          size_t const recvbufsize = c_datatype_size() * dstbox.size();
          void * const recvbuf =
            state.recv_buffer (c_datatype(), src->proc(), dstbox.size());
          gdata * const buf =
            make_typed (varindex, cent, transport_operator, tag);
          buf->allocate (dstbox, proc(), recvbuf, recvbufsize);
          state.commit_recv_space (c_datatype(), src->proc(), dstbox.size());
          copy_from_innerloop (buf, dstbox);
          delete buf;
        } else {
          gdata const * const null = 0;
          vector <gdata const *> bufs (timelevel0 + ntimelevels, null);
          for (int tl = timelevel0; tl < timelevel0 + ntimelevels; ++ tl) {
            size_t const recvbufsize = c_datatype_size() * srcbox.size();
            void * const recvbuf =
              state.recv_buffer (c_datatype(), src->proc(), srcbox.size());
            gdata * const buf =
              make_typed (varindex, cent, transport_operator, tag);
            buf->allocate (srcbox, proc(), recvbuf, recvbufsize);
            state.commit_recv_space (c_datatype(), src->proc(), srcbox.size());
            bufs.AT(tl) = buf;
          }
          transfer_from_innerloop
            (bufs, times, dstbox, time, order_space, order_time);
          for (int tl = timelevel0; tl < timelevel0 + ntimelevels; ++ tl) {
            delete bufs.AT(tl);
          }
        }
      }
    }
    break;
    
  default:
    assert (0);
  }
}



void
gdata::
find_source_timelevel (vector <CCTK_REAL> const & times,
                       CCTK_REAL const time,
                       int const order_time,
                       int & timelevel0,
                       int & ntimelevels)
  const
{
  // Ensure that the times are consistent
  assert (times.size() > 0);
  assert (order_time >= 0);
  
  CCTK_REAL const eps = 1.0e-12;
  CCTK_REAL const min_time = * min_element (times.begin(), times.end());
  CCTK_REAL const max_time = * max_element (times.begin(), times.end());
  CCTK_REAL const some_time = abs (min_time) + abs (max_time);
  if (transport_operator != op_copy) {
    if (time < min_time - eps * some_time or
        time > max_time + eps * some_time)
    {
      ostringstream buf;
      buf << setprecision (17)
          << "Internal error: extrapolation in time."
          << "  time=" << time
          << "  times=" << times;
      CCTK_WARN (0, buf.str().c_str());
    }
  }
  
  // Use this timelevel, or interpolate in time if set to -1
  int timelevel = -1;
  
  // Try to avoid time interpolation if possible
  if (timelevel == -1) {
    if (times.size() == 1) {
      timelevel = 0;
    }
  }
  if (timelevel == -1) {
    if (transport_operator == op_copy) {
      timelevel = 0;
    }
  }
  if (timelevel == -1) {
    for (size_t tl=0; tl<times.size(); ++tl) {
      static_assert (abs(0.1) > 0,
                     "Function CarpetLib::abs has wrong signature");
      if (abs (times.AT(tl) - time) < eps * some_time) {
        timelevel = tl;
        break;
      }
    }
  }
  
  if (timelevel >= 0) {
    timelevel0 = timelevel;
    ntimelevels = 1;
  } else {
    timelevel0 = 0;
    ntimelevels = order_time + 1;
  }
  
  assert (timelevel0 >= 0 and timelevel0 < (int)times.size());
  assert (ntimelevels > 0);
}
