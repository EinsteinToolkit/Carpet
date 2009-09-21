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



list<gdata*> gdata::allgdata;



// Constructors
gdata::gdata (const int varindex_,
              const centering cent_,
              const operator_type transport_operator_)
  : _storage(NULL),
    varindex(varindex_),
    cent(cent_),
    transport_operator(transport_operator_),
    _has_storage(false),
    comm_active(false)
{
  DECLARE_CCTK_PARAMETERS;
  
  allgdatai = allgdata.insert(allgdata.end(), this);
  
  if (barriers) {
    MPI_Barrier (dist::comm());
  }
}

// Destructors
gdata::~gdata ()
{
  DECLARE_CCTK_PARAMETERS;
  
  allgdata.erase(allgdatai);
  
  if (barriers) {
    MPI_Barrier (dist::comm());
  }
}



// Storage management

ivect
gdata::
allocated_memory_shape (ibbox const& extent)
{
  DECLARE_CCTK_PARAMETERS;
  ivect shape = max (ivect(0), extent.shape() / extent.stride());
  // Enlarge shape to avoid multiples of cache line colours
  if (avoid_arraysize_bytes > 0) {
    for (int d=0; d<dim; ++d) {
      if (shape[d] > 0 and
          shape[d] * sizeof(CCTK_REAL) % avoid_arraysize_bytes == 0)
      {
        ++shape[d];
      }
    }
  }
  return shape;
}



// Data manipulators

void
gdata::
copy_from (comm_state & state,
           gdata const * const src,
           ibbox const & box,
           int const dstproc,
           int const srcproc)
{
  vector <gdata const *> const srcs (1, src);
  CCTK_REAL const time = 0.0;
  vector <CCTK_REAL> const times (1, time);
  transfer_from (state,
                 srcs, times,
                 box, box,
                 dstproc, srcproc,
                 time, 0, 0);
}



void
gdata::
transfer_from (comm_state & state,
               vector<gdata const *> const & srcs,
               vector<CCTK_REAL>     const & times,
               ibbox const & dstbox,
               ibbox const & srcbox,
               int const dstproc,
               int const srcproc,
               CCTK_REAL const time,
               int const order_space,
               int const order_time)
{
  bool const is_dst = dist::rank() == dstproc;
  bool const is_src = dist::rank() == srcproc;
  // Return early if this communication does not concern us
  assert (is_dst or is_src);    // why should we be here?
  if (not is_dst and not is_src) return;

  if (is_dst) {
    assert (proc() == dstproc);
    assert (has_storage());
    assert (not dstbox.empty());
    assert (all(dstbox.lower() >= extent().lower()));
    assert (all(dstbox.upper() <= extent().upper()));
    assert (all(dstbox.stride() == extent().stride()));
    assert (all((dstbox.lower() - extent().lower()) % dstbox.stride() == 0));
  }
  
  if (is_src) {
    assert (not srcbox.empty());
    assert (srcs.size() == times.size() and srcs.size() > 0);
    for (int t=0; t<(int)srcs.size(); ++t) {
      assert (srcs.AT(t)->proc() == srcproc);
      assert (srcs.AT(t)->has_storage());
      assert (all(srcbox.lower() >= srcs.AT(t)->extent().lower()));
      assert (all(srcbox.upper() <= srcs.AT(t)->extent().upper()));
    }
  }
  gdata const * const src = is_src ? srcs.AT(0) : NULL;
  
  operator_type const my_transport_operator =
    is_dst ? transport_operator : src->transport_operator;
  assert (my_transport_operator != op_error);
  assert (my_transport_operator != op_none); // why should we be here?
  if (my_transport_operator == op_none) return;
  
  // Interpolate either on the source or on the destination processor,
  // depending on whether this increases or reduces the amount of data
  int timelevel0, ntimelevels;
  find_source_timelevel
    (times, time, order_time, my_transport_operator, timelevel0, ntimelevels);
  if (is_src) assert (int (srcs.size()) >= ntimelevels);
  int const dstpoints = dstbox.size();
  int const srcpoints = srcbox.size() * ntimelevels;
  bool const interp_on_src = dstpoints <= srcpoints;
  int const npoints = interp_on_src ? dstpoints : srcpoints;
  
  switch (state.thestate) {
    
  case state_get_buffer_sizes:
    // don't count processor-local copies
    if (not (is_dst and is_src)) {
      if (is_dst) {
        // increment the recv buffer size
        state.reserve_recv_space (c_datatype(), srcproc, npoints);
      }
      if (is_src) {
        // increment the send buffer size
        state.reserve_send_space (src->c_datatype(), dstproc, npoints);
      }
    }
    break;
    
  case state_fill_send_buffers:
    if (not (is_dst and is_src)) {
      if (is_src) {
        // copy the data into the send buffer
        if (interp_on_src) {
          size_t const sendbufsize = src->c_datatype_size() * dstbox.size();
          void * const sendbuf =
            state.send_buffer (src->c_datatype(), dstproc, dstbox.size());
          gdata * const buf =
            src->make_typed (src->varindex, src->cent, src->transport_operator);
          buf->allocate (dstbox, srcproc, sendbuf, sendbufsize);
          buf->transfer_from_innerloop
            (srcs, times, dstbox, time, order_space, order_time);
          delete buf;
          state.commit_send_space (src->c_datatype(), dstproc, dstbox.size());
        } else {
          for (int tl = timelevel0; tl < timelevel0 + ntimelevels; ++ tl) {
            size_t const sendbufsize = src->c_datatype_size() * srcbox.size();
            void * const sendbuf =
              state.send_buffer (src->c_datatype(), dstproc, srcbox.size());
            gdata * const buf =
              src->make_typed (src->varindex, src->cent,
                               src->transport_operator);
            buf->allocate (srcbox, srcproc, sendbuf, sendbufsize);
            buf->copy_from_innerloop (srcs.AT(tl), srcbox);
            delete buf;
            state.commit_send_space (src->c_datatype(), dstproc, srcbox.size());
          }
        }
      }
    }
    break;
    
  case state_do_some_work:
    // handle the processor-local case
    if (is_dst and is_src) {
      transfer_from_innerloop
        (srcs, times, dstbox, time, order_space, order_time);
    }
    break;
    
  case state_empty_recv_buffers:
    if (not (is_dst and is_src)) {
      if (is_dst) {
        // copy from the recv buffer
        if (interp_on_src) {
          size_t const recvbufsize = c_datatype_size() * dstbox.size();
          void * const recvbuf =
            state.recv_buffer (c_datatype(), srcproc, dstbox.size());
          gdata * const buf = make_typed (varindex, cent, transport_operator);
          buf->allocate (dstbox, dstproc, recvbuf, recvbufsize);
          state.commit_recv_space (c_datatype(), srcproc, dstbox.size());
          copy_from_innerloop (buf, dstbox);
          delete buf;
        } else {
          gdata const * const null = NULL;
          vector <gdata const *> bufs (ntimelevels, null);
          vector <CCTK_REAL> timebuf (ntimelevels);
          for (int tl = 0; tl < ntimelevels; ++ tl) {
            size_t const recvbufsize = c_datatype_size() * srcbox.size();
            void * const recvbuf =
              state.recv_buffer (c_datatype(), srcproc, srcbox.size());
            gdata * const buf = make_typed (varindex, cent, transport_operator);
            buf->allocate (srcbox, dstproc, recvbuf, recvbufsize);
            state.commit_recv_space (c_datatype(), srcproc, srcbox.size());
            bufs.AT(tl) = buf;
            timebuf.AT(tl) = times.AT(timelevel0 + tl);
          }
          transfer_from_innerloop
            (bufs, timebuf, dstbox, time, order_space, order_time);
          for (int tl = 0; tl < ntimelevels; ++ tl) {
            delete bufs.AT(tl);
          }
        }
      }
    }
    break;
    
  default:
    assert (0); abort();
  }
}



void
gdata::
find_source_timelevel (vector <CCTK_REAL> const & times,
                       CCTK_REAL const time,
                       int const order_time,
                       operator_type const transport_operator,
                       int & timelevel0,
                       int & ntimelevels)
{
  // Ensure that the times are consistent
  assert (times.size() > 0);
  assert (order_time >= 0);
  
  CCTK_REAL const eps = 1.0e-12;
  CCTK_REAL const min_time = * min_element (times.begin(), times.end());
  CCTK_REAL const max_time = * max_element (times.begin(), times.end());
  // TODO: Use a real delta-time from somewhere instead of 1.0
  CCTK_REAL const some_time = abs (min_time) + abs (max_time) + 1.0;
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



size_t
gdata::
allmemory ()
{
  size_t mem = memoryof(allgdata);
  for (list<gdata*>::const_iterator
         gdatai = allgdata.begin(); gdatai != allgdata.end(); ++ gdatai)
  {
    mem += memoryof(**gdatai);
  }
  return mem;
}
