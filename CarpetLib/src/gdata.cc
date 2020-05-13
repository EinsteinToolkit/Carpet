#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <vectors.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <typeinfo>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include "bbox.hh"
#include "cacheinfo.hh"
#include "commstate.hh"
#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"
#include "vect.hh"

#include "gdata.hh"

namespace CarpetLib {
using namespace std;

template <typename T, int D>
ostream &operator<<(ostream &os, slab<T, D> const &slabinfo) {
  return os << "slab<" << typeid(T).name() << "," << D << ">"
            << "{offset=" << slabinfo.offset << "}";
}

template ostream &operator<<(ostream &os, slab<int, dim> const &slabinfo);
template ostream &operator<<(ostream &os, slab<CCTK_REAL, dim> const &slabinfo);

namespace dist {
template <typename T, int D> MPI_Datatype mpi_datatype(slab<T, D> const &) {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    slab<T, D> const s;
    typedef vect<T, D> avect;
#define ENTRY(type, name)                                                      \
  {                                                                            \
      sizeof s.name / sizeof(type), /* count elements */                       \
      (const char *)&s.name -                                                  \
          (const char *)&s,       /* offsetof doesn't work (why?) */           \
      dist::mpi_datatype<type>(), /* find MPI datatype */                      \
      STRINGIFY(name),            /* field name */                             \
      STRINGIFY(type),            /* type name */                              \
  }
    dist::mpi_struct_descr_t const descr[] = {ENTRY(avect, offset)};
#undef ENTRY
    ostringstream buf;
    buf << "slab<" << typeid(T).name() << "," << D << ">";
    newtype = dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                        buf.str().c_str(), sizeof s);
    initialised = true;
  }
  return newtype;
}

template MPI_Datatype mpi_datatype(slab<int, dim> const &);
}

set<gdata *> gdata::allgdata;

// Constructors
gdata::gdata(const int varindex_, const centering cent_,
             const operator_type transport_operator_)
    : _storage(NULL), varindex(varindex_), cent(cent_),
      transport_operator(transport_operator_), _has_storage(false),
      comm_active(false) {
  DECLARE_CCTK_PARAMETERS;

  allgdata.insert(this);

  if (barriers) {
    dist::barrier(dist::comm(), 783988953, "CarpetLib::gdata::gdata");
  }
}

// Destructors
gdata::~gdata() {
  DECLARE_CCTK_PARAMETERS;

  allgdata.erase(this);

  if (barriers) {
    dist::barrier(dist::comm(), 109687805, "CarpetLib::gdata::~gdata");
  }
}

bool gdata::fence_is_energized() {
  DECLARE_CCTK_PARAMETERS;

  return electric_fence;
}

// Data manipulators

void gdata::copy_data(gdata * const dst, comm_state &state, gdata const *const src,
                      ibbox const &dstbox, ibbox const &srcbox,
                      islab const *restrict const slabinfo, int const dstproc,
                      int const srcproc) {
  assert(dst or src); // why should we be here?
  vector<gdata const *> const srcs(1, src);
  CCTK_REAL const time = 0.0;
  vector<CCTK_REAL> const times(1, time);
  int const order_space = (dst ? dst->cent : src->cent) == vertex_centered ? 1 : 0;
  int const order_time = 0;
  transfer_data(dst, state, srcs, times, dstbox, srcbox, slabinfo, dstproc, srcproc,
                time, order_space, order_time);
}

void gdata::transfer_data(gdata * const dst, comm_state &state, vector<gdata const *> const &srcs,
                          vector<CCTK_REAL> const &times, ibbox const &dstbox,
                          ibbox const &srcbox,
                          islab const *restrict const slabinfo,
                          int const dstproc, int const srcproc,
                          CCTK_REAL const time, int const order_space,
                          int const order_time) {
  bool const is_dst = dist::rank() == dstproc;
  bool const is_src = dist::rank() == srcproc;
  assert(not is_dst or dst);
  // Return early if this communication does not concern us
  assert(is_dst or is_src); // why should we be here?
  if (not is_dst and not is_src)
    return;

  if (is_dst) {
    assert(dst->proc() == dstproc);
    assert(dst->has_storage());
    assert(not dstbox.empty());
    assert(all(dstbox.lower() >= dst->extent().lower()));
    assert(all(dstbox.upper() <= dst->extent().upper()));
    assert(all(dstbox.stride() == dst->extent().stride()));
    assert(all((dstbox.lower() - dst->extent().lower()) % dstbox.stride() == 0));
  }

  if (is_src) {
    assert(not srcbox.empty());
    assert(srcs.size() == times.size() and srcs.size() > 0);
    for (int t = 0; t < (int)srcs.size(); ++t) {
      assert(srcs.AT(t)->proc() == srcproc);
      assert(srcs.AT(t)->has_storage());
    }
  }
  gdata const *const src = srcs.AT(0);

  operator_type const my_transport_operator =
      is_dst ? dst->transport_operator : src->transport_operator;
  assert(my_transport_operator != op_error);
  assert(my_transport_operator != op_none); // why should we be here?
  if (my_transport_operator == op_none)
    return;

  // Interpolate either on the source or on the destination process,
  // depending on whether this increases or reduces the amount of data
  int timelevel0, ntimelevels;
  dst->find_source_timelevel(times, time, order_time, my_transport_operator,
                             timelevel0, ntimelevels);
  if (is_src)
    assert(int(srcs.size()) >= ntimelevels);
  // Communication buffers are neither padded nor aligned
  int const dstpoints = dstbox.size();
  int const srcpoints = srcbox.size() * ntimelevels;
  bool const interp_on_src = dstpoints <= srcpoints;
  int const npoints = interp_on_src ? dstpoints : srcpoints;

  switch (state.thestate) {

  case state_get_buffer_sizes:
    // don't count process-local copies
    if (not(is_dst and is_src)) {
      if (is_dst) {
        // increment the recv buffer size
        state.reserve_recv_space(dst->c_datatype(), srcproc, npoints);
      }
      if (is_src) {
        // increment the send buffer size
        state.reserve_send_space(src->c_datatype(), dstproc, npoints);
      }
    }
    break;

  case state_fill_send_buffers:
    if (not(is_dst and is_src)) {
      if (is_src) {
        // copy the data into the send buffer
        if (interp_on_src) {
          size_t const sendbufsize = src->c_datatype_size() * dstbox.size();
          void *const sendbuf =
              state.send_buffer(src->c_datatype(), dstproc, dstbox.size());
          gdata *const buf = src->make_typed(src->varindex, src->cent,
                                             src->transport_operator);
          buf->allocate(dstbox, dstbox.sizes(), ivect(0), srcproc, sendbuf,
                        sendbufsize);
          buf->transfer_from_innerloop(srcs, times, dstbox, srcbox, slabinfo,
                                       time, order_space, order_time);
          delete buf;
          state.commit_send_space(src->c_datatype(), dstproc, dstbox.size());
        } else {
          for (int tl = timelevel0; tl < timelevel0 + ntimelevels; ++tl) {
            size_t const sendbufsize = src->c_datatype_size() * srcbox.size();
            void *const sendbuf =
                state.send_buffer(src->c_datatype(), dstproc, srcbox.size());
            gdata *const buf = src->make_typed(src->varindex, src->cent,
                                               src->transport_operator);
            buf->allocate(srcbox, srcbox.sizes(), ivect(0), srcproc, sendbuf,
                          sendbufsize);
            buf->copy_from_innerloop(srcs.AT(tl), srcbox, srcbox, NULL);
            delete buf;
            state.commit_send_space(src->c_datatype(), dstproc, srcbox.size());
          }
        }
      }
    }
    break;

  case state_do_some_work:
    // handle the process-local case
    if (is_dst and is_src) {
      dst->transfer_from_innerloop(srcs, times, dstbox, srcbox, slabinfo, time,
                                   order_space, order_time);
    }
    break;

  case state_empty_recv_buffers:
    if (not(is_dst and is_src)) {
      if (is_dst) {
        // copy from the recv buffer
        if (interp_on_src) {
          size_t const recvbufsize = dst->c_datatype_size() * dstbox.size();
          void *const recvbuf =
              state.recv_buffer(dst->c_datatype(), srcproc, dstbox.size());
          gdata *const buf = dst->make_typed(dst->varindex, dst->cent,
                                             dst->transport_operator);
          buf->allocate(dstbox, dstbox.sizes(), ivect(0), dstproc, recvbuf,
                        recvbufsize);
          state.commit_recv_space(dst->c_datatype(), srcproc, dstbox.size());
          dst->copy_from_innerloop(buf, dstbox, dstbox, NULL);
          delete buf;
        } else {
          gdata const *const null = NULL;
          vector<gdata const *> bufs(ntimelevels, null);
          vector<CCTK_REAL> timebuf(ntimelevels);
          for (int tl = 0; tl < ntimelevels; ++tl) {
            size_t const recvbufsize = dst->c_datatype_size() * srcbox.size();
            void *const recvbuf =
                state.recv_buffer(dst->c_datatype(), srcproc, srcbox.size());
            gdata *const buf = dst->make_typed(dst->varindex, dst->cent,
                                               dst->transport_operator);
            buf->allocate(srcbox, srcbox.sizes(), ivect(0), dstproc, recvbuf,
                          recvbufsize);
            state.commit_recv_space(dst->c_datatype(), srcproc, srcbox.size());
            bufs.AT(tl) = buf;
            timebuf.AT(tl) = times.AT(timelevel0 + tl);
          }
          dst->transfer_from_innerloop(bufs, timebuf, dstbox, srcbox, slabinfo,
                                       time, order_space, order_time);
          for (int tl = 0; tl < ntimelevels; ++tl)
            delete bufs.AT(tl);
        }
      }
    }
    break;

  default:
    assert(0);
    abort();
  }
}

void gdata::find_source_timelevel(vector<CCTK_REAL> const &times,
                                  CCTK_REAL const time, int const order_time,
                                  operator_type const op, int &timelevel0,
                                  int &ntimelevels) const {
  // Ensure that the times are consistent
  assert(times.size() > 0);
  assert(order_time >= 0);

  CCTK_REAL const eps = 1.0e-12;
  CCTK_REAL const min_time = *min_element(times.begin(), times.end());
  CCTK_REAL const max_time = *max_element(times.begin(), times.end());
  // TODO: Use a real delta-time from somewhere instead of 1.0
  CCTK_REAL const some_time = std::fabs(min_time) + std::fabs(max_time) + 1.0;
  if (op != op_copy) {
    if (time < min_time - eps * some_time or
        time > max_time + eps * some_time) {
      ostringstream buf;
      buf << setprecision(17) << "Internal error: extrapolation in time.";
      if (varindex >= 0) {
        char *const fullname = CCTK_FullName(varindex);
        buf << "  variable=" << fullname;
        ::free(fullname);
      }
      buf << "  time=" << time << "  times=" << times;
      CCTK_ERROR(buf.str().c_str());
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
    if (op == op_copy) {
      timelevel = 0;
    }
  }
  if (timelevel == -1) {
    for (size_t tl = 0; tl < times.size(); ++tl) {
      if (std::fabs(times.AT(tl) - time) < eps * some_time) {
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

  assert(timelevel0 >= 0 and timelevel0 < (int)times.size());
  assert(ntimelevels > 0);
}

size_t gdata::allmemory() {
  size_t mem = memoryof(allgdata);
  for (set<gdata *>::const_iterator gdatai = allgdata.begin();
       gdatai != allgdata.end(); ++gdatai) {
    mem += memoryof(**gdatai);
  }
  return mem;
}
}
