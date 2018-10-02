#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ostream>
#include <sstream>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>
#include <dist.hh>

#include "interpolate.h"

using namespace std;
using namespace Carpet;
using namespace CarpetLib;

// Divide rounding downwards (towards negative infinity)
static int divdn(int const a, int const b) {
  assert(b > 0);
  // Note: C++ division rounds towards zero
  return a >= 0 ? a / b : (a - (b - 1)) / b;
}

// Round downwards towards next multiple
static int rounddn(int const a, int const b) {
  assert(b > 0);
  return divdn(a, b) * b;
}

struct xferinfo_t {
  int m;                            // map
  sendrecv_pseudoregion_t sendrecv; // regions and components
  islab slab;                       // slabbing offset
  friend ostream &operator<<(ostream &os, const xferinfo_t &xi) {
    return os << "xferinfo_t{m=" << xi.m << ",sendrecv=" << xi.sendrecv
              << ",slab=" << xi.slab << "}";
  }
};

namespace CarpetLib {
namespace dist {
MPI_Datatype mpi_datatype(xferinfo_t const &);
template <> MPI_Datatype mpi_datatype<xferinfo_t>() {
  xferinfo_t dummy;
  return mpi_datatype(dummy);
}
} // namespace dist
} // namespace CarpetLib

#include <mpi_string.hh>

namespace CarpetLib {
namespace dist {
MPI_Datatype mpi_datatype(xferinfo_t const &) {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static xferinfo_t s;
#define ENTRY(type, name)                                                      \
  {                                                                            \
      sizeof s.name / sizeof(type), /* count elements */                       \
      (char *)&s.name - (char *)&s, /* offsetof doesn't work (why?) */         \
      dist::mpi_datatype<type>(),   /* find MPI datatype */                    \
      STRINGIFY(name),              /* field name */                           \
      STRINGIFY(type),              /* type name */                            \
  }
    dist::mpi_struct_descr_t const descr[] = {
        ENTRY(int, m), ENTRY(sendrecv_pseudoregion_t, sendrecv),
        ENTRY(islab, slab)};
#undef ENTRY
    newtype = dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                        "xferinfo_t", sizeof s);
    initialised = true;
  }
  return newtype;
}
} // namespace dist
} // namespace CarpetLib

struct levelinfo_t {
  int regridding_epoch;
  vector<xferinfo_t> xferinfos;
  levelinfo_t() : regridding_epoch(-1) {}
  friend ostream &operator<<(ostream &os, const levelinfo_t &li) {
    os << "levelinfo_t{regridding_epoch=" << li.regridding_epoch
       << ",xferinfos=vector<xferinfo_t>[";
    for (const auto &xi : li.xferinfos)
      os << xi << ",";
    os << "]}";
    return os;
  }
};
ostream &operator<<(ostream &os, const vector<levelinfo_t> &lis) {
  os << "vector<levelinfo_t>[";
  for (const auto &li : lis)
    os << li << ",";
  os << "]";
  return os;
}

static vector<levelinfo_t> levelinfos;

static void periodic_carpet(cGH const *restrict const cctkGH, int const size,
                            CCTK_INT const *restrict const stencil,
                            CCTK_INT const do_periodic[3],
                            CCTK_INT const *restrict const vars,
                            int const nvars) {
  DECLARE_CCTK_PARAMETERS;

  // Check arguments
  assert(cctkGH);
  assert(stencil);
  assert(vars);

  // Early return if possible, because the code below requires at
  // least one variable
  if (nvars == 0)
    return;

  // Get and check group info
  cGroup group;
  for (int i = 0; i < nvars; ++i) {
    int const vi = vars[i];
    assert(vi >= 0 and vi < CCTK_NumVars());
    int const gi = CCTK_GroupIndexFromVarI(vi);
    assert(gi >= 0 and gi < CCTK_NumGroups());

    int ierr = CCTK_GroupData(gi, &group);
    assert(not ierr);
    assert(group.grouptype == CCTK_GF);
    assert(group.disttype == CCTK_DISTRIB_DEFAULT);
    int const vartypesize = CCTK_VarTypeSize(group.vartype);
    assert(vartypesize > 0);

    if (group.dim != size)
      CCTK_VERROR("The group \"%s\" has dimension %d, but the given stencil "
                  "has size %d",
                  CCTK_GroupNameFromVarI(vi), group.dim, size);
  }

  // Ensure we are in level mode (and not global mode)
  assert(reflevel >= 0);

  // The cached data structures, containing the communication schedule
  // for each refinement level
  int const ml = mglevel;
  int const rl = reflevel;
  int const tl = 0;

  if (reflevels >= (int)levelinfos.size())
    levelinfos.resize(reflevels);
  levelinfo_t &levelinfo = levelinfos.at(reflevel);
  vector<xferinfo_t> &xferinfos = levelinfo.xferinfos;

  // Do we need a new levelinfo? (We need a new levelinfo after each
  // regridding.)
  if (levelinfo.regridding_epoch != level_regridding_epochs.at(reflevel)) {
    if (verbose)
      CCTK_INFO("New regridding epoch -- (re-)generating levelinfos");

    levelinfo.regridding_epoch = level_regridding_epochs.at(reflevel);
    xferinfos.clear();

    // Loop over all components and find out how to fill their
    // periodic boundaries

    gh const &hh = *vhh.at(Carpet::map);
    ibbox const &domain_exterior = hh.baseextent(mglevel, reflevel);
    i2vect const &boundary_width = hh.boundary_width;
    ibbox const domain_active = domain_exterior.expand(-boundary_width);
    dh const &dd = *vdd.at(Carpet::map);
    dh::light_cboxes const &light_level =
        dd.light_boxes.AT(mglevel).AT(reflevel);

    if (hh.local_components(reflevel) != 1)
      CCTK_ERROR("Cannot handle more than one local component");

    // Interior of the domain with respect to periodic boundaries,
    // i.e. periodic boundaries cut off, but all other boundaries
    // still included
    ibbox domain_perint;
    {
      ivect lo, hi, str;
      for (int d = 0; d < dim; ++d) {
        ibbox const &dom = do_periodic[d] ? domain_active : domain_exterior;
        lo[d] = dom.lower()[d];
        hi[d] = dom.upper()[d];
        str[d] = dom.stride()[d];
      }
      domain_perint = ibbox(lo, hi, str);
    }

    // CCTK_INT width[2*dim];
    // CCTK_INT is_internal[2*dim];
    // CCTK_INT is_staggered[2*dim];
    // CCTK_INT shiftout[2*dim];
    // int ierr = GetBoundarySpecification
    //   (2*dim, width, is_internal, is_staggered, shiftout);
    // assert(not ierr);

    // // Domain size
    // CCTK_REAL physical_min[dim];
    // CCTK_REAL physical_max[dim];
    // CCTK_REAL interior_min[dim];
    // CCTK_REAL interior_max[dim];
    // CCTK_REAL exterior_min[dim];
    // CCTK_REAL exterior_max[dim];
    // CCTK_REAL spacing[dim];
    // ierr = GetDomainSpecification(dim,
    //                               physical_min, physical_max,
    //                               interior_min, interior_max,
    //                               exterior_min, exterior_max,
    //                               spacing);
    // assert(not ierr);

    // Domain size on this level, in terms of grid points
    // ivect npoints;
    // for (int d=0; d<dim; ++d) {
    //   CCTK_REAL const rnpoints =
    //     (exterior_max[d] - exterior_min[d]) / spacing[d]
    //     + 1 - (width[2*d] + width[2*d+1]);
    //   npoints[d] = floor(rnpoints + 0.01);
    //   // Ensure the domain size is integer
    //   assert(fabs(npoints[d] - rnpoints) < 0.01);
    // }

    DECLARE_CCTK_ARGUMENTS;

    // Exterior of this component
    ibbox const &ext = light_level.at(component).exterior;

    // Collect all periodic boundaries for this component
    ibset dst_bset;
    for (int d = 0; d < dim; ++d) {
      for (int f = 0; f < 2; ++f) {
        if (cctk_bbox[2 * d + f] and do_periodic[d]) {
          // This is a periodic boundary
          ivect lo = ext.lower();
          ivect up = ext.upper();
          ivect const str = ext.stride();
          int const bnd = boundary_width[f][d];
          if (f == 0)
            up[d] = lo[d] + (bnd - 1) * str[d];
          else
            lo[d] = up[d] - (bnd - 1) * str[d];
          ibbox const dst_bbox(lo, up, str);
          assert(dst_bbox.shape()[d] / dst_bbox.stride()[d] == bnd);
          dst_bset |= dst_bbox;
        }
      }
    }

    // Determine the source bboxes
    while (not dst_bset.empty()) {

      // Pick one (arbitrary) dst_bbox
      ibbox const &dst_bbox1 = *dst_bset.begin();

      // Find the slabbing offset (a multiple of the domain size)
      // that moves dst_bbox1 into the domain
      islab slab;
      for (int d = 0; d < dim; ++d) {
        if (not do_periodic[d]) {
          slab.offset[d] = 0;
        } else {
          assert((dst_bbox1.lower()[d] - domain_perint.lower()[d]) %
                     domain_perint.stride()[d] ==
                 0);
          slab.offset[d] =
              -rounddn(dst_bbox1.lower()[d] - domain_perint.lower()[d],
                       domain_perint.shape()[d]) /
              domain_perint.stride()[d];
        }
      }
      assert(not all(slab.offset == 0)); // this would be trivial

      // Determine corresponding source bbox
      ibbox const src_bbox1 = dst_bbox1.shift(slab.offset);

      // Ensure we copy from within the domain, not from boundary
      // points
      ibbox const src_bbox = src_bbox1 & domain_perint;
      assert(not src_bbox.empty());
      ibbox const dst_bbox = src_bbox.shift(-slab.offset);
      assert(dst_bbox.is_contained_in(dst_bbox1));
      dst_bset -= dst_bbox;

      // Find source components for all points
      ibset src_bset(src_bbox);
      for (int oc = 0; oc < hh.components(reflevel); ++oc) {
        // We cannot copy only from the interior (the "good"
        // points), because this excludes the ghost zones, and we
        // may need to copy from prolongation boundaries. (This
        // assumes that synchronisation and prolongation occurs
        // before periodic boundaries are applied!)
        ibset const intersection = src_bset & light_level.at(oc).exterior;
        if (not intersection.empty()) {
          src_bset -= intersection;

          // Loop over the intersection and insert respective
          // xferinfos
          ibset::const_iterator const ib = intersection.begin();
          ibset::const_iterator const ie = intersection.end();
          for (ibset::const_iterator ireg = ib; ireg != ie; ++ireg) {
            ibbox const src_reg = *ireg;
            ibbox const dst_reg = src_reg.shift(-slab.offset);
            xferinfo_t xferinfo;
            xferinfo.m = Carpet::map;
            xferinfo.sendrecv.send.extent = src_reg;
            xferinfo.sendrecv.send.component = oc;
            xferinfo.sendrecv.recv.extent = dst_reg;
            xferinfo.sendrecv.recv.component = component;
            xferinfo.slab = slab;
            xferinfos.push_back(xferinfo);
          }
        }
      } // for oc
      // Ensure we know how to handle all points
      assert(src_bset.empty());

    } // while dst_bset

    // Tell the source processes what they have to send
    {
      int const p = dist::rank();
      vector<vector<xferinfo_t> > sends(dist::size());
      vector<xferinfo_t>::const_iterator const ib = xferinfos.begin();
      vector<xferinfo_t>::const_iterator const ie = xferinfos.end();
      for (vector<xferinfo_t>::const_iterator ixferinfo = ib; ixferinfo != ie;
           ++ixferinfo) {
        xferinfo_t const &xferinfo = *ixferinfo;
        gh const &hh = *vhh.at(xferinfo.m);
        int const oc = xferinfo.sendrecv.send.component;
        int const op = hh.processor(rl, oc);
        if (op != p)
          sends.at(op).push_back(xferinfo);
      }

      vector<xferinfo_t> const recv = alltoallv1(dist::comm(), sends);

      xferinfos.insert(xferinfos.end(), recv.begin(), recv.end());
    }

    if (verbose) {
      ostringstream buf;
      buf << "levelinfos=" << levelinfos;
      CCTK_INFO(buf.str().c_str());
    }

  } // if need new levelinfo

  // Transfer: Loop over all communication phases, all variables, and
  // all regions
  for (comm_state state; not state.done(); state.step()) {

    vector<xferinfo_t>::const_iterator const ib = xferinfos.begin();
    vector<xferinfo_t>::const_iterator const ie = xferinfos.end();
    for (vector<xferinfo_t>::const_iterator ixferinfo = ib; ixferinfo != ie;
         ++ixferinfo) {
      xferinfo_t const &xferinfo = *ixferinfo;

      gh const &hh = *vhh.at(xferinfo.m);

      // Determine components, local components, processes
      int const oc = xferinfo.sendrecv.send.component;
      int const c = xferinfo.sendrecv.recv.component;
      int const olc = hh.get_local_component(rl, oc);
      int const lc = hh.get_local_component(rl, c);
      int const op = hh.processor(rl, oc);
      int const p = hh.processor(rl, c);

      for (int n = 0; n < nvars; ++n) {
        int const vi = vars[n];
        int const gi = CCTK_GroupIndexFromVarI(vi);
        int const v0 = CCTK_FirstVarIndexI(gi);

        ggf &ff = *arrdata.at(gi).at(xferinfo.m).data.at(vi - v0);

        // Get pointers to the variable's data
        gdata *const src =
            hh.is_local(rl, oc) ? ff.data_pointer(tl, rl, olc, ml) : NULL;
        gdata *const dst =
            hh.is_local(rl, c) ? ff.data_pointer(tl, rl, lc, ml) : NULL;

        // Copy
        ibbox const &dstbox = xferinfo.sendrecv.recv.extent;
        ibbox const &srcbox = xferinfo.sendrecv.send.extent;
        islab const &slab = xferinfo.slab;
        gdata::copy_data(dst, state, src, dstbox, srcbox, &slab, p, op);
      }
    }
  }
}

extern "C" CCTK_INT BndPeriodicCarpetVI(CCTK_POINTER_TO_CONST const cctkGH_,
                                        CCTK_INT const size,
                                        CCTK_INT const *restrict const stencil,
                                        CCTK_INT const do_periodic[3],
                                        CCTK_INT const vi) {
  cGH const *restrict const cctkGH = static_cast<cGH const *>(cctkGH_);
  CCTK_INT const vi1 = vi;
  periodic_carpet(cctkGH, size, stencil, do_periodic, &vi1, 1);
  return 0;
}

extern "C" CCTK_INT BndPeriodicCarpetVN(CCTK_POINTER_TO_CONST const cctkGH_,
                                        CCTK_INT const size,
                                        CCTK_INT const *restrict const stencil,
                                        CCTK_INT const do_periodic[3],
                                        char const *restrict const vn) {
  CCTK_INT const vi = CCTK_VarIndex(vn);
  assert(vi >= 0 and vi < CCTK_NumVars());
  BndPeriodicCarpetVI(cctkGH_, size, stencil, do_periodic, vi);
  return 0;
}

extern "C" CCTK_INT BndPeriodicCarpetGI(CCTK_POINTER_TO_CONST const cctkGH_,
                                        CCTK_INT const size,
                                        CCTK_INT const *restrict const stencil,
                                        CCTK_INT const do_periodic[3],
                                        CCTK_INT const gi) {
  cGH const *restrict const cctkGH = static_cast<cGH const *>(cctkGH_);
  int const nvars = CCTK_NumVarsInGroupI(gi);
  vector<CCTK_INT> vis(nvars);
  int const v0 = CCTK_FirstVarIndexI(gi);
  for (int vi = 0; vi < nvars; ++vi)
    vis.at(vi) = v0 + vi;
  periodic_carpet(cctkGH, size, stencil, do_periodic, &vis.front(), nvars);
  return 0;
}

extern "C" CCTK_INT BndPeriodicCarpetGN(CCTK_POINTER_TO_CONST const cctkGH_,
                                        CCTK_INT const size,
                                        CCTK_INT const *restrict const stencil,
                                        CCTK_INT const do_periodic[3],
                                        char const *restrict const gn) {
  int const gi = CCTK_GroupIndex(gn);
  assert(gi >= 0 and gi < CCTK_NumGroups());
  BndPeriodicCarpetGI(cctkGH_, size, stencil, do_periodic, gi);
  return 0;
}

extern "C" void PeriodicCarpet_RegisterBC(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT faces[6];
  faces[0] = faces[1] = periodic or periodic_x;
  faces[2] = faces[3] = periodic or periodic_y;
  faces[4] = faces[5] = periodic or periodic_z;

  CCTK_INT width[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  int ierr =
      GetBoundarySpecification(6, width, is_internal, is_staggered, shiftout);
  if (ierr < 0)
    CCTK_ERROR("Could not get the boundary specification");

  CCTK_INT const handle = SymmetryRegister("periodic");
  if (handle < 0)
    CCTK_ERROR("Could not register periodicity boundary condition");

  ierr = SymmetryRegisterGrid(cctkGH, handle, faces, width);
  if (ierr < 0)
    CCTK_ERROR("Could not register the periodic boundaries -- probably some "
               "other thorn has already registered the same boundary faces for "
               "a different symmetry");

  ierr = SymmetryRegisterGridInterpolator(cctkGH, handle,
                                          PeriodicCarpet_Interpolate);
  if (ierr < 0)
    CCTK_ERROR("Could not register the symmetry interpolator");
}

extern "C" void PeriodicCarpet_ApplyBC(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH);

  CCTK_INT do_periodic[dim];
  assert(dim == 3);
  do_periodic[0] = periodic or periodic_x;
  do_periodic[1] = periodic or periodic_y;
  do_periodic[2] = periodic or periodic_z;

  CCTK_INT const nvars = Boundary_SelectedGVs(cctkGH, 0, 0, 0, 0, 0, 0);
  assert(nvars >= 0);
  if (nvars == 0)
    return;

  vector<CCTK_INT> indices(nvars);
  vector<CCTK_INT> faces(nvars);
  vector<CCTK_INT> widths(nvars);
  vector<CCTK_INT> tables(nvars);
  int const iret =
      Boundary_SelectedGVs(cctkGH, nvars, &indices.front(), &faces.front(),
                           &widths.front(), &tables.front(), 0);
  assert(iret == nvars);

  CCTK_INT width[2 * dim];
  CCTK_INT is_internal[2 * dim];
  CCTK_INT is_staggered[2 * dim];
  CCTK_INT shiftout[2 * dim];
  int ierr = GetBoundarySpecification(2 * dim, width, is_internal, is_staggered,
                                      shiftout);
  if (ierr < 0)
    CCTK_ERROR("Could not get the boundary specification");

  CCTK_INT stencil[dim];
  for (int d = 0; d < dim; ++d) {
    if (do_periodic[d]) {
      assert(width[2 * d] == width[2 * d + 1]);
      stencil[d] = width[2 * d];
    } else {
      stencil[d] = 0;
    }
  }

  for (int n = 0; n < nvars; ++n) {
    int const vi = indices.at(n);
    assert(vi >= 0 and vi < CCTK_NumVars());

    if (verbose) {
      char *const fullname = CCTK_FullName(vi);
      assert(fullname);
      CCTK_VINFO("Applying periodicity boundary conditions to \"%s\"",
                 fullname);
      free(fullname);
    }
  } // for n

  periodic_carpet(cctkGH, dim, stencil, do_periodic, &indices.front(), nvars);
}

namespace CarpetLib {
template vector<xferinfo_t> alltoallv1(MPI_Comm comm,
                                       vector<vector<xferinfo_t> > const &data);
}
