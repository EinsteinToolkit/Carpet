#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

using namespace std;
using namespace Carpet;



static void
periodic_carpet(cGH const *restrict const cctkGH,
                 int const size,
                 CCTK_INT const *restrict const stencil,
                 CCTK_INT const do_periodic[3],
                 int const *restrict const vars,
                 int const nvars)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Check arguments
  assert(cctkGH);
  assert(stencil);
  assert(vars);
  
  // Early return if possible, because the code below requires at
  // least one variable
  if (nvars==0) return;
  
  // Get and check group info
  cGroup group;
  cGroupDynamicData data;
  for (int i=0; i<nvars; ++i) {
    int const vi = vars[i];
    assert(vi>=0 and vi<CCTK_NumVars());
    int const gi = CCTK_GroupIndexFromVarI(vi);
    assert(gi>=0 and gi<CCTK_NumGroups());
    
    int ierr = CCTK_GroupData(gi, &group);
    assert(not ierr);
    assert(group.grouptype == CCTK_GF);
    assert(group.disttype == CCTK_DISTRIB_DEFAULT);
    assert(group.stagtype == 0);
    int const vartypesize = CCTK_VarTypeSize(group.vartype);
    assert(vartypesize>0);
    
    if (group.dim != size) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "The group \"%s\" has dimension %d, but the given stencil has "
                 "size %d.", CCTK_GroupNameFromVarI(vi), group.dim, size);
    }
    
    ierr = CCTK_GroupDynamicData(cctkGH, gi, &data);
    assert(not ierr);
  }
  
  // Ensure we are in level mode (and not global mode)
  assert(reflevel>=0);
  
  // The cached data structures, containing the communication schedule
  // for each refinement level
  int const ml = mglevel;
  int const rl = reflevel;
  int const tl = 0;
  
  struct xferinfo_t {
    int m;                            // map
    sendrecv_pseudoregion_t sendrecv; // regions and components
    islab slab;                       // slabbing offset
  };
  
  struct levelinfo_t {
    int regridding_epoch;
    vector<xferinfo_t> xferinfos;
    levelinfo_t(): regridding_epoch(-1) {}
  };
  
  static vector<levelinfo_t> levelinfos;
  
  if (reflevels >= (int)levelinfos.size()) {
    levelinfos.resize(reflevels);
  }
  levelinfo_t& levelinfo = levelinfos.at(reflevel);
  
  // Do we need a new levelinfo? (We need a new levelinfo after each
  // regridding.)
  if (levelinfo.regridding_epoch != level_regridding_epochs.at(reflevel)) {
    levelinfo.regridding_epoch = level_regridding_epochs.at(reflevel);
    levelinfo.xferinfos.clear();
    
    // Loop over all components and find out how to fill their
    // periodic boundaries
    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      
      CCTK_INT width[2*dim];
      CCTK_INT is_internal[2*dim];
      CCTK_INT is_staggered[2*dim];
      CCTK_INT shiftout[2*dim];
      int ierr = GetBoundarySpecification
        (2*dim, width, is_internal, is_staggered, shiftout);
      assert(not ierr);
      
      // Domain size
      CCTK_REAL physical_min[dim];
      CCTK_REAL physical_max[dim];
      CCTK_REAL interior_min[dim];
      CCTK_REAL interior_max[dim];
      CCTK_REAL exterior_min[dim];
      CCTK_REAL exterior_max[dim];
      CCTK_REAL spacing[dim];
      ierr = GetDomainSpecification(dim, 
                                    physical_min, physical_max,
                                    interior_min, interior_max,
                                    exterior_min, exterior_max,
                                    spacing);
      assert(not ierr);
      
      // Domain size on this level, in terms of grid points
      ivect npoints;
      for (int d=0; d<dim; ++d) {
        CCTK_REAL const rnpoints =
          (exterior_max[d] - exterior_min[d]) / spacing[d]
          + 1 - (width[2*d] + width[2*d+1]);
        npoints[d] = floor(rnpoints + 0.01);
        // Ensure the domain size is integer
        assert(fabs(npoints[d] - rnpoints) < 0.01);
      }
      
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        
        // Loop over all edges/faces/corners (26 in 3D)
        int const nefc = prod(ivect(3));
        for (int efc=0; efc<nefc; ++efc) {
          ivect face;        // {-1; 0; +1} for {lower, inside, upper}
          {
            int efc1 = efc;
            for (int d=0; d<dim; ++d) {
              face[d] = efc1 % 3 - 1;
              efc1 /= 3;
            }
            assert(efc1==0);
          }
          
          // Does this efc have a periodic symmetry?
          bool have_periodic_face = false;
          for (int d=0; d<dim; ++d) {
            if (face[d] !=0 ) {
              int const f = face[d]>0;
              // We have a periodic symmetry if we have one for any
              // face of this efc
              have_periodic_face |= do_periodic[d] and cctk_bbox[2*d+f];
            }
          }
          
          // We have a periodic efc -- treat it
          if (have_periodic_face) {
            
            // Slabbing offset
            islab slab;
            slab.offset = -face * npoints;
            assert(not all(slab.offset == 0)); // this would be trivial
            
            gh const& hh = *vhh.at(Carpet::map);
            dh const& dd = *vdd.at(Carpet::map);
            dh::light_cboxes const& light_level
              = dd.light_boxes.AT(mglevel).AT(reflevel);
            ibbox const& ext = light_level.at(component).exterior;
            
            // Calculate the bbox for this efc
            ivect lo = ext.lower();
            ivect up = ext.upper();
            ivect const str = ext.stride();
            for (int d=0; d<dim; ++d) {
              switch (face[d]) {
              case -1: up[d] = lo[d] + stencil[d] * str[d]; break;
              case +1: lo[d] = up[d] - stencil[d] * str[d]; break;
              }
            }
            ibbox const dst_bbox(lo, up, str);
            assert(dst_bbox.is_contained_in(ext));
            
            // Determine source bbox
            ibbox const src_bbox = dst_bbox.shift(slab.offset);
            // Ensure we copy from within the domain, not from
            // boundary points
            assert(src_bbox.is_contained_in(hh.baseextent(mglevel, reflevel)));
            
            // Find source components for all points
            ibset src_points(src_bbox);
            for (int oc=0; oc<hh.components(reflevel); ++oc) {
              // We cannot copy only from the interior (the "good"
              // points), because this excludes the ghost zones, and
              // we may need to copy from prolongation boundaries.
              // (This assumes that synchronisation and prolongation
              // occurs before periodic boundaries are applied!)
              ibset const intersection =
                src_points & light_level.at(oc).exterior;
              if (not intersection.empty()) {
                src_points -= intersection;
                
                // Loop over the intersection and insert respective
                // xferinfos
                ibset::const_iterator const ib = intersection.begin();
                ibset::const_iterator const ie = intersection.end();
                for (ibset::const_iterator ireg = ib; ireg != ie; ++ireg) {
                  ibbox const src_reg = *ireg;
                  ibbox const dst_reg = src_reg.shift(-slab.offset);
                  xferinfo_t xferinfo;
                  xferinfo.m = Carpet::map;
                  xferinfo.sendrecv.send.extent    = src_reg;
                  xferinfo.sendrecv.send.component = oc;
                  xferinfo.sendrecv.recv.extent    = dst_reg;
                  xferinfo.sendrecv.recv.component = component;
                  xferinfo.slab = slab;
                  levelinfo.xferinfos.push_back(xferinfo);
                }
              }
            } // for oc
            // Ensure we know how to handle all points
            assert(src_points.empty());
            
          } // if have_periodic_face
        } // for efc
        
      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
    
    // Tell the source processes what they have to send
    
    
  } // if need new levelinfo
  
  
  
  // Transfer: Loop over all communication phases, all variables, and
  // all regions
  for (comm_state state; not state.done(); state.step()) {
    for (int n=0; n<nvars; ++n) {
      int const vi = vars[n];
      int const gi = CCTK_GroupIndexFromVarI(vi);
      int const v0 = CCTK_FirstVarIndexI(gi);
      
      vector<xferinfo_t>::const_iterator const ib = levelinfo.xferinfos.begin();
      vector<xferinfo_t>::const_iterator const ie = levelinfo.xferinfos.end();
      for (vector<xferinfo_t>::const_iterator
             ixferinfo = ib; ixferinfo != ie; ++ixferinfo)
      {
        xferinfo_t const& xferinfo = *ixferinfo;
        gh const& hh = *vhh.at(xferinfo.m);
        ggf& ff = *arrdata.at(gi).at(xferinfo.m).data.at(vi-v0);
        
        // Determine components, local components, processes
        int const oc = xferinfo.sendrecv.send.component;
        int const c  = xferinfo.sendrecv.recv.component;
        int const olc = hh.get_local_component(rl, oc);
        int const lc  = hh.get_local_component(rl, c );
        int const op = hh.processor(rl, oc);
        int const p  = hh.processor(rl, c );
        
        // Get pointers to the variable's data
        gdata *const src = ff(tl, rl, olc, ml);
        gdata *const dst = ff(tl, rl, lc , ml);
        
        // Copy
        ibbox const& box = xferinfo.sendrecv.recv.extent;
        islab const& slab = xferinfo.slab;
        dst->copy_from(state, src, box, &slab, p, op);
      }
    }
  }
}



extern "C" CCTK_INT
BndPeriodicCarpetVI(CCTK_POINTER_TO_CONST const cctkGH_,
                     CCTK_INT const size,
                     CCTK_INT const *restrict const stencil,
                     CCTK_INT const do_periodic[3],
                     CCTK_INT const vi)
{
  cGH const *restrict const cctkGH = static_cast<cGH const*>(cctkGH_);
  int const vi1 = vi;
  periodic_carpet(cctkGH, size, stencil, do_periodic, &vi1, 1);
  return 0;
}

extern "C" CCTK_INT
BndPeriodicCarpetVN(CCTK_POINTER_TO_CONST const cctkGH_,
                     CCTK_INT const size,
                     CCTK_INT const *restrict const stencil,
                     CCTK_INT const do_periodic[3],
                     char const *restrict const vn)
{
  int const vi = CCTK_VarIndex(vn);
  assert(vi>=0 and vi<CCTK_NumVars());
  BndPeriodicCarpetVI(cctkGH_, size, stencil, do_periodic, vi);
  return 0;
}

extern "C" CCTK_INT
BndPeriodicCarpetGI(CCTK_POINTER_TO_CONST const cctkGH_,
                     CCTK_INT const size,
                     CCTK_INT const *restrict const stencil,
                     CCTK_INT const do_periodic[3],
                     CCTK_INT const gi)
{
  cGH const *restrict const cctkGH = static_cast<cGH const*>(cctkGH_);
  int const nvars = CCTK_NumVarsInGroupI(gi);
  vector<int> vis(nvars);
  int const v0 = CCTK_FirstVarIndexI(gi);
  for (int vi=0; vi<nvars; ++vi) vis.at(vi) = v0+vi;
  periodic_carpet(cctkGH, size, stencil, do_periodic, &vis.front(), nvars);
  return 0;
}

extern "C" CCTK_INT
BndPeriodicCarpetGN(CCTK_POINTER_TO_CONST const cctkGH_,
                     CCTK_INT const size,
                     CCTK_INT const *restrict const stencil,
                     CCTK_INT const do_periodic[3],
                     char const *restrict const gn)
{
  int const gi = CCTK_GroupIndex(gn);
  assert(gi>=0 and gi<CCTK_NumGroups());
  BndPeriodicCarpetGI(cctkGH_, size, stencil, do_periodic, gi);
  return 0;
}



extern "C" void
PeriodicCarpet_RegisterBC(CCTK_ARGUMENTS)
{
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
  int ierr = GetBoundarySpecification
    (6, width, is_internal, is_staggered, shiftout);
  if (ierr < 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "Could not get the boundary specification");
  }
  
  CCTK_INT const handle = SymmetryRegister("periodic");
  if (handle < 0) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Could not register periodicity boundary condition");
  }
  
  ierr = SymmetryRegisterGrid(cctkGH, handle, faces, width);
  if (ierr < 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "Could not register the periodic boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
  }
}



extern "C" void
PeriodicCarpet_ApplyBC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  assert(cctkGH);
  
  CCTK_INT do_periodic[dim];
  assert(dim==3);
  do_periodic[0] = periodic or periodic_x;
  do_periodic[1] = periodic or periodic_y;
  do_periodic[2] = periodic or periodic_z;
  
  CCTK_INT const nvars = Boundary_SelectedGVs(cctkGH, 0, 0, 0, 0, 0, 0);
  assert(nvars>=0);
  if (nvars==0) return;
  
  vector<int> indices(nvars);
  vector<int> faces(nvars);
  vector<int> widths(nvars);
  vector<int> tables(nvars);
  int const iret = Boundary_SelectedGVs
    (cctkGH, nvars,
     &indices.front(), &faces.front(), &widths.front(), &tables.front(), 0);
  assert(iret == nvars);
  
  CCTK_INT width[2*dim];
  CCTK_INT is_internal[2*dim];
  CCTK_INT is_staggered[2*dim];
  CCTK_INT shiftout[2*dim];
  int ierr = GetBoundarySpecification
    (2*dim, width, is_internal, is_staggered, shiftout);
  if (ierr < 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "Could not get the boundary specification");
  }
  
  CCTK_INT stencil[dim];
  for (int d=0; d<dim; ++d) {
    assert(width[d] == width[d+1]);
    stencil[d] = width[2*d];
  }
  
  for (int n=0; n<nvars; ++n) {
    int const vi = indices.at(n);
    assert(vi>=0 and vi<CCTK_NumVars());
    
    if (verbose) {
      char *const fullname = CCTK_FullName(vi);
      assert(fullname);
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Applying periodicity boundary conditions to \"%s\"",
                 fullname);
      free(fullname);
    }
  } // for n
  
  periodic_carpet(cctkGH, dim, stencil, do_periodic, &indices.front(), nvars);
}
