#include <cassert>
#include <istream>
#include <sstream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <dh.hh>

#include <carpet.hh>

#include <loopcontrol.h>

#include "bits.h"



#define LOOP_OVER_NEIGHBOURS(dir)               \
{                                               \
  ivect dir_(-1);                               \
  do {                                          \
    ivect const& dir = dir_;                    \
    {
#define END_LOOP_OVER_NEIGHBOURS                \
    }                                           \
    for (int d_=0; d_<dim; ++d_) {              \
      if (dir_[d_] < +1) {                      \
        ++dir_[d_];                             \
        break;                                  \
      }                                         \
      dir_[d_] = -1;                            \
    }                                           \
  } while (not all (dir_ == -1));               \
}



namespace CarpetMask {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  /**
   * Reduce the weight on the current and the next coarser level to
   * make things consistent. Set the weight to 0 inside the active
   * region of the next coarser level, maybe to 1/2 near the boundary
   * of that region, and also to 1/2 near the prolongation boundary of
   * this level.
   */
  
  extern "C" {
    void
    CarpetMaskSetup (CCTK_ARGUMENTS);
  }
  
  void
  CarpetMaskSetup (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (not is_singlemap_mode()) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "This routine may only be called in singlemap mode");
    }
    
    gh const& hh = *vhh.AT(Carpet::map);
    dh const& dd = *vdd.AT(Carpet::map);
    ibbox const& base = hh.baseextent(mglevel, reflevel);
    
    if (reflevel > 0) {
      ivect const reffact =
        spacereffacts.AT(reflevel) / spacereffacts.AT(reflevel-1);
      assert (all (reffact == 2));
    }
    
    // Set prolongation boundaries and restricted region of this level
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      
      ibbox const& ext =
        dd.light_boxes.AT(mglevel).AT(reflevel).AT(component).exterior;
      ibbox const& owned =
        dd.light_boxes.AT(mglevel).AT(reflevel).AT(component).owned;
      ibbox const world = owned;
      ibset const& active =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).active;
      // There are no prolongation boundaries on the coarsest level;
      // the outer boundary is treated elsewhere
      ibset const not_active = reflevel==0 ? ibset() : world - active;
      ibset const& fine_active =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).fine_active;
      if (verbose) {
        ostringstream buf;
        buf << "Setting prolongation region " << not_active << " on level " << reflevel;
        CCTK_INFO (buf.str().c_str());
        buf << "Setting restriction region " << fine_active << " on level " << reflevel;
        CCTK_INFO (buf.str().c_str());
      }
      
      // Set the weight in the interior of the not_active and the
      // fine_active regions to zero, and set the weight on the
      // boundary of the not_active and fine_active regions to 1/2.
      //
      // For the prolongation region, the "boundary" is the first
      // layer outside of the region. For the restricted region, the
      // "boundary" is the outermost layer of grid points if this
      // layer is aligned with the next coarser (i.e. the current)
      // grid; otherwise, the boundary is empty.
      ibset test_boxes;
      ibset test_cfboxes;
      
      LOOP_OVER_NEIGHBOURS (shift) {
        // In this loop, shift [1,1,1] denotes a convex corner of the
        // region which should be masked out, i.e. a region where only
        // a small bit (1/8) of the region should be masked out.
        // Concave edges are treated implicitly (sequentially), i.e.
        // larger chunks are cut out multiple times: a concave xy edge
        // counts as both an x face and a y face.
        
        ibset boxes  = not_active;
        ibset fboxes = fine_active;
        for (int d=0; d<dim; ++d) {
          ivect const dir = ivect::dir(d);
          fboxes = fboxes.shift(-dir) & fboxes & fboxes.shift(+dir);
        }
        for (int d=0; d<dim; ++d) {
          // Calculate the boundary in direction d
          ivect const dir = ivect::dir(d);
          switch (shift[d]) {
          case -1: {
            // left boundary
            boxes  = boxes.shift (-dir) - boxes;
            fboxes = fboxes.shift(-dir) - fboxes;
            break;
          }
          case 0: {
            // interior
            // do nothing
            break;
          }
          case +1: {
            // right boundary
            boxes  = boxes.shift (+dir) - boxes;
            fboxes = fboxes.shift(+dir) - fboxes;
            break;
          }
          default:
            assert (0);
          }
        }
        boxes &= ext;
        ibset const cfboxes = fboxes.contracted_for(base) & ext;
        test_boxes   |= boxes;
        test_cfboxes |= cfboxes;
        
        if (verbose) {
          ostringstream buf;
          buf << "Setting boundary " << shift << ": prolongation region " << boxes;
          buf << "Setting boundary " << shift << ": restriction region " << fboxes;
          CCTK_INFO (buf.str().c_str());
        }
        
        // Set up a bit mask that keeps the upper (when dir[d]=-1) or
        // lower (when dir[d]=+1) half of the bits in each direction d
        unsigned const bits = BMSK(dim);
        unsigned bmask = 0;
        for (int d=0; d<dim; ++d) {
          for (unsigned b=0; b<bits; ++b) {
            if ((shift[d]==-1 and not BGET(b,d)) or
                (shift[d]==+1 and     BGET(b,d)))
            {
              bmask = BSET(bmask, b);
            }
          }
        }
        
        // Handle prolongation region
        LOOP_OVER_BSET (cctkGH, boxes, box, imin, imax) {
          if (verbose) {
            ostringstream buf;
            buf << "Setting prolongation region "<< imin << ":" << imax-ivect(1) << " on level " << reflevel << " boundary " << shift << " to bmask " << bmask;
            CCTK_INFO (buf.str().c_str());
          }
#pragma omp parallel
          CCTK_LOOP3(CarpetMaskSetup_prolongation,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            iweight[ind] &= bmask;
          } CCTK_ENDLOOP3(CarpetMaskSetup_prolongation);
        } END_LOOP_OVER_BSET;
        
        // Handle restricted region
        LOOP_OVER_BSET (cctkGH, cfboxes, box, imin, imax) {
          if (verbose) {
            ostringstream buf;
            buf << "Setting restricted region "<< imin << ":" << imax-ivect(1) << " on level " << reflevel << " boundary " << shift << " to bmask " << bmask;
            CCTK_INFO (buf.str().c_str());
          }
#pragma omp parallel
          CCTK_LOOP3(CarpetMaskSetup_restriction,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            iweight[ind] &= bmask;
          } CCTK_ENDLOOP3(CarpetMaskSetup_restriction);
        } END_LOOP_OVER_BSET;
        
      } END_LOOP_OVER_NEIGHBOURS;
      
      {
        ibset const boxes = not_active.expand(ivect(1), ivect(1)) & ext;
        ibset fboxes = fine_active;
        for (int d=0; d<dim; ++d) {
          ivect const dir = ivect::dir(d);
          fboxes = fboxes.shift(-dir) & fboxes & fboxes.shift(+dir);
        }
        fboxes = fboxes.expand(ivect(1), ivect(1));
        ibset const cfboxes = fboxes.contracted_for(base) & ext;
        if (not (test_boxes   == boxes  ) or
            not (test_cfboxes == cfboxes))
        {
          cout << "boxes=" << boxes << "\n"
               << "test_boxes=" << test_boxes << "\n"
               << "b-t=" << (boxes-test_boxes) << "\n"
               << "t-b=" << (test_boxes-boxes) << "\n"
               << "cfboxes=" << cfboxes << "\n"
               << "test_cfboxes=" << test_cfboxes << "\n"
               << "b-t=" << (cfboxes-test_cfboxes) << "\n"
               << "t-b=" << (test_cfboxes-cfboxes) << "\n";
        }
        assert (test_boxes   == boxes  );
        assert (test_cfboxes == cfboxes);
      }
      
      
      
    } END_LOCAL_COMPONENT_LOOP;
  }
  
} // namespace CarpetMask
