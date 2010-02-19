#include <cassert>
#include <sstream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <dh.hh>

#include <carpet.hh>

#include <loopcontrol.h>

#include "mask_carpet.hh"



namespace CarpetMask {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  void ibbox2iminimax (ibbox const& ext, // component extent
                       ibbox const& box, // this bbox
                       ivect const& lsh, // component's lsh 
                       ivect& imin, ivect& imax)
  {
    ivect const izero = ivect(0);
    
    assert (all ((box.lower() - ext.lower()               ) >= 0));
    assert (all ((box.upper() - ext.lower() + ext.stride()) >= 0));
    assert (all ((box.lower() - ext.lower()               ) % ext.stride() == 0));
    assert (all ((box.upper() - ext.lower() + ext.stride()) % ext.stride() == 0));
    
    imin = (box.lower() - ext.lower()               ) / ext.stride();
    imax = (box.upper() - ext.lower() + ext.stride()) / ext.stride();
    
    assert (all (izero <= imin));
    assert (box.empty() xor all (imin <= imax));
    assert (all (imax <= lsh));
  }
  
  
  
  /**
   * Reduce the weight on the current and the next coarser level to
   * make things consistent.  Set the weight to 0 inside the
   * restriction region of the next coarser level, maybe to 1/2 near
   * the boundary of that region, and also to 1/2 near the
   * prolongation boundary of this level.
   */
  
  void
  CarpetMaskSetup (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (not is_singlemap_mode()) {
      CCTK_WARN (0, "This routine may only be called in singlemap mode");
    }
    
    if (reflevel > 0) {
      
      ivect const ione  = ivect(1);
      
      gh const & hh = *vhh.at(Carpet::map);
      dh const & dd = *vdd.at(Carpet::map);
      
      ibbox const & base = hh.baseextents.at(mglevel).at(reflevel);
      
      ivect const reffact
        = spacereffacts.at(reflevel) / spacereffacts.at(reflevel-1);
      assert (all (reffact == 2));
      
      
      
      // Calculate the union of all refined regions
      ibset active;
      for (int c=0; c<hh.components(reflevel); ++c) {
        // refined |= hh.extent(mglevel,reflevel,c);
        ibset this_active;
        dh::dboxes const& box = dd.boxes.AT(mglevel).AT(reflevel).AT(c);
        dh::dboxes::ibboxs2ibset (box.active, box.numactive, this_active);
        active += this_active;
      }
      active.normalize();
      
      // Calculate the union of all coarse regions
      ibset parent;
      for (int c=0; c<hh.components(reflevel-1); ++c) {
        parent |= hh.extent(mglevel,reflevel-1,c)
          .expand(ivect(reffact-1),ivect(reffact-1)).contracted_for(base);
      }
      parent.normalize();
      
      // Subtract the refined region
      ibset notrefined = parent - active;
      notrefined.normalize();
      
      // Enlarge this set
      ibset enlarged[dim];
      for (int d=0; d<dim; ++d) {
        for (ibset::const_iterator
               bi = notrefined.begin(); bi != notrefined.end(); ++bi)
        {
          if (hh.refcent == vertex_centered) {
            enlarged[d] |= (*bi).expand(ivect::dir(d), ivect::dir(d));
          } else {
            enlarged[d] |= *bi;
          }
        }
        enlarged[d].normalize();
      }
      
      // Intersect with the union of refined regions
      ibset boundaries[dim];
      for (int d=0; d<dim; ++d) {
        boundaries[d] = active & enlarged[d];
        boundaries[d].normalize();
      }
      
      // Subtract the boundaries from the refined region
      ibset refined = active;
      for (int d=0; d<dim; ++d) {
        refined -= boundaries[d];
      }
      refined.normalize();
      
      
      
      // Set prolongation boundaries of this level
      {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          
          DECLARE_CCTK_ARGUMENTS;
          
          ibbox const & ext
            = dd.boxes.at(mglevel).at(reflevel).at(component).exterior;
          
          ibset notactive = ext - active;
          notactive.normalize();
          
          for (int d=0; d<dim; ++d) {
            assert ((notactive & boundaries[d]).empty());
          }
          
          for (ibset::const_iterator bi = notactive.begin();
               bi != notactive.end();
               ++bi)
          {
            ibbox const& box = *bi;
            assert (box <= ext);
            assert (not box.empty());
            
            ivect imin, imax;
            ibbox2iminimax (ext, box, ivect::ref(cctk_lsh), imin, imax);
            
            if (verbose) {
              ostringstream buf;
              buf << "Setting buffer region on level " << reflevel << " to weight 0: " << imin << ":" << imax-ione;
              CCTK_INFO (buf.str().c_str());
            }
            
            // Set weight on the boundary to 0
            assert (dim == 3);
#pragma omp parallel
            LC_LOOP3(CarpetMaskSetup_buffers,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
            {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              weight[ind] = 0.0;
            } LC_ENDLOOP3(CarpetMaskSetup_buffers);
            
          } // for box
          
          for (int d=0; d<dim; ++d) {
            for (ibset::const_iterator bi = boundaries[d].begin();
                 bi != boundaries[d].end();
                 ++bi)
            {
              
              ibbox const & box = (*bi) & ext;
              if (not box.empty()) {
                
                ivect imin, imax;
                ibbox2iminimax (ext, box, ivect::ref(cctk_lsh), imin, imax);
                
                if (verbose) {
                  ostringstream buf;
                  buf << "Setting prolongation boundary on level " << reflevel << " direction " << d << " to weight 1/2: " << imin << ":" << imax-ione;
                  CCTK_INFO (buf.str().c_str());
                }
                
                // Set weight on the boundary to 1/2
                assert (dim == 3);
#pragma omp parallel
                LC_LOOP3(CarpetMaskSetup_prolongation,
                         i,j,k,
                         imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
                {
                  int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                  weight[ind] *= 0.5;
                } LC_ENDLOOP3(CarpetMaskSetup_prolongation);
                
              } // if box not empty
              
            } // for box
          } // for d
          
        } END_LOCAL_COMPONENT_LOOP;
      }
      
      
      
      // Set restriction region on next coarser level
      {
        int const oldreflevel = reflevel;
        int const oldgrouptype = mc_grouptype;
        int const oldmap = Carpet::map;
        leave_singlemap_mode (cctkGH);
        leave_level_mode (cctkGH);
        enter_level_mode (cctkGH, oldreflevel-1);
        enter_singlemap_mode (cctkGH, oldmap, oldgrouptype);
        
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          
          DECLARE_CCTK_ARGUMENTS;
          
          ibbox const & ext
            = dd.boxes.at(mglevel).at(reflevel).at(component).exterior;
          
          for (ibset::const_iterator bi = refined.begin();
               bi != refined.end();
               ++bi)
          {
            
            ibbox const & box = (*bi).contracted_for(ext) & ext;
            if (not box.empty()) {
              
              ivect imin, imax;
              ibbox2iminimax (ext, box, ivect::ref(cctk_lsh), imin, imax);
              
              if (verbose) {
                ostringstream buf;
                buf << "Setting restricted region on level " << reflevel << " to weight 0: " << imin << ":" << imax-ione;
                CCTK_INFO (buf.str().c_str());
              }
              
              // Set weight in the restricted region to 0
              assert (dim == 3);
#pragma omp parallel
              LC_LOOP3(CarpetMaskSetup_restriction,
                       i,j,k,
                       imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                       cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                weight[ind] = 0;
              } LC_ENDLOOP3(CarpetMaskSetup_restriction);
              
            } // if box not empty
              
          } // for box
          
          assert (dim == 3);
          vector<int> mask (cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]);
          
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(CarpetMaskSetup_restriction_boundary_init,
                   i,j,k,
                   0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            mask[ind] = 0;
          } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_init);
          
          for (int d=0; d<dim; ++d) {
            for (ibset::const_iterator bi = boundaries[d].begin();
                 bi != boundaries[d].end();
                 ++bi)
            {
              
              ibbox const & box = (*bi).contracted_for(ext) & ext;
              if (not box.empty()) {
                
                ivect imin, imax;
                ibbox2iminimax (ext, box, ivect::ref(cctk_lsh), imin, imax);
                
                if (verbose) {
                  ostringstream buf;
                  buf << "Setting restriction boundary on level " << reflevel << " direction " << d << " to weight 1/2: " << imin << ":" << imax-ione;
                  CCTK_INFO (buf.str().c_str());
                }
                
                // Set weight on the boundary to 1/2
                assert (dim == 3);
#pragma omp parallel
                LC_LOOP3(CarpetMaskSetup_restriction_boundary_partial,
                         i,j,k,
                         imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
                {
                  int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                  if (mask[ind] == 0) {
                    mask[ind] = 1;
                  }
                  mask[ind] *= 2;
                } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_partial);
                
              } // if box not empty
              
            } // for box
          } // for d
          
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(CarpetMaskSetup_restriction_boundary_apply,
                   i,j,k,
                   0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (mask[ind] > 0) {
              weight[ind] *= 1.0 - 1.0 / mask[ind];
            }
          } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_apply);
          
        } END_LOCAL_COMPONENT_LOOP;
        
        leave_singlemap_mode (cctkGH);
        leave_level_mode (cctkGH);
        enter_level_mode (cctkGH, oldreflevel);
        enter_singlemap_mode (cctkGH, oldmap, oldgrouptype);
      }
      
    } // if reflevel>0
  }
  
} // namespace CarpetMask
