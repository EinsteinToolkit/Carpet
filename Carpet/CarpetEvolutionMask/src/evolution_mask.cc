#include <cassert>
#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

#include "evolution_mask.hh"



namespace CarpetEvolutionMask {
  
  using namespace std;
  using namespace Carpet;
  
  
  void
  CarpetEvolutionMaskSetup (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // ask the flesh what the value of carpet::buffer_width is...
    // this should be fixed at some point...

    CCTK_INT use_outer_buffer_zones = *(CCTK_INT*) CCTK_ParameterGet
      ("use_outer_buffer_zones","Carpet",0);

    CCTK_INT buffer_width;
    if (use_outer_buffer_zones) {
      buffer_width = 0;
    } else {
      buffer_width = *(CCTK_INT*) CCTK_ParameterGet
        ("buffer_width","Carpet",0);
    }
    // cout << "buffer width: " << buffer_width << endl;
    
    if (! is_singlemap_mode()) {
      CCTK_WARN (0, "This routine may only be called in singlemap mode");
    }
    
    if (reflevel > 0) {
      
      ivect const izero = ivect(0);
      ivect const ione  = ivect(1);
      
      gh const & hh = *vhh.at(Carpet::map);
      dh const & dd = *vdd.at(Carpet::map);
      
      ibbox const base =
        hh.baseextents.at(mglevel).at(reflevel).expand(100,100);
      ibbox const coarsebase =
        hh.baseextents.at(mglevel).at(reflevel-1).expand(100,100);
      
      ivect const reffact
        = spacereffacts.at(reflevel) / spacereffacts.at(reflevel-1);
      assert (all (reffact == 2));
      
      // cout << "base: " << base << endl;
      // cout << "coarsebase: " << coarsebase << endl;
      
      // Calculate the union of all refined regions
      ibset refined;
      for (int c=0; c<hh.components(reflevel); ++c) {
	ibbox refcomp = hh.extent(mglevel,reflevel,c);
	b2vect outer_boundary = hh.outer_boundaries(reflevel,c);
	ivect expand_right, expand_left;
	for (int d=0;d<dim;d++) {
	  if (outer_boundary[0][d]) {
	    expand_left[d] = 50;
	  } else {
	    expand_left[d] = 0;
	  }
	  if (outer_boundary[1][d]) {
	    expand_right[d] = 50;
	  } else {
	    expand_right[d] = 0;
	  }
	}
	//	cout << "rl,c,expand_left,expand_right" << " " << reflevel << " " << c
	//	     << " " << expand_left << " " << expand_right << endl;
	refcomp = refcomp.expand(expand_left,expand_right);
	refined |= refcomp;
      }
      refined.normalize();
      

      //      cout << "refined: " << refined << endl;

      // calculate inverse set of refined regions on the current level
      ibset antirefined = base - refined;


      //      cout << "antirefined: " << antirefined << endl;

      // now make antibase larger
      // make refined regions SMALLER
      ivect antishrinkby;
      for (int d=0; d<dim;d++) {
	antishrinkby[d] = cctkGH->cctk_nghostzones[d] + buffer_width;
      }
      ibset antishrunk;
      for (ibset::const_iterator bi = antirefined.begin();
             bi != antirefined.end();
             ++bi)
        {
          antishrunk |= (*bi).expand(antishrinkby, antishrinkby).expanded_for(coarsebase);
        }
      antishrunk.normalize();

      //      cout << "antishrunk1: " << antishrunk << endl;

      // now cut away dangling edges
      antishrunk &= coarsebase;

      //      cout << "antishrunk2: " << antishrunk << endl;

      // cut holes into coarsebase
      ibset shrunk = coarsebase - antishrunk;

      //      cout << "shrunk: " << shrunk << endl;

      // Calculate the union of all coarse regions
      ibset parent;
      for (int c=0; c<hh.components(reflevel-1); ++c) {
        parent |= hh.extent(mglevel,reflevel-1,c);
      }
      parent.normalize();

      //      cout << "parent: " << parent << endl;
      
      // Subtract the refined region
      ibset notrefined = parent - shrunk;
      notrefined.normalize();

      //      cout << "notrefined: " << notrefined << endl;

      ivect enlargeby;
      int stencil_size = dd.prolongation_stencil_size();
      for (int d=0; d<dim;d++) {
	enlargeby[d] = cctkGH->cctk_nghostzones[d] + buffer_width +
	  stencil_size;
      }
      ibset enlarged;
      for (ibset::const_iterator bi = notrefined.begin();
             bi != notrefined.end();
             ++bi)
        {
          enlarged |= (*bi).expand(enlargeby, enlargeby);
        }
      enlarged.normalize();

      //      cout << "enlarged: " << enlarged << endl;
      
      // Intersect with the original union      
      ibset evolveon = parent & enlarged;
      evolveon.normalize();

      ibset notevolveon = parent - evolveon;
      notevolveon.normalize();

      //      cout << "notevolveon: " << notevolveon << endl;
      
      // Set not evolved region on next coarser level
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
          
          for (ibset::const_iterator bi = notevolveon.begin();
               bi != notevolveon.end();
               ++bi)
          {
            
            ibbox const & box = (*bi) & ext;
            if (! box.empty()) {
              
              assert (all ((box.lower() - ext.lower()               ) >= 0));
              assert (all ((box.upper() - ext.lower() + ext.stride()) >= 0));
              assert (all ((box.lower() - ext.lower()               ) % ext.stride() == 0));
              assert (all ((box.upper() - ext.lower() + ext.stride()) % ext.stride() == 0));
              ivect const imin = (box.lower() - ext.lower()               ) / ext.stride();
              ivect const imax = (box.upper() - ext.lower() + ext.stride()) / ext.stride();
              assert (all (izero <= imin));
              assert (box.empty() || all (imin <= imax));
              assert (all (imax <= ivect::ref(cctk_lsh)));
              
              if (verbose) {
                ostringstream buf;
                buf << "Setting restricted region on level " << reflevel << " to weight 0: " << imin << ":" << imax-ione;
                CCTK_INFO (buf.str().c_str());
              }
              
              // Set mask in the restricted region to 0
              assert (dim == 3);
              for (int k=imin[2]; k<imax[2]; ++k) {
                for (int j=imin[1]; j<imax[1]; ++j) {
                  for (int i=imin[0]; i<imax[0]; ++i) {
                    int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                    evolution_mask[ind] = 0;
                  }
                }
              }
              
            } // if box not empty
              
          } // for box
          
        } END_LOCAL_COMPONENT_LOOP;
        
        leave_singlemap_mode (cctkGH);
        leave_level_mode (cctkGH);
        enter_level_mode (cctkGH, oldreflevel);
        enter_singlemap_mode (cctkGH, oldmap, oldgrouptype);
      }
      
    } // if reflevel>0
  }
  
} // namespace CarpetEvolutionMask
