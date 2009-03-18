#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "dh.hh"
#include "gh.hh"
#include "region.hh"
#include "vect.hh"

#include "carpet.hh"
#include "modes.hh"

#define DEBUG false             // false or true



namespace Carpet {
  
  using namespace std;
  
  
  
  // Helper routines for spliting regions automatically
  
  // The cost for a region, assuming a cost of 1 per interior point
  static rvect
  cost (region_t const & reg)
  {
    DECLARE_CCTK_PARAMETERS;
    static rvect costfactor;
    static bool initialised = false;
    if (not initialised) {
      costfactor = rvect(1.0);
      if (dim > 0) costfactor[0] = 1.0 / aspect_ratio_x;
      if (dim > 1) costfactor[1] = 1.0 / aspect_ratio_y;
      if (dim > 2) costfactor[2] = 1.0 / aspect_ratio_z;
    }
    if (reg.extent.empty()) return rvect(0);
    return rvect (reg.extent.shape() / reg.extent.stride()) * costfactor;
  }
  
  
  
  static void
  SplitRegionsMaps_Automatic_Recursively (bvect    const   & dims,
                                          int      const     firstproc,
                                          int      const     nprocs,
                                          region_t         & superreg,
                                          vector<region_t> & newregs);
  static void
  SplitRegions_AsSpecified (cGH const * cctkGH,
                            vector<region_t> & superregs,
                            vector<region_t> & regs);
  
  
  
  void
  CheckRegions (gh::mregs const & regsss)
  {
    // At least one multigrid level
    if (regsss.size() == 0) {
      CCTK_WARN (0, "I cannot set up a grid hierarchy with zero multigrid levels.");
    }
    assert (regsss.size() > 0);
    // At least one refinement level
    if (regsss.at(0).size() == 0) {
      CCTK_WARN (0, "I cannot set up a grid hierarchy with zero refinement levels.");
    }
    assert (regsss.at(0).size() > 0);
    // At most maxreflevels levels
    if ((int)regsss.at(0).size() > maxreflevels) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "I cannot set up a grid hierarchy with more than Carpet::max_refinement_levels refinement levels.  I found Carpet::max_refinement_levels=%d, while %d levels were requested.",
                  (int)maxreflevels, (int)regsss.at(0).size());
    }
    assert ((int)regsss.at(0).size() <= maxreflevels);
    for (int ml=0; ml<(int)regsss.size(); ++ml) {
      int num_regions = 0;
      for (int rl=0; rl<(int)regsss.at(0).size(); ++rl) {
        // No empty levels
        // (but allow some empty maps)
        // assert (regsss.at(ml).at(rl).size() > 0);
        num_regions += regsss.at(ml).at(rl).size();
        for (int c=0; c<(int)regsss.at(ml).at(rl).size(); ++c) {
          // Check sizes
          // (but allow processors with zero grid points)
          // assert (all(regsss.at(rl).at(c).at(ml).extent.lower() <=
          //             regsss.at(rl).at(c).at(ml).extent.upper()));
          // Check strides
          ivect const str =
            (maxspacereflevelfact / spacereffacts.at(rl) * ipow(mgfact, ml));
          assert (all(regsss.at(ml).at(rl).at(c).extent.stride() % str == 0));
          // Check alignments
          assert (all(regsss.at(ml).at(rl).at(c).extent.lower() % str == 0));
          assert (all(regsss.at(ml).at(rl).at(c).extent.upper() % str == 0));
        }
      }
      // No empty levels
      assert (num_regions > 0);
    }
  }
  
  

  bool
  Regrid (cGH const * const cctkGH,
          bool const force_recompose)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("Regridding level %d...", reflevel);
    
    assert (is_level_mode());
    
    bool const have_regrid     = CCTK_IsFunctionAliased ("Carpet_Regrid");
    bool const have_regridmaps = CCTK_IsFunctionAliased ("Carpet_RegridMaps");
    bool const use_regridmaps = regrid_in_level_mode and have_regridmaps;
    
    if (not use_regridmaps and not have_regrid) {
      static bool didtell = false;
      if (maxreflevels > 1 and not didtell) {
	CCTK_WARN (2, "The regridding routine Carpet_Regrid has not been provided.  There will be no regridding.  Maybe you forgot to activate a regridding thorn?");
	didtell = true;
      }
      return false;
    }
    
    if (regrid_in_level_mode and not have_regridmaps) {
      static bool didtell = false;
      if (maps > 1 and maxreflevels > 1 and not didtell) {
	CCTK_WARN (2, "The regridding routine Carpet_RegridMaps has not been provided.  Regridding will be performed in singlemap mode instead of level mode.");
	didtell = true;
      }
    }
    
    
    
    bool did_change = false;
    
    if (not use_regridmaps) {
      
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        
        gh::rregs superregss = vhh.at(map)->superregions;
        gh::mregs regsss;
        
        // Check whether to recompose
        CCTK_INT const do_recompose =
          Carpet_Regrid (cctkGH, & superregss, & regsss, force_recompose);
        assert (do_recompose >= 0);
        did_change = did_change or do_recompose;
        
        if (do_recompose) {
          RegridMap (cctkGH, map, superregss, regsss);
        }
        
      } END_MAP_LOOP;
      
    } else {
      
      vector<gh::rregs> superregsss (maps);
      vector<gh::mregs> regssss (maps);
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        superregsss.at(map) = vhh.at(map)->superregions;
      } END_MAP_LOOP;
      
      // Check whether to recompose
      CCTK_INT const do_recompose =
        Carpet_RegridMaps (cctkGH, & superregsss, & regssss, force_recompose);
      assert (do_recompose >= 0);
      did_change = did_change or do_recompose;
      
      if (do_recompose) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          gh::rregs const & superregss = superregsss.at(map);
          gh::mregs const & regsss = regssss.at(map);
          RegridMap (cctkGH, map, superregss, regsss);
        } END_MAP_LOOP;
      }
      
    } // if regrid in level mode
    
    
    
    if (did_change) {
      
      PostRegrid (cctkGH);
      
    } // if did change
    
    
    
    return did_change;
  }
  
  
  
  void
  RegridMap (cGH const * const cctkGH,
             int const m,
             gh::rregs const & superregss,
             gh::mregs const & regsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Regridding map %d...", m);

    
#warning "TODO: keep levels fixed here"
#if 0
    //
    // Keep this level fixed if it is not evolved
    //
    if (regrid_only_evolved_levels) {
      int type;
      bool const use_tapered_grids =
        * static_cast<CCTK_INT const *>
        (CCTK_ParameterGet ("use_tapered_grids", "Carpet", &type));
      assert (type == PARAMETER_BOOLEAN);
      
      int const do_every =
        use_tapered_grids ?
        maxtimereflevelfact / timereffacts.at(max(0,rl-1)):
        maxtimereflevelfact / timereffacts.at(      rl   );
      
      bool const regrid_this_level =
        (cctkGH->cctk_iteration - 1) % do_every == 0;
      
      if (not regrid_this_level) {
        // Set regions from current grid structure
        regions.at(rl) = ...;
      }
    }
#endif
    
    // Check the regions
    CheckRegions (regsss);
    // TODO: check also that the current and all coarser levels did
    // not change
    
    // Regrid
    vhh.at(m)->regrid (superregss, regsss);
    
    // Write grid structure to file
    OutputGridStructure (cctkGH, m, regsss);
    OutputGridCoordinates (cctkGH, m, regsss);
#warning "TODO: output superregss"
    
    if (verbose or veryverbose) {
      OutputGrids (cctkGH, m, * vhh.at(m), * vdd.at(m));
    }
    
    Waypoint ("Done regridding map %d.", m);
  }
  
  
  
  void
  PostRegrid (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Calculate new number of levels
    int const oldreflevels = reflevels;
    reflevels = vhh.at(0)->reflevels();
    for (int m=0; m<maps; ++m) {
      assert (vhh.at(m)->reflevels() == reflevels);
    }
    
    // One cannot switch off the current level
    assert (reflevels > reflevel);
    
    // Set new number of active time levels
    for (int n=0; n<CCTK_NumGroups(); ++n) {
      int const grouptype = CCTK_GroupTypeI (n);
      if (grouptype == CCTK_GF) {
        for (int ml=0; ml<mglevels; ++ml) {
          groupdata.at(n).activetimelevels.at(ml).resize
            (reflevels, groupdata.at(n).activetimelevels.at(ml).at(0));
        }
      }
    }
    
    // Set new level times
    for (int ml=0; ml<mglevels; ++ml) {
      leveltimes.at(ml).resize
        (reflevels, leveltimes.at(ml).at(oldreflevels-1));
    }
    
    OutputGridStatistics (cctkGH);
  }
  
  
  
  bool
  Recompose (cGH const * const cctkGH,
             int const rl,
             bool const do_init)
  {
    bool did_recompose = false;
    
    for (int m=0; m<maps; ++m) {
      Waypoint ("Recomposing the grid hierarchy for map %d level %d...", m, rl);
      
      assert (rl>=0 and rl<vhh.at(m)->reflevels());
      did_recompose |= vhh.at(m)->recompose (rl, do_init);
      
      Waypoint ("Done recomposing the grid hierarchy for map %d level %d.",
                m, rl);
    }
    return did_recompose;
  }
  
  
  
  void
  OutputGrids (cGH const * const cctkGH,
               int const m,
               gh const & hh,
               dh const & dd)
  {
    CCTK_INFO ("Grid structure (grid points):");
    for (int ml=0; ml<hh.mglevels(); ++ml) {
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        for (int c=0; c<hh.components(rl); ++c) {
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.at(rl);
          const ibbox   ext = hh.extent(ml,rl,c);
          const ivect & lower  = ext.lower();
          const ivect & upper  = ext.upper();
          const ivect & stride = ext.stride();
          assert (all(lower * levfact % maxspacereflevelfact == 0));
          assert (all(upper * levfact % maxspacereflevelfact == 0));
          assert (all(((upper - lower) * levfact / maxspacereflevelfact)
                      % convfact == 0));
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << "proc "
               << hh.processor(rl,c)
               << "   "
               << lower / stride
               << " : "
               << upper / stride
               << "   ("
               << (upper - lower) / stride + 1
               << ") "
               << prod ((upper - lower) / stride + 1)
               << endl;
        }
      }
    }
    
    CCTK_INFO ("Grid structure (boundaries):");
    for (int rl=0; rl<hh.reflevels(); ++rl) {
      for (int c=0; c<hh.components(rl); ++c) {
        cout << "   [" << rl << "][" << m << "][" << c << "]"
             << "   bbox: "
             << hh.outer_boundaries(rl,c)
             << endl;
      }
    }
    
    CCTK_INFO ("Grid structure (coordinates):");
    for (int ml=0; ml<hh.mglevels(); ++ml) {
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        for (int c=0; c<hh.components(rl); ++c) {
          const rvect origin = domainspecs.at(m).exterior_min;
          const rvect delta  = (domainspecs.at(m).exterior_max - domainspecs.at(m).exterior_min) / rvect (domainspecs.at(m).npoints - 1);
          const ibbox ext = hh.extent(ml,rl,c);
          const ivect & lower = ext.lower();
          const ivect & upper = ext.upper();
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.at(rl);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << origin + delta * rvect(lower) / rvect(maxspacereflevelfact)
               << " : "
               << origin + delta * rvect(upper) / rvect(maxspacereflevelfact)
               << " : "
               << delta * rvect(convfact) / rvect(levfact) << endl;
        }
      }
    }
    
    CCTK_INFO ("Grid structure (coordinates, including ghosts):");
    for (int ml=0; ml<hh.mglevels(); ++ml) {
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        for (int c=0; c<hh.components(rl); ++c) {
          const rvect origin = domainspecs.at(m).exterior_min;
          const rvect delta  = (domainspecs.at(m).exterior_max - domainspecs.at(m).exterior_min) / rvect (domainspecs.at(m).npoints - 1);
          const ivect lower = dd.boxes.at(ml).at(rl).at(c).exterior.lower();
          const ivect upper = dd.boxes.at(ml).at(rl).at(c).exterior.upper();
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.at(rl);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << origin + delta * rvect(lower) / rvect(maxspacereflevelfact)
               << " : "
               << origin + delta * rvect(upper) / rvect(maxspacereflevelfact)
               << " : "
               << delta * rvect(convfact) / rvect(levfact) << endl;
        }
      }
    }
    
    CCTK_INFO ("Grid statistics:");
    const int oldprecision = cout.precision();
    const ios_base::fmtflags oldflags = cout.flags();
    cout.setf (ios::fixed);
    for (int ml=0; ml<hh.mglevels(); ++ml) {
      CCTK_REAL coarsevolume = 0;
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        
        const CCTK_REAL basevolume = hh.baseextents.AT(0).AT(0).size();
        CCTK_REAL countvolume = 0;
        CCTK_REAL totalvolume = 0;
        CCTK_REAL totalvolume2 = 0;
        
        for (int c=0; c<hh.components(rl); ++c) {
          const CCTK_REAL volume = hh.extent(ml,rl,c).size();
          ++ countvolume;
          totalvolume += volume;
          totalvolume2 += ipow(volume, 2);
        }
        
        const CCTK_REAL avgvolume = totalvolume / countvolume;
        const CCTK_REAL stddevvolume
          = sqrt (max (CCTK_REAL(0),
                       totalvolume2 / countvolume - ipow (avgvolume, 2)));
        
        for (int c=0; c<hh.components(rl); ++c) {
          const CCTK_REAL volume = hh.extent(ml,rl,c).size();
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   volume: " << setprecision(0) << volume
               << "   of parent: " << setprecision(1) << 100 * volume / totalvolume << "%"
               << "   of domain: " << setprecision(1) << 100 * volume / basevolume << "%"
               << endl;
        }
        
        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   average volume: " << setprecision(0) << avgvolume
             << "   of parent: " << setprecision(1) << 100 * avgvolume / totalvolume << "%"
             << "   of domain: " << setprecision(1) << 100 * avgvolume / basevolume << "%"
             << endl;
        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   standard deviation: " << setprecision(0) << stddevvolume
             << "   of parent: " << setprecision(1) << 100 * stddevvolume / totalvolume << "%"
             << "   of domain: " << setprecision(1) << 100 * stddevvolume / basevolume << "%"
             << endl;
        
        // TODO: ghost points vs. interior points (and boundary
        // points)
        
        CCTK_REAL countquadrupole = 0;
        CCTK_REAL minquadrupole = 1;
        CCTK_REAL totalquadrupole = 0;
        CCTK_REAL totalquadrupole2 = 0;
        
        for (int c=0; c<hh.components(rl); ++c) {
          const ibbox ext = hh.extent(ml,rl,c);
          const ivect shape = ext.shape();
          const ivect stride = ext.stride();
          const CCTK_REAL minlength = minval (rvect (shape) / rvect (stride));
          const CCTK_REAL maxlength = maxval (rvect (shape) / rvect (stride));
          const CCTK_REAL quadrupole = minlength / maxlength;
          ++ countquadrupole;
          minquadrupole = min (minquadrupole, quadrupole);
          totalquadrupole += quadrupole;
          totalquadrupole2 += ipow (quadrupole, 2);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   length ratio: " << setprecision(3) << quadrupole
               << endl;
        }
        
        const CCTK_REAL avgquadrupole = totalquadrupole / countquadrupole;
        const CCTK_REAL stddevquadrupole
          = sqrt (max (CCTK_REAL(0),
                       (totalquadrupole2 / countquadrupole
                        - ipow (avgquadrupole, 2))));
        
        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   average length ratio: " << setprecision(3) << avgquadrupole
             << "   standard deviation: " << setprecision(3) << stddevquadrupole
             << endl;
        
        // TODO: processor distribution, average load, std deviation
        
        coarsevolume = totalvolume * prod (rvect (spacereflevelfact));
      }
    }
    cout.precision (oldprecision);
    cout.setf (oldflags);
    
  }
  
  
  
  void
  OutputGridStructure (cGH const * const cctkGH,
                       int const m,
                       gh::mregs const & regsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Output only on the root processor
    if (CCTK_MyProc(cctkGH) != 0) return;
    
    // Output only if output is desired
    if (CCTK_EQUALS (grid_structure_filename, "")) return;
    
    // Create the output directory
    CCTK_CreateDirectory (0755, out_dir);
    
    ostringstream filenamebuf;
    filenamebuf << out_dir << "/" << grid_structure_filename;
    // we need a persistent temporary here
    string filenamestr = filenamebuf.str();
    const char * filename = filenamestr.c_str();
    
    ofstream file;
    
    static bool do_truncate = true;
    
    if (do_truncate) {
      do_truncate = false;
      struct stat fileinfo;
      if (IO_TruncateOutputFiles (cctkGH)
	  or stat(filename, &fileinfo)!=0) {
	file.open (filename, ios::out | ios::trunc);
	assert (file.good());
	file << "# grid structure" << eol
	     << "# format: map reflevel component mglevel   processor bounding-box is-outer-boundary" << eol;
	assert (file.good());
      }
    }
    if (not file.is_open()) {
      file.open (filename, ios::app);
      assert (file.good());
    }
    
    file << "iteration " << cctkGH->cctk_iteration << eol;
    file << "maps " << maps << eol;
    file << m << " mglevels " << regsss.size() << eol;
    for (int ml=0; ml<(int)regsss.size(); ++ml) {
      file << m << " " << ml << " reflevels " << regsss.at(ml).size() << eol;
      for (int rl=0; rl<(int)regsss.at(ml).size(); ++rl) {
        file << m << " " << ml << " " << rl << " components " << regsss.at(ml).at(rl).size() << eol;
        for (int c=0; c<(int)regsss.at(ml).at(rl).size(); ++c) {
          file << m << " " << ml << " " << rl << " " << c << "   "
               << regsss.at(ml).at(rl).at(c).processor << " "
               << regsss.at(ml).at(rl).at(c).extent << " "
               << regsss.at(ml).at(rl).at(c).outer_boundaries << eol;
        }
      }
    }
    file << eol;
    
    file.close();
    assert (file.good());
  }
  
  
  
  void
  OutputGridCoordinates (cGH const * const cctkGH,
                         int const m,
                         gh::mregs const & regsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Output only on the root processor
    if (CCTK_MyProc(cctkGH) != 0) return;
    
    // Output only if output is desired
    if (CCTK_EQUALS (grid_coordinates_filename, "")) return;
    
    // Create the output directory
    CCTK_CreateDirectory (0755, out_dir);
    
    ostringstream filenamebuf;
    filenamebuf << out_dir << "/" << grid_coordinates_filename;
    // we need a persistent temporary here
    string filenamestr = filenamebuf.str();
    const char * filename = filenamestr.c_str();
    
    ofstream file;
    
    static bool do_truncate = true;
    
    if (do_truncate) {
      do_truncate = false;
      struct stat fileinfo;
      if (IO_TruncateOutputFiles (cctkGH)
	  or stat(filename, &fileinfo)!=0) {
	file.open (filename, ios::out | ios::trunc);
	assert (file.good());
	file << "# grid coordinates" << eol
	     << "# format: map reflevel region mglevel   bounding-box" << eol;
	assert (file.good());
      }
    }
    if (not file.is_open()) {
      file.open (filename, ios::app);
      assert (file.good());
    }
    
    //
    // Find extent of domain
    //
    
    // This requires that CoordBase is used
#warning "TODO: check this (for Carpet, and maybe also for CartGrid3D)"
#warning "TODO: (the check that these two are consistent should be in Carpet)"
    
    typedef vect<vect<CCTK_INT,2>,3> jjvect;
    jjvect nboundaryzones;
    jjvect is_internal;
    jjvect is_staggered;
    jjvect shiftout;
    check (not GetBoundarySpecification
           (2*dim,
            & nboundaryzones[0][0],
            & is_internal[0][0],
            & is_staggered[0][0],
            & shiftout[0][0]));
    
    rvect physical_lower, physical_upper;
    rvect interior_lower, interior_upper;
    rvect exterior_lower, exterior_upper;
    rvect spacing;
    check (not GetDomainSpecification
           (dim,
            & physical_lower[0], & physical_upper[0],
            & interior_lower[0], & interior_upper[0],
            & exterior_lower[0], & exterior_upper[0],
            & spacing[0]));
    
    // Adapt spacing for convergence level
#warning "TODO: take ml into account"
    spacing *= ipow ((CCTK_REAL) mgfact, basemglevel);
    
    check (not ConvertFromPhysicalBoundary
           (dim,
            & physical_lower[0], & physical_upper[0],
            & interior_lower[0], & interior_upper[0],
            & exterior_lower[0], & exterior_upper[0],
            & spacing[0]));
    
    // Affine transformation between index space and coordinate space
    rvect const origin = exterior_lower;
    rvect const scale =
      rvect (vhh.at(m)->baseextents.at(0).at(0).stride()) / spacing;
    
    
    
    file << "iteration " << cctkGH->cctk_iteration << eol;
    file << "maps " << maps << eol;
    file << m << " mglevels " << regsss.size() << eol;
    for (int ml=0; ml<(int)regsss.size(); ++ml) {
      file << m << " " << ml << " reflevels " << regsss.at(ml).size() << eol;
      for (int rl=0; rl<(int)regsss.at(ml).size(); ++rl) {
        ibset extents;
        for (int c=0; c<(int)regsss.at(ml).at(rl).size(); ++c) {
          extents += regsss.at(ml).at(rl).at(c).extent;
        }
        extents.normalize();
        file << m << " " << ml << " " << rl << " regions " << extents.setsize() << eol;
        int c=0;
        for (ibset::const_iterator
               bi = extents.begin(); bi != extents.end(); ++ bi)
        {
#if 0
          ibbox const & ext       = * bi;
          ibbox const & baseext   = vhh.at(m)->baseextents.at(ml).at(rl);
          ibbox const & coarseext = vhh.at(m)->baseextents.at(ml).at(0 );
          
          // This is nice, but is wrong since CartGrid3D has not yet
          // initialised the coordinates
          ivect const cctk_levfac       = spacereffacts.at(rl);
          ivect const cctk_levoff       = baseext.lower() - coarseext.lower();
          ivect const cctk_levoffdenom  = baseext.stride();
          
          ivect const cctk_lbnd =
            (ext.lower() - baseext.lower()) / ext.stride();
          ivect const cctk_ubnd =
            (ext.upper() - baseext.lower()) / ext.stride();
          
          rvect const cctk_origin_space =
            origin_space.at(m).at(ml);
          rvect const cctk_delta_space  =
            delta_space.at(m) * rvect (ipow (mgfact, ml));
          
          rvect const CCTK_ORIGIN_SPACE =
            cctk_origin_space +
            cctk_delta_space / rvect (cctk_levfac) *
            rvect (cctk_levoff) / rvect (cctk_levoffdenom);
          rvect const CCTK_DELTA_SPACE =
            cctk_delta_space / rvect (cctk_levfac);
          
          rvect const rlower =
            CCTK_ORIGIN_SPACE + rvect (cctk_lbnd) * CCTK_DELTA_SPACE;
          rvect const rupper =
            CCTK_ORIGIN_SPACE + rvect (cctk_ubnd) * CCTK_DELTA_SPACE;
          rvect const rdelta =
            CCTK_DELTA_SPACE;
#endif
          
          ibbox const & ext       = * bi;
          
          rvect const rlower = origin + rvect (ext.lower()) / scale;
          rvect const rupper = origin + rvect (ext.upper()) / scale;
          rvect const rdelta = rvect (ext.stride()) / scale;
          
          rbbox const rb (rlower, rupper, rdelta);
          file << m << " " << ml << " " << rl << " " << c << "   " << rb << eol;
          ++c;
        }
      }
    }
    file << eol;
    
    file.close();
    assert (file.good());
  }
  
  
  
  // TODO: this routine should go into CarpetRegrid (except maybe
  // SplitRegions_AlongZ for grid arrays)
  void
  SplitRegions (cGH const * const cctkGH,
                vector<region_t> & superregs,
                vector<region_t> & regs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (regs.empty());
    if (CCTK_EQUALS (processor_topology, "along-z")) {
      SplitRegions_AlongZ (cctkGH, superregs, regs);
    } else if (CCTK_EQUALS (processor_topology, "along-dir")) {
      SplitRegions_AlongDir (cctkGH, superregs, regs, split_direction);
    } else if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegions_Automatic (cctkGH, superregs, regs);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      SplitRegions_AsSpecified (cctkGH, superregs, regs);
    } else {
      assert (0);
    }
  }
  
  
  
  void
  OutputGridStatistics (cGH const * const cctkGH)
  {
    // Grid array statistics
    int num_gfs = 0;
    int num_arrays = 0;
    CCTK_REAL num_active_array_points = 0;
    CCTK_REAL num_total_array_points  = 0;
    for (int g=0; g<CCTK_NumGroups(); ++g) {
      int const num_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
      switch (CCTK_GroupTypeI(g)) {
      case CCTK_SCALAR:
      case CCTK_ARRAY: {
        int const num_vars = CCTK_NumVarsInGroupI (g);
        num_arrays += num_tl * num_vars;
        gh const * const hh = arrdata.AT(g).AT(0).hh;
        dh const * const dd = arrdata.AT(g).AT(0).dd;
        for (int c=0; c<hh->components(0); ++c) {
          dh::dboxes const & b = dd->boxes.AT(0).AT(0).AT(c);
          num_active_array_points += num_tl * num_vars * b.active.size();
          num_total_array_points  += num_tl * num_vars * b.exterior.size();
        }
        break;
      }
      case CCTK_GF:
        num_gfs += num_tl;
        break;
      default:
        assert (0);
      }
    }
    
    // Grid function statistics
    int num_comps = 0;
    int num_steps = 0;
    CCTK_REAL num_active_mem_points = 0;
    CCTK_REAL num_owned_mem_points  = 0;
    CCTK_REAL num_total_mem_points  = 0;
    CCTK_REAL num_active_cpu_points = 0;
    CCTK_REAL num_owned_cpu_points  = 0;
    CCTK_REAL num_total_cpu_points  = 0;
    for (int m=0; m<maps; ++m) {
      gh const * const hh = vhh.AT(m);
      dh const * const dd = vdd.AT(m);
      for (int ml=0; ml<mglevels; ++ml) {
        for (int rl=0; rl<reflevels; ++rl) {
          int const trf = timereffacts.AT(rl);
          if (m==0 and ml==0) {
            num_steps += trf;
          }
          for (int c=0; c<hh->components(rl); ++c) {
            ++ num_comps;
            dh::dboxes const & b = dd->boxes.AT(ml).AT(rl).AT(c);
            num_active_mem_points += num_gfs * b.active.size();
            num_owned_mem_points  += num_gfs * b.owned.size();
            num_total_mem_points  += num_gfs * b.exterior.size();
            num_active_cpu_points += trf * b.active.size();
            num_owned_cpu_points  += trf * b.owned.size();
            num_total_cpu_points  += trf * b.exterior.size();
          }
        }
      }
    }
    num_active_cpu_points /= num_steps * delta_time;
    num_owned_cpu_points  /= num_steps * delta_time;
    num_total_cpu_points  /= num_steps * delta_time;
    
    // Output
    CCTK_VInfo (CCTK_THORNSTRING,
                "Grid structure statistics:");
    CCTK_VInfo (CCTK_THORNSTRING,
                "GF: rhs: %.0fk active, %.0fk owned (+%.0f%%), %.0fk total (+%.0f%%), %.3g steps/time",
                double (num_active_cpu_points / 1000),
                double (num_owned_cpu_points / 1000),
                double (num_owned_cpu_points / num_active_cpu_points * 100 - 100),
                double (num_total_cpu_points / 1000),
                double (num_total_cpu_points / num_owned_cpu_points * 100 - 100),
                double (num_steps / delta_time));
    CCTK_VInfo (CCTK_THORNSTRING,
                "GF: vars: %d, pts: %.0fM active, %.0fM owned (+%.0f%%), %.0fM total (+%.0f%%), %.1f comp/proc",
                num_gfs,
                double (num_active_mem_points / 1000000),
                double (num_owned_mem_points / 1000000),
                double (num_owned_mem_points / num_active_mem_points * 100 - 100),
                double (num_total_mem_points / 1000000),
                double (num_total_mem_points / num_owned_mem_points * 100 - 100),
                double (num_comps / (reflevels * CCTK_nProcs (cctkGH))));
    CCTK_VInfo (CCTK_THORNSTRING,
                "GA: vars: %d, pts: %.0fM active, %.0fM total (+%.0f%%)",
                num_arrays,
                double (num_active_array_points / 1000000),
                double (num_total_array_points / 1000000),
                double (num_total_array_points / num_active_array_points * 100 - 100));
    
    // After this, we will begin to allocate memory for the grid
    // structure.  If we run out of memory, ensure that this output
    // still makes it to disk.
    fflush (stdout);
    fflush (stderr);
  }
  
  
  
  void
  SplitRegions_AlongZ (cGH const * const cctkGH,
                       vector<region_t> & superregs,
                       vector<region_t> & regs)
  {
    assert (regs.empty());
    SplitRegions_AlongDir (cctkGH, superregs, regs, dim - 1);
  }
  
  
  
  void
  SplitRegions_AlongDir (cGH const * const cctkGH,
                         vector<region_t> & superregs,
                         vector<region_t> & regs,
                         int const dir)
  {
    assert (regs.empty());
    
    // Something to do?
    if (superregs.size() == 0) {
      return;
    }
    
    const int nprocs = CCTK_nProcs (cctkGH);
    
    assert (superregs.size() == 1);
    regs = superregs;
    
    if (nprocs == 1) {
      regs.at(0).processor = 0;
      pseudoregion_t const preg (regs.at(0).extent, regs.at(0).processor);
      assert (superregs.at(0).processors == NULL);
      superregs.at(0).processors = new ipfulltree (preg);
      return;
    }
    
    assert (dir>=0 and dir<dim);
    
    region_t const & reg0 = regs.at(0);
    const ivect  rstr0  = reg0.extent.stride();
    const ivect  rlb0   = reg0.extent.lower();
    const ivect  rub0   = reg0.extent.upper() + rstr0;
    const b2vect obnd0 = reg0.outer_boundaries;
    
    regs.resize (nprocs);
    vector<int> bounds (nprocs+1);
    vector<ipfulltree *> subtrees (nprocs);
    for (int c=0; c<nprocs; ++c) {
      ivect cstr = rstr0;
      ivect clb  = rlb0;
      ivect cub  = rub0;
      const int glonpz = (rub0[dir] - rlb0[dir]) / cstr[dir];
      const int locnpz = (glonpz + nprocs - 1) / nprocs;
      const int zstep  = locnpz * cstr[dir];
      clb[dir] = rlb0[dir] + zstep *  c;
      cub[dir] = rlb0[dir] + zstep * (c+1);
      if (clb[dir] > rub0[dir]) clb[dir] = rub0[dir];
      if (cub[dir] > rub0[dir]) cub[dir] = rub0[dir];
      assert (clb[dir] <= cub[dir]);
      assert (cub[dir] <= rub0[dir]);
      region_t & reg = regs.at(c);
      ibbox    & ext  = reg.extent;
      b2vect   & obnd = reg.outer_boundaries;
      int      & proc = reg.processor;
      ext = ibbox(clb, cub-cstr, cstr);
      obnd = obnd0;
      obnd[0][dir] &= clb[dir] == rlb0[dir];
      obnd[1][dir] &= cub[dir] == rub0[dir];
      proc = c;
      pseudoregion_t const preg (reg.extent, c);
      bounds.at(c) = reg.extent.lower()[dir];
      subtrees.at(c) = new ipfulltree (preg);
    }
    bounds.at(nprocs) = rub0[dir] + rstr0[dir];
    
    assert (superregs.at(0).processors == NULL);
    superregs.at(0).processors = new ipfulltree (dir, bounds, subtrees);
    
    for (int c=0; c<(int)regs.size(); ++c) {
      assert (regs.at(c).processor == c);
    }
  }
  
  
  
  void
  SplitRegions_Automatic (cGH const * const cctkGH,
                          vector<region_t> & superregs,
                          vector<region_t> & regs)
  {
    assert (regs.empty());
    vector<vector<region_t> > superregss (1, superregs);
    vector<vector<region_t> > regss (1);
    SplitRegionsMaps_Automatic (cctkGH, superregss, regss);
    assert (superregss.size() == 1);
    superregs = superregss.at(0);
    assert (regss.size() == 1);
    regs = regss.at(0);
  }
  
  
  
  static void
  SplitRegions_AsSpecified (cGH const * const cctkGH,
                            vector<region_t> & superregs,
                            vector<region_t> & regs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (regs.empty());
    
    // Something to do?
    if (superregs.size() == 0) {
      return;
    }
    
    const int nprocs = CCTK_nProcs (cctkGH);
    
    assert (superregs.size() == 1);
    
    region_t const & reg0 = superregs.at(0);
    const ivect rstr0  = reg0.extent.stride();
    const ivect rlb0   = reg0.extent.lower();
    const ivect rub0   = reg0.extent.upper() + rstr0;
    const b2vect obnd0 = reg0.outer_boundaries;
    
    const ivect nprocs_dir
      (processor_topology_3d_x,
       processor_topology_3d_y, 
       processor_topology_3d_z);
    assert (all (nprocs_dir > 0));
    if (prod (nprocs_dir) != nprocs) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The specified processor topology [%d,%d,%d] requires %d processors, but there are %d processors",
                  nprocs_dir[0], nprocs_dir[1], nprocs_dir[2],
                  prod (nprocs_dir),
                  nprocs);
    }
    assert (prod (nprocs_dir) == nprocs);
    
    regs.resize (nprocs);
    const ivect cstr = rstr0;
    const ivect glonp = (rub0 - rlb0) / cstr;
//     const ivect locnp = (glonp + nprocs_dir - 1) / nprocs_dir;
    const ivect locnp = glonp / nprocs_dir;
    const ivect rem = glonp % nprocs_dir;
    const ivect step = locnp * cstr;
    assert (dim==3);
    
    vector<int> boundsz(nprocs_dir[2]);
    vector<ipfulltree *> subtreesz(nprocs_dir[2]+1);
    for (int k=0; k<nprocs_dir[2]; ++k) {
      vector<int> boundsy(nprocs_dir[2]);
      vector<ipfulltree *> subtreesy(nprocs_dir[2]+1);
      for (int j=0; j<nprocs_dir[1]; ++j) {
        vector<int> boundsx(nprocs_dir[2]);
        vector<ipfulltree *> subtreesx(nprocs_dir[2]+1);
	for (int i=0; i<nprocs_dir[0]; ++i) {
          
	  const int c = i + nprocs_dir[0] * (j + nprocs_dir[1] * k);
	  const ivect ipos (i, j, k);
	  ivect clb = rlb0 + step *  ipos;
	  ivect cub = rlb0 + step * (ipos+1);
// 	  clb = min (clb, rub);
// 	  cub = min (cub, rub);
          for (int d=0; d<dim; ++d) {
            if (ipos[d]<rem[d]) {
              clb[d] += cstr[d] * ipos[d];
              cub[d] += cstr[d] * (ipos[d]+1);
            } else {
              clb[d] += cstr[d] * rem[d];
              cub[d] += cstr[d] * rem[d];
            }
          }
	  assert (all (clb >= 0));
	  assert (all (clb <= cub));
	  assert (all (cub <= rub0));
          assert (all (not (ipos==0) or clb==rlb0));
          assert (all (not (ipos==nprocs_dir-1) or cub==rub0));
          region_t & reg  = regs.at(c);
          ibbox    & ext  = reg.extent;
          b2vect   & obnd = reg.outer_boundaries;
          int      & proc = reg.processor;
          ext = ibbox(clb, cub-cstr, cstr);
	  obnd = obnd0;
          obnd[0] &= clb == rlb0;
          obnd[1] &= cub == rub0;
          proc = c;
          
          pseudoregion_t preg (reg.extent, c);
          subtreesx.at(i) = new ipfulltree (preg);
	}
        boundsx.at(nprocs_dir[0]) = rub0[0] + rstr0[0];
        subtreesy.at(j) = new ipfulltree (0, boundsx, subtreesx);
      }
      boundsy.at(nprocs_dir[1]) = rub0[1] + rstr0[1];
      subtreesz.at(k) = new ipfulltree (1, boundsy, subtreesy);
    }
    boundsz.at(nprocs_dir[2]) = rub0[2] + rstr0[2];
    
    assert (superregs.at(0).processors == NULL);
    superregs.at(0).processors = new ipfulltree (2, boundsz, subtreesz);
    
    for (int c=0; c<(int)regs.size(); ++c) {
      assert (regs.at(c).processor == c);
    }
  }
  
  
  
  // TODO: this routine should go into CarpetRegrid (except maybe
  // SplitRegions_AlongZ for grid arrays)
  void
  SplitRegionsMaps (cGH const * const cctkGH,
                    vector<vector<region_t> > & superregss,
                    vector<vector<region_t> > & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (regss.size() == superregss.size());
    for (size_t m=0; m<regss.size(); ++m) {
      assert (regss.AT(m).empty());
    }
    if (CCTK_EQUALS (processor_topology, "along-z")) {
      assert (0);
//       SplitRegionsMaps_AlongZ (cctkGH, superregss, regss);
    } else if (CCTK_EQUALS (processor_topology, "along-dir")) {
      assert (0);
//       SplitRegionsMaps_AlongDir (cctkGH, superregss, regss, split_direction);
    } else if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegionsMaps_Automatic (cctkGH, superregss, regss);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      assert (0);
//       SplitRegionsMaps_AsSpecified (cctkGH, superregss, regss);
    } else {
      assert (0);
    }
  }
  
  
  
  static void
  SplitRegionsMaps_Automatic_Recursively (bvect    const   & dims,
                                          int      const     firstproc,
                                          int      const     nprocs,
                                          region_t         & superreg,
                                          vector<region_t> & newregs)
  {
    if (DEBUG) cout << "SRMAR enter" << endl;
    // Check preconditions
    assert (nprocs >= 1);
    
    // Remember old size
    size_t const oldsize = newregs.size();
    
    // Are we done?
    if (all(dims)) {
      if (DEBUG) cout << "SRMAR bottom" << endl;
      
      // Check precondition
      assert (nprocs == 1);
      
      // Return argument
      region_t newreg = superreg;
      newreg.processor = firstproc;
      newregs.push_back (newreg);
      pseudoregion_t const preg (newreg.extent, newreg.processor);
      superreg.processors = new ipfulltree (preg);
      
      // Check postcondition
      assert (newregs.size() == oldsize + nprocs);
      
      if (DEBUG) cout << "SRMAR exit" << endl;
      return;
    }
    
    // Special case
    if (superreg.extent.empty()) {
      if (DEBUG) cout << "SRMAR empty" << endl;
      
      // Create a new region
      region_t newreg (superreg);
      newreg.outer_boundaries = b2vect(false);
      if (DEBUG) cout << "SRMAR newreg " << newreg << endl;
      
      // Store
      for (int p=0; p<nprocs; ++p) {
        newreg.processor = firstproc + p;
        newregs.push_back (newreg);
      }
      superreg.processors = new ipfulltree ();
      
      // Check postcondition
      assert (newregs.size() == oldsize + nprocs);
      
      if (DEBUG) cout << "SRMAR exit" << endl;
      return;
    }
    
    // Calculate cost factors
    rvect rcost = cost (superreg);
    CCTK_REAL const rfact = pow (nprocs / prod(rcost), CCTK_REAL(1)/dim);
    rcost *= rfact;
    assert (abs (prod (rcost) - nprocs) <= 1.0e-6);
    if (DEBUG) cout << "SRMA shapes " << rcost << endl;
    
    // Choose a direction
    int mydim = -1;
    int alldims = 0;
    CCTK_REAL mycost = 0;
    CCTK_REAL totalcost = 1;
    // Prefer to split in the z direction
    for (int d=dim-1; d>=0; --d) {
      if (not dims[d]) {
        ++ alldims;
        CCTK_REAL const thiscost = rcost[d];
        if (thiscost >= 0.999999 * mycost) {
          mydim = d;
          mycost = thiscost;
        }
        totalcost *= thiscost;
      }
    }
    assert (mydim>=0 and mydim<dim);
    assert (mycost>=0);
    if (DEBUG) cout << "SRMAR mydim " << mydim << endl;
    if (DEBUG) cout << "SRMAR mycost " << mycost << endl;
    
    // Mark this direction as done
    assert (not dims[mydim]);
    bvect const newdims = dims.replace(mydim, true);
    
    // Choose a number of slices for this direction
    int const nslices
      = min (nprocs,
             int (floor (mycost * pow(nprocs / totalcost,
                                      CCTK_REAL(1) / alldims) +
                         CCTK_REAL(0.5))));
    assert (nslices <= nprocs);
    if (DEBUG) cout << "SRMAR " << mydim << " nprocs " << nprocs << endl;
    if (DEBUG) cout << "SRMAR " << mydim << " nslices " << nslices << endl;
    
    // Split the remaining processors
    vector<int> mynprocs(nslices);
    int const mynprocs_base = nprocs / nslices;
    int const mynprocs_left = nprocs - nslices * mynprocs_base;
    int sum_mynprocs = 0;
    for (int n=0; n<nslices; ++n) {
      mynprocs.at(n) = mynprocs_base + int (n < mynprocs_left);
      sum_mynprocs += mynprocs.at(n);
    }
    assert (sum_mynprocs == nprocs);
    if (DEBUG) cout << "SRMAR " << mydim << " mynprocs " << mynprocs << endl;
    
    // Split the region
    vector<int> mynpoints(nslices);
    int const npoints =
      (superreg.extent.shape() / superreg.extent.stride())[mydim];
    
    // Keep track of how many points and processors we have left to
    // distribute
    int npoints_left = npoints;
    int nprocs_left  = nprocs;
    for (int n=0; n<nslices; ++n) {
      mynpoints.at(n) = int (floor (CCTK_REAL(1) * npoints_left * mynprocs.at(n)
                                    / nprocs_left + CCTK_REAL(0.5)));
      assert (mynprocs .at(n) > 0);
      npoints_left -= mynpoints.at(n);
      nprocs_left  -= mynprocs.at(n);
      assert (npoints_left >= 0);
      assert (nprocs_left  >= 0);
    }
    assert (npoints_left == 0);
    assert (nprocs_left  == 0);
    if (DEBUG) cout << "SRMAR " << mydim << " mynpoints " << mynpoints << endl;
    
    // Create the regions and recurse
    if (DEBUG) cout << "SRMAR " << mydim << ": create bboxes" << endl;
    ivect lo = superreg.extent.lower();
    ivect up = superreg.extent.upper();
    ivect const str = superreg.extent.stride();
    int p = firstproc;
    vector<int> bounds (nslices+1);
    vector<ipfulltree *> subtrees (nslices);
    for (int n=0; n<nslices; ++n) {
      if (DEBUG) cout << "SRMAR " << mydim << " n " << n << endl;
      
      // Create a new region
      up[mydim] = lo[mydim] + (mynpoints.at(n) - 1) * str[mydim];
      b2vect newob (superreg.outer_boundaries);
      if (n > 0) {
        newob[0][mydim] = false;
      }
      if (n < nslices-1) {
        up[mydim] = lo[mydim] + (mynpoints.at(n) - 1) * str[mydim];
        newob[1][mydim] = false;
      }
      region_t newreg (superreg);
      newreg.extent = ibbox(lo, up, str);
      newreg.outer_boundaries = newob;
      if (DEBUG) cout << "SRMAR " << mydim << " newreg " << newreg << endl;
      
      // Recurse
      bounds.at(n) = lo[mydim];
      SplitRegionsMaps_Automatic_Recursively
        (newdims, p, mynprocs.at(n), newreg, newregs);
      if (DEBUG) cout << "SRMAR " << mydim << " newregs " << newregs << endl;
      
      assert (newreg.processors != NULL);
      subtrees.at(n) = newreg.processors;
      newreg.processors = NULL;
      newreg.processor = p;
      
      // Next slice
      lo[mydim] = up[mydim] + str[mydim];
      p += mynprocs.at(n);
    }
    assert (up[mydim] == superreg.extent.upper()[mydim]);
    assert (p == firstproc + nprocs);
    bounds.at(nslices) = up[mydim] + str[mydim];
    
    // Create tree
    assert (superreg.processors == NULL);
    superreg.processors = new ipfulltree (mydim, bounds, subtrees);
    
    // Check postcondition
    assert (newregs.size() == oldsize + nprocs);
    
    if (DEBUG) cout << "SRMAR exit" << endl;
  }
  
  
  
  void
  SplitRegionsMaps_Automatic (cGH const * const cctkGH,
                              vector<vector<region_t> > & superregss,
                              vector<vector<region_t> > & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (DEBUG) cout << "SRMA enter" << endl;
    
    int const nmaps = superregss.size();
    int minmap = 1000000000;
    for (int m=0; m<nmaps; ++m) {
      for (int r=0; r<int(superregss.at(m).size()); ++r) {
        minmap = min (minmap, superregss.at(m).at(r).map);
      }
    }
    int nregs = 0;
    for (int m=0; m<nmaps; ++m) {
      nregs += superregss.at(m).size();
    }
    if (DEBUG) cout << "SRMA nregs " << nregs << endl;
    
    // Something to do?
    if (nregs == 0) {
      return;
    }
    
    // Collect slices
    vector<region_t> superregs;
    {
      for (int m=0; m<nmaps; ++m) {
        combine_regions (superregss.at(m), superregs);
      }
      nregs = superregs.size();
      
      // If the last region was removed, add a new empty region again.
      // A set of regions (corresponding to a refinement level or a
      // grid array) cannot be empty.
      if (nregs == 0) {
        assert (nmaps == 1);    // we should only be here for grid
                                // arrays
        region_t reg;
        reg.extent           = ibbox (ivect (0), ivect (-1), ivect (1));
        reg.outer_boundaries = b2vect (bvect (true), bvect (true));
        reg.map              = 0;
        reg.processor        = -1;
        superregs.push_back (reg);
        nregs = superregs.size();
      }
    }
    
    int const real_nprocs = CCTK_nProcs (cctkGH);
    if (DEBUG) cout << "SRMA real_nprocs " << real_nprocs << endl;
    
    // Deactivate some processors if there are too many
    int nprocs;
    if (min_points_per_proc < 0) {
      nprocs = real_nprocs;
    } else {
      CCTK_REAL mycost = 0;
      for (int r=0; r<nregs; ++r) {
        mycost += prod (cost (superregs.at(r)));
      }
      int const goodnprocs = int (floor (mycost / min_points_per_proc));
      nprocs = max (1, min (real_nprocs, goodnprocs));
    }
    if (DEBUG) cout << "SRMA nprocs " << nprocs << endl;
    
    // ncomps: number of components per processor
    int const ncomps = (nregs + nprocs - 1) / nprocs;
    if (DEBUG) cout << "SRMA ncomps " << ncomps << endl;
    assert (ncomps > 0);
    
    int const newnregs = nprocs * ncomps;
    if (DEBUG) cout << "SRMA newnregs " << newnregs << endl;
    
    if (DEBUG) cout << "SRMA: distributing processors to regions" << endl;
    vector<CCTK_REAL> mycosts(nregs);
    for (int r=0; r<nregs; ++r) {
      mycosts.at(r) = prod (cost (superregs.at(r)));
    }
    int nregs_left = newnregs;
    vector<int> mynprocs(nregs);
    for (int r=0; r<nregs; ++r) {
      mynprocs.at(r) = 1;
      -- nregs_left;
    }
#warning "TODO: split regions if necessary"
    while (nregs_left > 0) {
      if (DEBUG) cout << "SRMA nregs_left " << nregs_left << endl;
      int maxr = -1;
      CCTK_REAL maxratio = -1;
      for (int r=0; r<nregs; ++r) {
        CCTK_REAL const ratio = mycosts.at(r) / mynprocs.at(r);
        if (ratio > maxratio) { maxr=r; maxratio=ratio; }
      }
      assert (maxr>=0 and maxr<nregs);
      ++ mynprocs.at(maxr);
      if (DEBUG) cout << "SRMA maxr " << maxr << endl;
      if (DEBUG) cout << "SRMA mynprocs[maxr] " << mynprocs.at(maxr) << endl;
      -- nregs_left;
    }
    assert (nregs_left == 0);
    int sum_nprocs = 0;
    for (int r=0; r<nregs; ++r) {
      sum_nprocs += mynprocs.at(r);
    }
    assert (sum_nprocs == newnregs);
    if (DEBUG) cout << "SRMA mynprocs " << mynprocs << endl;
    
    if (DEBUG) cout << "SRMA: splitting work units" << endl;
#warning "TODO: rename newregs to regs"
#warning "TODO: rename nregs to nsuperregs"
#warning "TODO: rename newnregs to nregs"
    vector<region_t> newregs;
    newregs.reserve (newnregs);
    for (int r=0, p=0; r<nregs; p+=mynprocs.at(r), ++r) {
      if (DEBUG) cout << "SRMA superreg[" << r << "] " << superregs.at(r) << endl;
      bvect const dims = false;
      SplitRegionsMaps_Automatic_Recursively
        (dims, p, mynprocs.at(r), superregs.at(r), newregs);
    } // for r
    if (DEBUG) cout << "SRMA newregs " << newregs << endl;
    assert (int(newregs.size()) == newnregs);
    
    // Count components per map
    vector<int> myncomps(nmaps, 0);
    vector<int> empty_comps(nmaps, 0);
    for (int r=0; r<newnregs; ++r) {
      int const m = newregs.at(r).map - minmap;
      assert (m>=0 and m<nmaps);
      if (not newregs.at(r).extent.empty()) {
        // Ignore empty regions, which may come from empty grid arrays
        ++ myncomps.at(m);
      } else {
        ++ empty_comps.at(m);
      }
    }
    vector<int> mynregs(nmaps, 0);
    for (int r=0; r<nregs; ++r) {
      int const m = superregs.at(r).map - minmap;
      assert (m>=0 and m<nmaps);
      ++ mynregs.at(m);
    }
    
    // Convert virtual to real processors
    for (int r=0; r<newnregs; ++r) {
      newregs.at(r).processor /= ncomps;
      assert (newregs.at(r).processor >= 0 and
              newregs.at(r).processor < nprocs);
    }
    {
      vector<int> tmpncomps(nmaps, 0);
      for (int r=0; r<nregs; ++r) {
        int const m = superregs.at(r).map - minmap;
        ipfulltree * const regf = superregs.at(r).processors;
        assert (regf != NULL);
        for (ipfulltree::iterator fti (* regf); not fti.done(); ++ fti) {
          pseudoregion_t & preg = (* fti).payload();
          preg.component = tmpncomps.at(m)++;
        }
      }
      for (int m=0; m<nmaps; ++m) {
        assert (tmpncomps.at(m) == myncomps.at(m));
      }
    }
    
    // Distribute regions
    // Allocate regions
    assert ((int)regss.size() == nmaps);
    for (int m=0; m<nmaps; ++m) {
      assert (regss.at(m).empty());
      regss.at(m).reserve (myncomps.at(m));
      superregss.at(m).clear();
      superregss.at(m).reserve (mynregs.at(m));
    }
    // Assign regions
    for (int r=0; r<newnregs; ++r) {
      int const m = newregs.at(r).map - minmap;
      assert (m>=0 and m<nmaps);
      regss.at(m).push_back (newregs.at(r));
    }
    for (int r=0; r<nregs; ++r) {
      int const m = superregs.at(r).map - minmap;
      assert (m>=0 and m<nmaps);
      superregss.at(m).push_back (superregs.at(r));
    }
    // Output regions
    if (DEBUG) {
      cout << "SRMA superregss " << superregss << endl;
      cout << "SRMA regss " << regss << endl;
    }
    // Check sizes
    for (int m=0; m<nmaps; ++m) {
      assert (int(regss.at(m).size()) == myncomps.at(m) + empty_comps.at(m));
      assert (int(superregss.at(m).size()) == mynregs.at(m));
    }
    
    if (DEBUG) cout << "SRMA exit" << endl;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  static void
  MakeMultigridBoxes (cGH const * const cctkGH,
                      int const m,
                      ibbox    const & base,
                      region_t const & reg,
                      vector<region_t> & regs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    regs.resize (mglevels, reg);
    if (mglevels > 1) {
      // boundary offsets
      jjvect nboundaryzones, is_internal, is_staggered, shiftout;
      if (domain_from_multipatch and
          CCTK_IsFunctionAliased ("MultiPatch_GetBoundarySpecification"))
      {
        check (not MultiPatch_GetBoundarySpecification
               (m, 2*dim, &nboundaryzones[0][0], &is_internal[0][0],
                &is_staggered[0][0], &shiftout[0][0]));
      } else {
        check (not GetBoundarySpecification
               (2*dim, &nboundaryzones[0][0], &is_internal[0][0],
                &is_staggered[0][0], &shiftout[0][0]));
      }
      // (distance in grid points between the exterior and the physical boundary)
      iivect offset;
      for (int d=0; d<dim; ++d) {
        for (int f=0; f<2; ++f) {
          assert (not is_staggered[d][f]);
          offset[d][f] = (+ (is_internal[d][f] ? 0 : nboundaryzones[d][f] - 1)
                          + shiftout[d][f]);
        }
      }
      vector<ibbox> bases(mglevels);
      bases.at(0) = base;
      for (int ml=1; ml<mglevels; ++ml) {
        // next finer base
        ivect const fbaselo = bases.at(ml-1).lower();
        ivect const fbasehi = bases.at(ml-1).upper();
        ivect const fbasestr = bases.at(ml-1).stride();
        // this base
        ivect const basestr = fbasestr * mgfact;
        ivect const baselo = fbaselo + (xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0]) * fbasestr;
        ivect const basehi = fbasehi + (xpose(offset)[1] - ivect(mgfact) * xpose(offset)[1]) * fbasestr;
        ivect const baselo1 = baselo;
        ivect const basehi1 = baselo1 + (basehi - baselo1) / basestr * basestr;
        bases.at(ml) = ibbox(baselo1, basehi1, basestr);
        // next finer grid
        ivect const flo = regs.at(ml-1).extent.lower();
        ivect const fhi = regs.at(ml-1).extent.upper();
        ivect const fstr = regs.at(ml-1).extent.stride();
        // this grid
        ivect const str = fstr * mgfact;
        ivect const lo = flo + either (reg.outer_boundaries[0],   (xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0]) * fstr, ivect(0));
        ivect const hi = fhi + either (reg.outer_boundaries[1], - (xpose(offset)[1] - ivect(mgfact) * xpose(offset)[1]) * fstr, ivect(0));
        ivect const lo1 = baselo1 + (lo - baselo1 + str - 1) / str * str;
        ivect const hi1 = lo1 + (hi - lo1) / str * str;
        regs.at(ml).extent = ibbox(lo1, hi1, str);
      }
    }
  }
  
  void
  MakeMultigridBoxes (cGH const * const cctkGH,
                      int const m,
                      vector<vector<region_t> > const & regss,
                      vector<vector<vector<region_t> > > & regsss)
  {
    regsss.resize (mglevels);
    for (int ml=0; ml<mglevels; ++ml) {
      regsss.at(ml).resize (regss.size());
      for (int rl=0; rl<(int)regss.size(); ++rl) {
        regsss.at(ml).at(rl).resize (regss.at(rl).size());
      }
    }
    
    for (int rl=0; rl<(int)regss.size(); ++rl) {
      
      ibbox base;
      for (int c=0; c<(int)regss.at(rl).size(); ++c) {
        base = base.expanded_containing(regss.at(rl).at(c).extent);
      }
      
      for (int c=0; c<(int)regss.at(rl).size(); ++c) {
        vector<region_t> mg_regs;
        MakeMultigridBoxes (cctkGH, m, base, regss.at(rl).at(c), mg_regs);
        assert ((int)mg_regs.size() == mglevels);
        for (int ml=0; ml<mglevels; ++ml) {
          regsss.at(ml).at(rl).at(c) = mg_regs.at(ml);
        }
        
      } // for c
      
    } // for rl
  }
  
  void
  MakeMultigridBoxesMaps (cGH const * const cctkGH,
                          vector<vector<vector<region_t> > > const & regsss,
                          vector<vector<vector<vector<region_t> > > > & regssss)
  {
    for (int m = 0; m < maps; ++m) {
      MakeMultigridBoxes (cctkGH, m, regsss.at(m), regssss.at(m));
    } // for m
  }
  
} // namespace Carpet
