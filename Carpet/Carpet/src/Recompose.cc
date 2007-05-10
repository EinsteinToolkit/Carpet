#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
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
    if (reg.extent.empty()) return rvect(0);
    return rvect (reg.extent.shape() / reg.extent.stride());
  }
  
  
  
  static void
  SplitRegionsMaps_Automatic_Recursively (bvect    const & dims,
                                          int      const   nprocs,
                                          region_t const & reg,
                                          vector<region_t> & newregs);
  static void
  SplitRegions_AsSpecified (cGH const * cctkGH,
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
      for (int rl=0; rl<(int)regsss.at(0).size(); ++rl) {
        // No empty levels
        assert (regsss.at(ml).at(rl).size() > 0);
        for (int c=0; c<(int)regsss.at(ml).at(rl).size(); ++c) {
          // Check sizes
          // Do allow processors with zero grid points
//           assert (all(regsss.at(rl).at(c).at(ml).extent.lower() <= regsss.at(rl).at(c).at(ml).extent.upper()));
          // Check strides
          const ivect str
            = (maxspacereflevelfact / spacereffacts.at(rl) * ipow(mgfact, ml));
          assert (all(regsss.at(ml).at(rl).at(c).extent.stride() % str == 0));
          // Check alignments
          assert (all(regsss.at(ml).at(rl).at(c).extent.lower() % str == 0));
          assert (all(regsss.at(ml).at(rl).at(c).extent.upper() % str == 0));
        }
      }
    }
  }
  
  
  
  bool
  Regrid (cGH const * const cctkGH,
          bool const force_recompose)
  {
    DECLARE_CCTK_PARAMETERS;
    
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
      if (maxreflevels > 1 and not didtell) {
	CCTK_WARN (2, "The regridding routine Carpet_RegridMaps has not been provided.  Regridding will be performed in singlemap mode instead of level mode.");
	didtell = true;
      }
    }
    
    
    
    bool did_change = false;
    
    if (not use_regridmaps) {
      
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        
        gh::mregs regsss = vhh.at(map)->regions;
        
        // Check whether to recompose
        CCTK_INT const do_recompose =
          Carpet_Regrid (cctkGH, & regsss, force_recompose);
        assert (do_recompose >= 0);
        did_change = did_change or do_recompose;
        
        if (do_recompose) {
          RegridMap (cctkGH, map, regsss);
        }
        
      } END_MAP_LOOP;
      
    } else {
      
      vector<gh::mregs> regssss (maps);
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        regssss.at(map) = vhh.at(map)->regions;
      } END_MAP_LOOP;
      
      // Check whether to recompose
      CCTK_INT const do_recompose =
        Carpet_RegridMaps (cctkGH, & regssss, force_recompose);
      assert (do_recompose >= 0);
      did_change = did_change or do_recompose;
      
      if (do_recompose) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          gh::mregs const & regsss = regssss.at(map);
          RegridMap (cctkGH, map, regsss);
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
             gh::mregs const & regsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Regridding map %d...", m);
    
    // Check the regions
    CheckRegions (regsss);
    // TODO: check also that the current and all coarser levels did
    // not change
    
    // Regrid
    vhh.at(m)->regrid (regsss);
    
    // Write grid structure to file
    OutputGridStructure (cctkGH, m, regsss);
    
    if (verbose) OutputGrids (cctkGH, m, * vhh.at(m), * vdd.at(m));
    
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
               << ")"
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
    if (strcmp(grid_structure_filename, "") == 0) return;
    
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
	file << "# grid structure" << endl
	     << "# format: map reflevel component mglevel   processor bounding-box is-outer-boundary" << endl;
	assert (file.good());
      }
    }
    if (not file.is_open()) {
      file.open (filename, ios::app);
      assert (file.good());
    }
    
    file << "iteration " << cctkGH->cctk_iteration << endl;
    file << "maps " << maps << endl;
    file << m << " mglevels " << regsss.size() << endl;
    for (int ml=0; ml<(int)regsss.size(); ++ml) {
      file << m << " " << ml << " reflevels " << regsss.at(ml).size() << endl;
      for (int rl=0; rl<(int)regsss.at(ml).size(); ++rl) {
        file << m << " " << ml << " " << rl << " components " << regsss.at(ml).at(rl).size() << endl;
        for (int c=0; c<(int)regsss.at(ml).at(rl).size(); ++c) {
          file << m << " " << ml << " " << rl << " " << c << "   "
               << regsss.at(ml).at(rl).at(c).processor << " "
               << regsss.at(ml).at(rl).at(c).extent
               << regsss.at(ml).at(rl).at(c).outer_boundaries << endl;
        }
      }
    }
    file << endl;
    
    file.close();
    assert (file.good());
  }
  
  
  
  // TODO: this routine should go into CarpetRegrid (except maybe
  // SplitRegions_AlongZ for grid arrays)
  void
  SplitRegions (cGH const * const cctkGH,
                vector<region_t> & regs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (processor_topology, "along-z")) {
      SplitRegions_AlongZ (cctkGH, regs);
    } else if (CCTK_EQUALS (processor_topology, "along-dir")) {
      SplitRegions_AlongDir (cctkGH, regs, split_direction);
    } else if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegions_Automatic (cctkGH, regs);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      SplitRegions_AsSpecified (cctkGH, regs);
    } else {
      assert (0);
    }
  }
  
  
  
  void
  SplitRegions_AlongZ (cGH const * const cctkGH,
                       vector<region_t> & regs)
  {
    SplitRegions_AlongDir (cctkGH, regs, dim - 1);
  }
  
  
  
  void
  SplitRegions_AlongDir (cGH const * const cctkGH,
                         vector<region_t> & regs,
                         int const dir)
  {
    // Something to do?
    if (regs.size() == 0) {
      return;
    }
    
    const int nprocs = CCTK_nProcs (cctkGH);
    
    assert (regs.size() == 1);
    
    if (nprocs == 1) {
      regs.at(0).processor = 0;
      return;
    }
    
    assert (dir>=0 and dir<dim);
    
    region_t const & reg0 = regs.at(0);
    const ivect  rstr0  = reg0.extent.stride();
    const ivect  rlb0   = reg0.extent.lower();
    const ivect  rub0   = reg0.extent.upper() + rstr0;
    const b2vect obnd0 = reg0.outer_boundaries;
    
    regs.resize (nprocs);
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
    }
    
    for (int c=0; c<(int)regs.size(); ++c) {
      assert (regs.at(c).processor == c);
    }
  }
  
  
  
  void
  SplitRegions_Automatic (cGH const * const cctkGH,
                          vector<region_t> & regs)
  {
    vector<vector<region_t> > regss (1, regs);
    SplitRegionsMaps_Automatic (cctkGH, regss);
    assert (regss.size() == 1);
    regs = regss.at(0);
  }
  
  
  
  static void
  SplitRegions_AsSpecified (cGH const * const cctkGH,
                            vector<region_t> & regs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Something to do?
    if (regs.size() == 0) {
      return;
    }
    
    const int nprocs = CCTK_nProcs (cctkGH);
    
    assert (regs.size() == 1);
    
    region_t const & reg0 = regs.at(0);
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
    for (int k=0; k<nprocs_dir[2]; ++k) {
      for (int j=0; j<nprocs_dir[1]; ++j) {
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
	}
      }
    }
    
    for (int c=0; c<(int)regs.size(); ++c) {
      assert (regs.at(c).processor == c);
    }
  }
  
  
  
  // TODO: this routine should go into CarpetRegrid (except maybe
  // SplitRegions_AlongZ for grid arrays)
  void
  SplitRegionsMaps (cGH const * const cctkGH,
                    vector<vector<region_t> > & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (processor_topology, "along-z")) {
      assert (0);
//       SplitRegionsMaps_AlongZ (cctkGH, regss);
    } else if (CCTK_EQUALS (processor_topology, "along-dir")) {
      assert (0);
//       SplitRegionsMaps_AlongDir (cctkGH, regss, split_direction);
    } else if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegionsMaps_Automatic (cctkGH, regss);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      assert (0);
//       SplitRegionsMaps_AsSpecified (cctkGH, regss);
    } else {
      assert (0);
    }
  }
  
  
  
  static void
  SplitRegionsMaps_Automatic_Recursively (bvect    const & dims,
                                          int      const   nprocs,
                                          region_t const & reg,
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
      newregs.push_back (reg);
      
      // Check postcondition
      assert (newregs.size() == oldsize + nprocs);
      
      if (DEBUG) cout << "SRMAR exit" << endl;
      return;
    }
    
    // Special case
    if (reg.extent.empty()) {
      if (DEBUG) cout << "SRMAR empty" << endl;
      
      // Create a new region
      region_t newreg (reg);
      newreg.outer_boundaries      = b2vect(false);
      if (DEBUG) cout << "SRMAR newreg " << newreg << endl;
      
      // Store
      for (int p=0; p<nprocs; ++p) {
        newreg.processor = reg.processor + p;
        newregs.push_back (newreg);
      }
      
      // Check postcondition
      assert (newregs.size() == oldsize + nprocs);
      
      if (DEBUG) cout << "SRMAR exit" << endl;
      return;
    }
    
    // Calculate cost factors
    rvect rcost = cost (reg);
    CCTK_REAL const rfact = pow (nprocs / prod(rcost), CCTK_REAL(1)/dim);
    rcost *= rfact;
    assert (abs (prod (rcost) - nprocs) < 1.0e-6);
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
    int const npoints = (reg.extent.shape() / reg.extent.stride())[mydim];
    
    // Keep track of how many points and processors we have left to
    // distribute
    int npoints_left = npoints;
    int nprocs_left  = nprocs;
    for (int n=0; n<nslices; ++n) {
      mynpoints.at(n) = int (floor (CCTK_REAL(1) * npoints_left * mynprocs.at(n)
                                    / nprocs_left + CCTK_REAL(0.5)));
      assert (mynpoints.at(n) > 0);
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
    ivect lo = reg.extent.lower();
    ivect up = reg.extent.upper();
    ivect const str = reg.extent.stride();
    int p = reg.processor;
    for (int n=0; n<nslices; ++n) {
      if (DEBUG) cout << "SRMAR " << mydim << " n " << n << endl;
      
      // Create a new region
      up[mydim] = lo[mydim] + (mynpoints.at(n) - 1) * str[mydim];
      b2vect newob (reg.outer_boundaries);
      if (n > 0) {
        newob[0][mydim] = false;
      }
      if (n < nslices-1) {
        up[mydim] = lo[mydim] + (mynpoints.at(n) - 1) * str[mydim];
        newob[1][mydim] = false;
      }
      region_t newreg (reg);
      newreg.extent = ibbox(lo, up, str);
      newreg.outer_boundaries      = newob;
      newreg.processor = p;
      if (DEBUG) cout << "SRMAR " << mydim << " newreg " << newreg << endl;
      
      // Recurse
      SplitRegionsMaps_Automatic_Recursively
        (newdims, mynprocs.at(n), newreg, newregs);
      if (DEBUG) cout << "SRMAR " << mydim << " newregs " << newregs << endl;
      
      // Next slice
      lo[mydim] = up[mydim] + str[mydim];
      p += mynprocs.at(n);
    }
    assert (up[mydim] == reg.extent.upper()[mydim]);
    assert (p == reg.processor + nprocs);
    
    // Check postcondition
    assert (newregs.size() == oldsize + nprocs);
    
    if (DEBUG) cout << "SRMAR exit" << endl;
  }
  
  
  
  void
  SplitRegionsMaps_Automatic (cGH const * const cctkGH,
                              vector<vector<region_t> > & regss)
  {
    if (DEBUG) cout << "SRMA enter" << endl;
    
    int const nmaps = regss.size();
    int nregs = 0;
    for (int m=0; m<nmaps; ++m) {
      nregs += regss.at(m).size();
    }
    if (DEBUG) cout << "SRMA nregs " << nregs << endl;
    
    // Something to do?
    if (nregs == 0) {
      return;
    }
    
    // Collect slices
    vector<region_t> regs(nregs);
    {
      int r=0;
      for (int m=0; m<nmaps; ++m) {
        for (int c=0; c<int(regss.at(m).size()); ++c, ++r) {
          regs.at(r) = regss.at(m).at(c);
          regs.at(r).map       = m;
          regs.at(r).processor = -1;
        }
      }
      assert (r == nregs);
    }
    
    const int nprocs = CCTK_nProcs (cctkGH);
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
      mycosts.at(r) = prod (cost (regs.at(r)));
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
    vector<region_t> newregs;
    newregs.reserve (newnregs);
    for (int r=0, p=0; r<nregs; p+=mynprocs.at(r), ++r) {
      regs.at(r).processor = p;
      if (DEBUG) cout << "SRMA reg[" << r << "] " << regs.at(r) << endl;
      bvect const dims = false;
      SplitRegionsMaps_Automatic_Recursively
        (dims, mynprocs.at(r), regs.at(r), newregs);
    } // for r
    if (DEBUG) cout << "SRMA newregs " << newregs << endl;
    assert (int(newregs.size()) == newnregs);
    
    // Convert virtual to real processors
    for (int r=0; r<newnregs; ++r) {
      newregs.at(r).processor /= ncomps;
      assert (newregs.at(r).processor >= 0 and
              newregs.at(r).processor < nprocs);
    }
    
    // Distribute regions
    // Count components per map
    vector<int> myncomps(nmaps);
    for (int r=0; r<newnregs; ++r) {
      int const m = newregs.at(r).map;
      assert (m>=0 and m<nmaps);
      ++ myncomps.at(m);
    }
    // Allocate regions
    for (int m=0; m<nmaps; ++m) {
      regss.at(m).resize (0);
      regss.at(m).reserve (myncomps.at(m));
    }
    // Assign regions
    for (int r=0; r<newnregs; ++r) {
      int const m = newregs.at(r).map;
      assert (m>=0 and m<nmaps);
      regss.at(m).push_back (newregs.at(r));
    }
    // Check sizes
    for (int m=0; m<nmaps; ++m) {
      assert (int(regss.at(m).size()) == myncomps.at(m));
    }
    
    if (DEBUG) cout << "SRMA exit" << endl;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  static void
  MakeMultigridBoxes (cGH const * const cctkGH,
                      ibbox    const & base,
                      region_t const & reg,
                      vector<region_t> & regs)
  {
    regs.resize (mglevels, reg);
    if (mglevels > 1) {
      // boundary offsets
      jjvect nboundaryzones, is_internal, is_staggered, shiftout;
      const int ierr = GetBoundarySpecification
        (2*dim, &nboundaryzones[0][0], &is_internal[0][0],
         &is_staggered[0][0], &shiftout[0][0]);
      assert (not ierr);
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
                      vector<vector<region_t> >  const & regss,
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
        MakeMultigridBoxes (cctkGH, base, regss.at(rl).at(c), mg_regs);
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
      MakeMultigridBoxes (cctkGH, regsss.at(m), regssss.at(m));
    } // for m
  }
  
} // namespace Carpet
