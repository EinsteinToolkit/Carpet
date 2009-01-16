#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <defs.hh>
#include <ggf.hh>

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  //
  // Mode indicators
  //
  
  bool is_meta_mode ()
  {
    if (mglevel==-1) assert (reflevel==-1 and map==-1 and component==-1);
    return mglevel==-1;
  }
  
  bool is_global_mode ()
  {
    if (mglevel==-1 and reflevel!=-1) assert (map==-1 and component==-1);
    return mglevel!=-1 and reflevel==-1;
  }
  
  bool is_level_mode ()
  {
    if (mglevel!=-1 and reflevel!=-1 and map==-1) assert (component==-1);
    return mglevel!=-1 and reflevel!=-1 and map==-1;
  }
  
  bool is_singlemap_mode ()
  {
    return mglevel!=-1 and reflevel!=-1 and map!=-1 and component==-1;
  }
  
  bool is_local_mode ()
  {
    return mglevel!=-1 and reflevel!=-1 and map!=-1 and component!=-1;
//       assert (mglevel>=0 and mglevel<mglevels);
//       assert (reflevel>=0 and reflevel<reflevels);
//       assert (map>=0 and map<maps);
//       assert (vhh.at(map)->local_components(reflevel)==1 or component==-1);
  }
  
  
  
  //
  // Mode setting
  //
  
  // Set global mode
  
  void enter_global_mode (cGH * const cctkGH, int const ml)
  {
    assert (is_meta_mode());
    assert (ml>=0 and ml<mglevels);
    Checkpoint ("Entering global mode");
    
    mglevel = ml;
    mglevelfact = ipow(mgfact, mglevel);
    // TODO: this could also just be "mglevel" instead
    cctkGH->cctk_convlevel = basemglevel + mglevel;
    
    // Set time delta
    cctkGH->cctk_delta_time = delta_time * mglevelfact;
    if (maps == 1) {
      // Set space delta
      for (int d=0; d<dim; ++d) {
        cctkGH->cctk_origin_space[d] = origin_space.at(0).at(mglevel)[d];
        cctkGH->cctk_delta_space[d] = delta_space.at(0)[d] * mglevelfact;
      }
    }
    
    // Set array information
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      cGroup gp;
      check (not CCTK_GroupData (group, &gp));
      if (gp.grouptype != CCTK_GF) {
        
        const int rl = 0;
        const int m = 0;
        const int c = CCTK_MyProc(cctkGH);
        
        const ibbox& baseext = arrdata.at(group).at(m).hh->baseextents.at(ml).at(rl);
	const ibbox& ext = arrdata.at(group).at(m).dd->boxes.at(ml).at(rl).at(c).exterior;
        const b2vect& obnds = arrdata.at(group).at(m).hh->outer_boundaries(rl,c);
        
        ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
          = arrdata.at(group).at(m).dd->ghost_width[0];
        ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
          = baseext.shape() / baseext.stride();
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
          = ext.shape() / ext.stride();
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
          = (ext.lower() - baseext.lower()) / ext.stride();
        ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
          = (ext.upper() - baseext.lower()) / ext.stride();
        if (gp.disttype == CCTK_DISTRIB_CONSTANT) {
          int const d = gp.dim==0 ? 0 : gp.dim-1;
          ivect & gsh = ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh));
          ivect & lsh = ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh));
          ivect & lbnd = ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd));
          ivect & ubnd = ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd));
          gsh[d] = lsh[d];
          lbnd[d] = 0;
          ubnd[d] = lsh[d] - 1;
        }
        for (int d=0; d<dim; ++d) {
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ] = obnds[0][d];
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1] = obnds[1][d];
        }
        groupdata.at(group).info.activetimelevels
          = groupdata.at(group).activetimelevels.at(mglevel).at(0);
        
        for (int d=0; d<dim; ++d) {
          assert (groupdata.at(group).info.lsh[d]>=0);
          assert (groupdata.at(group).info.lsh[d]<=groupdata.at(group).info.gsh[d]);
          assert (groupdata.at(group).info.lbnd[d]>=0);
          assert (groupdata.at(group).info.lbnd[d]<=groupdata.at(group).info.ubnd[d]+1);
          assert (groupdata.at(group).info.ubnd[d]<groupdata.at(group).info.gsh[d]);
          assert (groupdata.at(group).info.lbnd[d] + groupdata.at(group).info.lsh[d] - 1
                  == groupdata.at(group).info.ubnd[d]);
          assert (groupdata.at(group).info.lbnd[d]<=groupdata.at(group).info.ubnd[d]+1);
        }
        
        const int numvars = CCTK_NumVarsInGroupI (group);
        if (numvars>0) {
          const int firstvar = CCTK_FirstVarIndexI (group);
          assert (firstvar>=0);
          const int max_tl = CCTK_MaxTimeLevelsGI (group);
          assert (max_tl>=0);
          const int active_tl = groupdata.at(group).info.activetimelevels;
          assert (active_tl>=0 and active_tl<=max_tl);
          
          assert (arrdata.at(group).at(m).hh->is_local(rl,c));
          
          for (int var=0; var<numvars; ++var) {
            assert (firstvar+var<CCTK_NumVars());
            ggf * const ff = arrdata.at(group).at(m).data.at(var);
            for (int tl=0; tl<max_tl; ++tl) {
              if (ff and tl<active_tl) {
                gdata * const data = (*ff) (tl, rl, c, ml);
                assert (data);
                cctkGH->data[firstvar+var][tl] = data->storage();
              } else {
                cctkGH->data[firstvar+var][tl] = NULL;
              }
            }
          }
        }
        
      } // if grouptype
    } // for group
    
    assert (is_global_mode());
  }
  
  void leave_global_mode (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_global_mode() or is_meta_mode());
    
    if (mglevel == -1) return;  // early return
    
    Checkpoint ("Leaving global mode");

    // Save and unset time delta
    delta_time = cctkGH->cctk_delta_time / mglevelfact;
    cctkGH->cctk_delta_time = 0.0;
    if (maps == 1) {
      // Save and unset space delta
      for (int d=0; d<dim; ++d) {
        origin_space.at(0).at(mglevel)[d] = cctkGH->cctk_origin_space[d];
        delta_space.at(mglevel)[d] = cctkGH->cctk_delta_space[d] / mglevelfact;
        cctkGH->cctk_origin_space[d] = -424242.0;
        cctkGH->cctk_delta_space[d] = -424242.0;
      }
    }
    
    // Set array information
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) != CCTK_GF) {
        
        const int m = 0;
        
//         ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
//           = deadbeef;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
          = arrdata.at(group).at(m).dd->ghost_width[0];
        ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
          = deadbeef;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
          = deadbeef;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
          = -deadbeef;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
          = deadbeef;
        for (int d=0; d<dim; ++d) {
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ] = deadbeef;
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1] = deadbeef;
        }
        groupdata.at(group).info.activetimelevels = deadbeef;
        
        const int numvars = CCTK_NumVarsInGroupI (group);
        if (numvars>0) {
          const int firstvar = CCTK_FirstVarIndexI (group);
          assert (firstvar>=0);
          const int max_tl = CCTK_MaxTimeLevelsGI (group);
          assert (max_tl>=0);
          
          assert (group<(int)arrdata.size());
          for (int var=0; var<numvars; ++var) {
            assert (firstvar+var<CCTK_NumVars());
            for (int tl=0; tl<max_tl; ++tl) {
              cctkGH->data[firstvar+var][tl] = NULL;
            }
          }
        }
        
      } // if grouptype
    } // for group
    
    mglevel = -1;
    mglevelfact = -deadbeef;
    cctkGH->cctk_convlevel = -deadbeef;
    
    assert (is_meta_mode());
  }
  
  
  
  // Set level mode
  
  void enter_level_mode (cGH * const cctkGH, int const rl)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_global_mode());
    assert (rl>=0 and rl<reflevels);
    Checkpoint ("Entering level mode");
    
    reflevel = rl;
    timereflevelfact = timereffacts.at (reflevel);
    spacereflevelfact = spacereffacts.at (reflevel);
    ivect::ref(cctkGH->cctk_levfac) = spacereflevelfact;
    cctkGH->cctk_timefac = timereflevelfact;
    
    // Set number of time levels
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        groupdata.at(group).info.activetimelevels
          = groupdata.at(group).activetimelevels.at(mglevel).at(reflevel);
      }
    }
    
    // Set current time
    assert (mglevel>=0 and mglevel<(int)leveltimes.size());
    assert (reflevel>=0 and reflevel<(int)leveltimes.at(mglevel).size());
    if (not adaptive_stepsize) {
      cctkGH->cctk_time = leveltimes.at(mglevel).at(reflevel);
    } else {
      leveltimes.at(mglevel).at(reflevel) = cctkGH->cctk_time;
    }
    
    assert (is_level_mode());
  }
  
  void leave_level_mode (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_level_mode() or is_global_mode());
    
    if (reflevel == -1) return; // early return
    
    Checkpoint ("Leaving level mode");

    // Save and unset current time
    assert (mglevel>=0 and mglevel<(int)leveltimes.size());
    assert (reflevel>=0 and reflevel<(int)leveltimes.at(mglevel).size());
    leveltimes.at(mglevel).at(reflevel) = cctkGH->cctk_time;
    if (not adaptive_stepsize) {
      cctkGH->cctk_time = global_time;
    } else {
      global_time = cctkGH->cctk_time;
    }
    
    // Unset number of time levels
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        groupdata.at(group).info.activetimelevels = deadbeef;
      }
    }
    
    reflevel = -1;
    timereflevelfact = timereffacts.at (reflevels - 1);
    // TODO: use spacereffacts.at (reflevel - 1) instead?
    spacereflevelfact = ivect(-deadbeef);
    ivect::ref(cctkGH->cctk_levfac) = spacereflevelfact;
    cctkGH->cctk_timefac = timereflevelfact;
    
    assert (is_global_mode());
  }
  
  
  
  // Set singlemap mode
  
  void enter_singlemap_mode (cGH * const cctkGH,
                             int const m, int const grouptype)
  {
    assert (is_level_mode());
    assert (mc_grouptype == -1);
    assert (m>=0 and m<maps);
    assert (grouptype == CCTK_SCALAR or grouptype == CCTK_ARRAY
            or grouptype == CCTK_GF);
    Checkpoint ("Entering singlemap mode");
    
    assert (mc_grouptype == -1);
    mc_grouptype = grouptype;
    carpetGH.map = map = m;
    
    if (mc_grouptype == CCTK_GF) {
      
      if (maps > 1) {
        // Set space delta
        for (int d=0; d<dim; ++d) {
          cctkGH->cctk_origin_space[d] = origin_space.at(map).at(mglevel)[d];
          cctkGH->cctk_delta_space[d] = delta_space.at(map)[d] * mglevelfact;
        }
      }
      
      // Set grid shape
      const ibbox& coarseext = vhh.at(map)->baseextents.at(mglevel).at(0       );
      const ibbox& baseext   = vhh.at(map)->baseextents.at(mglevel).at(reflevel);
//       assert (all (baseext.lower() % baseext.stride() == 0));
      ivect::ref(cctkGH->cctk_levoff) = baseext.lower() - coarseext.lower();
      ivect::ref(cctkGH->cctk_levoffdenom) = baseext.stride();
      ivect::ref(cctkGH->cctk_gsh) = baseext.shape() / baseext.stride();
      assert (all (vdd.at(map)->ghost_width[0] == vdd.at(map)->ghost_width[1]));
      ivect::ref(cctkGH->cctk_nghostzones) = vdd.at(map)->ghost_width[0];
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
            = ivect::ref(cctkGH->cctk_gsh);
          ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
            = ivect::ref(cctkGH->cctk_nghostzones);
        }
      }
      
    } // if mc_grouptype
    
    assert (is_singlemap_mode());
  }
  
  void leave_singlemap_mode (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_singlemap_mode() or is_level_mode());
    
    if (map == -1) return;      // early return
    
    Checkpoint ("Leaving singlemap mode");

    assert (mc_grouptype == CCTK_SCALAR or mc_grouptype == CCTK_ARRAY
            or mc_grouptype == CCTK_GF);
    
    if (mc_grouptype == CCTK_GF) {
      
      // Save space delta
      // (Do this early and often, so that interpolation has access to
      // the correct values right away.)
      for (int d=0; d<dim; ++d) {
        origin_space.at(map).at(mglevel)[d] = cctkGH->cctk_origin_space[d];
        delta_space.at(map)[d] = cctkGH->cctk_delta_space[d] / mglevelfact;
      }
      if (maps > 1) {
        // Unset space delta
        for (int d=0; d<dim; ++d) {
          cctkGH->cctk_origin_space[d] = -424242.0;
          cctkGH->cctk_delta_space[d] = -424242.0;
        }
      }
      
      // Unset grid shape
      ivect::ref(cctkGH->cctk_levoff) = deadbeef;
      ivect::ref(cctkGH->cctk_levoffdenom) = 0;
      ivect::ref(cctkGH->cctk_gsh) = deadbeef;
//       ivect::ref(cctkGH->cctk_nghostzones) = deadbeef;
      ivect::ref(cctkGH->cctk_nghostzones) = vdd.at(map)->ghost_width[0];
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
            = ivect::ref(cctkGH->cctk_gsh);
          ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
            = ivect::ref(cctkGH->cctk_nghostzones);
        }
      }
      
    } // if mc_grouptype
    
    mc_grouptype = -1;
    carpetGH.map = map = -1;
    
    assert (is_level_mode());
  }
  
  
  
  // Set local mode
  
  void enter_local_mode (cGH * const cctkGH, int const c, int const grouptype)
  {
    assert (is_singlemap_mode());
    if (mc_grouptype == CCTK_GF) {
      assert (c>=0 and c<vhh.at(map)->components(reflevel));
    } else {
      assert (c>=0 and c<CCTK_nProcs(cctkGH));
    }
    Checkpoint ("Entering local mode");
    
    assert (grouptype == mc_grouptype);
    component = c;
    
    if (mc_grouptype == CCTK_GF) {
      
      // Set cGH fields
      const ibbox& baseext = vhh.at(map)->baseextents.at(mglevel).at(reflevel);
      const ibbox& ext = vdd.at(map)->boxes.at(mglevel).at(reflevel).at(component).exterior;
      const b2vect& obnds = vhh.at(map)->outer_boundaries(reflevel,component);
      
      ivect::ref(cctkGH->cctk_lsh) = ext.shape() / ext.stride();
      ivect::ref(cctkGH->cctk_lbnd)
        = (ext.lower() - baseext.lower()) / ext.stride();
      ivect::ref(cctkGH->cctk_ubnd)
        = (ext.upper() - baseext.lower()) / ext.stride();
      ivect::ref(cctkGH->cctk_from) = 0;
      ivect::ref(cctkGH->cctk_to) = ivect::ref(cctkGH->cctk_lsh);
      
      for (int d=0; d<dim; ++d) {
        cctkGH->cctk_bbox[2*d  ] = obnds[0][d];
        cctkGH->cctk_bbox[2*d+1] = obnds[1][d];
      }
      
      for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
        for (int d=0; d<dim; ++d) {
          // TODO: support staggering
          cctkGH->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cctkGH->cctk_lsh[d];
        }
      }
      
      for (int d=0; d<dim; ++d) {
        assert (cctkGH->cctk_lsh[d] >= 0);
        assert (cctkGH->cctk_lsh[d] <= cctkGH->cctk_gsh[d]);
        assert (cctkGH->cctk_lbnd[d] >= 0);
        assert (cctkGH->cctk_lbnd[d] <= cctkGH->cctk_ubnd[d] + 1);
        assert (cctkGH->cctk_ubnd[d] < cctkGH->cctk_gsh[d]);
        assert (cctkGH->cctk_lbnd[d] + cctkGH->cctk_lsh[d] - 1 == cctkGH->cctk_ubnd[d]);
        assert (cctkGH->cctk_lbnd[d] <= cctkGH->cctk_ubnd[d]+1);
        assert (cctkGH->cctk_from[d] >= 0);
        assert (cctkGH->cctk_from[d] <= cctkGH->cctk_to[d]);
        assert (cctkGH->cctk_to[d] <= cctkGH->cctk_lsh[d]);
      }
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          
          ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
            = ivect::ref(cctkGH->cctk_lsh);
          ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
            = ivect::ref(cctkGH->cctk_lbnd);
          ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
            = ivect::ref(cctkGH->cctk_ubnd);
          
          for (int d=0; d<dim; ++d) {
            const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ]
              = cctkGH->cctk_bbox[2*d  ];
            const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1]
              = cctkGH->cctk_bbox[2*d+1];
          }
          
          const int numvars = CCTK_NumVarsInGroupI (group);
          if (numvars>0) {
            const int firstvar = CCTK_FirstVarIndexI (group);
            assert (firstvar>=0);
            const int max_tl = CCTK_MaxTimeLevelsGI (group);
            assert (max_tl>=0);
            const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, group);
            assert (active_tl>=0 and active_tl<=max_tl);
            const int available_tl =
              do_allow_past_timelevels ? active_tl : min (1, active_tl);
            
            // assert (vhh.at(map)->is_local(reflevel,component));
            
            assert (group<(int)arrdata.size());
            for (int var=0; var<numvars; ++var) {
              assert (firstvar+var<CCTK_NumVars());
              ggf * const ff = arrdata.at(group).at(map).data.at(var);
              for (int tl=0; tl<max_tl; ++tl) {
                if (ff and tl<available_tl) {
                  gdata * const data = (*ff) (tl, reflevel, component, mglevel);
                  assert (data);
                  cctkGH->data[firstvar+var][tl] = data->storage();
                } else {
                  cctkGH->data[firstvar+var][tl] = NULL;
                }
              }
            }
          }
          
        } // if grouptype
      } // for group
      
    } // if mc_grouptype
    
    assert (is_local_mode());
  }
  
  void leave_local_mode (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_local_mode() or is_singlemap_mode());
    
    if (component == -1) return; // early return

    Checkpoint ("Leaving local mode");
    
    if (mc_grouptype == CCTK_GF) {
      
      // Unset cGH fields
      ivect::ref(cctkGH->cctk_lsh) = deadbeef;
      ivect::ref(cctkGH->cctk_lbnd) = -deadbeef;
      ivect::ref(cctkGH->cctk_ubnd) = deadbeef;
      ivect::ref(cctkGH->cctk_from) = -deadbeef;
      ivect::ref(cctkGH->cctk_to) = deadbeef;
      
      for (int d=0; d<dim; ++d) {
        cctkGH->cctk_bbox[2*d  ] = deadbeef;
        cctkGH->cctk_bbox[2*d+1] = deadbeef;
      }
      
      for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
        for (int d=0; d<dim; ++d) {
          // TODO: support staggering
          cctkGH->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cctkGH->cctk_lsh[d];
        }
      }
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          
          ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
            = ivect::ref(cctkGH->cctk_lsh);
          ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
            = ivect::ref(cctkGH->cctk_lbnd);
          ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
            = ivect::ref(cctkGH->cctk_ubnd);
          
          for (int d=0; d<dim; ++d) {
            const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ]
              = cctkGH->cctk_bbox[2*d  ];
            const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1]
              = cctkGH->cctk_bbox[2*d+1];
          }
          
          const int numvars = CCTK_NumVarsInGroupI (group);
          if (numvars>0) {
            const int firstvar = CCTK_FirstVarIndexI (group);
            assert (firstvar>=0);
            const int max_tl = CCTK_MaxTimeLevelsGI (group);
            assert (max_tl>=0);
            
            assert (group<(int)arrdata.size());
            for (int var=0; var<numvars; ++var) {
              assert (firstvar+var<CCTK_NumVars());
              for (int tl=0; tl<max_tl; ++tl) {
                cctkGH->data[firstvar+var][tl] = NULL;
              }
            }
          }
          
        } // if grouptype
      } // for group
      
    } // if mc_grouptype
    
    component = -1;
    
    assert (is_singlemap_mode());
  }
  
  
  
  //
  // Mode iterators
  //
  
  // mglevel iterator
  
  mglevel_iterator::mglevel_iterator (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), ml(mglevels-1)
  {
    enter_global_mode (cctkGH, ml);
  }
  
  mglevel_iterator::~mglevel_iterator ()
  {
    leave_global_mode (cctkGH);
  }
  
  bool mglevel_iterator::done () const
  {
    return ml < 0;
  }
  
  void mglevel_iterator::step ()
  {
    -- ml;
    if (not done()) {
      leave_global_mode (cctkGH);
      enter_global_mode (cctkGH, ml);
    }
  }
  
  
  
  // reverse mglevel iterator
  
  reverse_mglevel_iterator::reverse_mglevel_iterator (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), ml(0)
  {
    enter_global_mode (cctkGH, ml);
  }
  
  reverse_mglevel_iterator::~reverse_mglevel_iterator ()
  {
    leave_global_mode (cctkGH);
  }
  
  bool reverse_mglevel_iterator::done () const
  {
    return ml >= mglevels;
  }
  
  void reverse_mglevel_iterator::step ()
  {
    ++ ml;
    if (not done()) {
      leave_global_mode (cctkGH);
      enter_global_mode (cctkGH, ml);
    }
  }
  
  
  
  // reflevel iterator
  
  reflevel_iterator::reflevel_iterator (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), rl(0)
  {
    enter_level_mode (cctkGH, rl);
  }
  
  reflevel_iterator::~reflevel_iterator ()
  {
    leave_level_mode (cctkGH);
  }
  
  bool reflevel_iterator::done () const
  {
    return rl >= reflevels;
  }
  
  void reflevel_iterator::step ()
  {
    ++ rl;
    if (not done()) {
      leave_level_mode (cctkGH);
      enter_level_mode (cctkGH, rl);
    }
  }
  
  
  
  // reverse reflevel iterator
  
  reverse_reflevel_iterator::reverse_reflevel_iterator (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), rl(reflevels-1)
  {
    enter_level_mode (cctkGH, rl);
  }
  
  reverse_reflevel_iterator::~reverse_reflevel_iterator ()
  {
    leave_level_mode (cctkGH);
  }
  
  bool reverse_reflevel_iterator::done () const
  {
    return rl < 0;
  }
  
  void reverse_reflevel_iterator::step ()
  {
    -- rl;
    if (not done()) {
      leave_level_mode (cctkGH);
      enter_level_mode (cctkGH, rl);
    }
  }
  
  
  
  // map iterator
  
  map_iterator::map_iterator (cGH const * const cctkGH_,
                              int const grouptype_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), grouptype(grouptype_), m(0)
  {
    assert (grouptype == CCTK_GF
            or grouptype == CCTK_ARRAY or grouptype == CCTK_SCALAR);
    enter_singlemap_mode (cctkGH, m, grouptype);
  }
  
  map_iterator::~map_iterator ()
  {
    leave_singlemap_mode (cctkGH);
  }
  
  bool map_iterator::done () const
  {
    return grouptype == CCTK_GF ? m >= maps : m >= 1;
  }
  
  void map_iterator::step ()
  {
    ++ m;
    if (not done()) {
      leave_singlemap_mode (cctkGH);
      enter_singlemap_mode (cctkGH, m, grouptype);
    }
  }
  
  
  
  // Local mode iterator
  
  component_iterator::component_iterator (cGH const * const cctkGH_,
                                          int const grouptype_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), grouptype(grouptype_), c(-1)
  {
    assert (grouptype == CCTK_GF
            or grouptype == CCTK_ARRAY or grouptype == CCTK_SCALAR);
    assert (is_singlemap_mode());
    step ();
  }
  
  component_iterator::~component_iterator ()
  {
    leave_local_mode (cctkGH);
  }
  
  bool component_iterator::done () const
  {
    int const maxc = (grouptype == CCTK_GF
                      ? vhh.AT(map)->components(reflevel)
                      : dist::size());
    return c >= maxc;
  }
  
  void component_iterator::step ()
  {
    ++ c;
    if (not done()) {
      leave_local_mode (cctkGH);
      enter_local_mode (cctkGH, c, grouptype);
    }
  }
  
  
  
  // Processor local mode iterator
  
  local_component_iterator::local_component_iterator (cGH const * const cctkGH_,
                                                      int const grouptype_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), grouptype(grouptype_), c(-1)
  {
    assert (grouptype == CCTK_GF
            or grouptype == CCTK_ARRAY or grouptype == CCTK_SCALAR);
    assert (is_singlemap_mode());
    step ();
  }
  
  local_component_iterator::~local_component_iterator ()
  {
    leave_local_mode (cctkGH);
  }
  
  bool local_component_iterator::done () const
  {
    int const maxc = (grouptype == CCTK_GF
                      ? vhh.AT(map)->components(reflevel)
                      : dist::size());
    return c >= maxc;
  }
  
  void local_component_iterator::step ()
  {
    do {
      ++ c;
    } while (not done() and not (grouptype == CCTK_GF
                                 ? vhh.at(map)->is_local(reflevel, c)
                                 : c == CCTK_MyProc(cctkGH)));
    
    if (not done()) {
      leave_local_mode (cctkGH);
      enter_local_mode (cctkGH, c, grouptype);
    }
  }
  
  
  
  //
  // Mode escapes
  //
  
  // Singlemap escape
  
  singlemap_escape::singlemap_escape (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)), c(component)
  {
    assert (not is_meta_mode());
    assert (not is_global_mode());
    assert (not is_level_mode());
    if (not is_singlemap_mode()) {
      leave_local_mode (cctkGH);
    }
  }
  
  singlemap_escape::~singlemap_escape ()
  {
    assert (is_singlemap_mode());
    if (c != -1) {
      enter_local_mode (cctkGH, c, mc_grouptype);
    }
  }
  
  
  
  // Level escape
  
  level_escape::level_escape (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)),
      grouptype(mc_grouptype), m(map), c(component)
  {
    assert (not is_meta_mode());
    assert (not is_global_mode());
    if (not is_level_mode()) {
      if (not is_singlemap_mode()) {
        leave_local_mode (cctkGH);
      }
      leave_singlemap_mode (cctkGH);
    }
  }
  
  level_escape::~level_escape ()
  {
    assert (is_level_mode());
    if (m != -1) {
      enter_singlemap_mode (cctkGH, m, grouptype);
      if (c != -1) {
        enter_local_mode (cctkGH, c, grouptype);
      }
    }
  }
  
  
  
  // Global escape
  
  global_escape::global_escape (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)),
      rl(reflevel), grouptype(mc_grouptype), m(map), c(component)
  {
    assert (not is_meta_mode());
    if (not is_global_mode()) {
      if (not is_level_mode()) {
        if (not is_singlemap_mode()) {
          leave_local_mode (cctkGH);
        }
        leave_singlemap_mode (cctkGH);
      }
      leave_level_mode (cctkGH);
    }
  }
  
  global_escape::~global_escape ()
  {
    assert (is_global_mode());
    if (rl != -1) {
      enter_level_mode (cctkGH, rl);
      if (m != -1) {
        enter_singlemap_mode (cctkGH, m, grouptype);
        if (c != -1) {
          enter_local_mode (cctkGH, c, grouptype);
        }
      }
    }
  }
  
  
  
  // Meta escape
  
  meta_escape::meta_escape (cGH const * const cctkGH_)
    : cctkGH(const_cast<cGH*>(cctkGH_)),
      ml(mglevel), rl(reflevel), grouptype(mc_grouptype), m(map), c(component)
  {
    if (not is_meta_mode()) {
      if (not is_global_mode()) {
        if (not is_level_mode()) {
          if (not is_singlemap_mode()) {
            leave_local_mode (cctkGH);
          }
          leave_singlemap_mode (cctkGH);
        }
        leave_level_mode (cctkGH);
      }
      leave_global_mode (cctkGH);
    }
  }
  
  meta_escape::~meta_escape ()
  {
    assert (is_meta_mode());
    if (ml != -1) {
      enter_global_mode (cctkGH, ml);
      if (rl != -1) {
        enter_level_mode (cctkGH, rl);
        if (m != -1) {
          enter_singlemap_mode (cctkGH, m, grouptype);
          if (c != -1) {
            enter_local_mode (cctkGH, c, grouptype);
          }
        }
      }
    }
  }
  
  
  
  //
  // Call functions in specific modes
  //
  
  int CallLocalFunction (cGH * const cctkGH,
                         void (* const function) (cGH * const cctkGH))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        BEGIN_REFLEVEL_LOOP(cctkGH) {
          BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
            BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              function (cctkGH);
            } END_LOCAL_COMPONENT_LOOP;
          } END_MAP_LOOP;
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    } else if (is_global_mode()) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            function (cctkGH);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } END_REFLEVEL_LOOP;
    } else if (is_level_mode()) {
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          function (cctkGH);
        } END_LOCAL_COMPONENT_LOOP;
      } END_MAP_LOOP;
    } else if (is_singlemap_mode()) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        function (cctkGH);
      } END_LOCAL_COMPONENT_LOOP;
    } else {
      function (cctkGH);
    }
    return 0;
  }
  
  int CallSinglemapFunction (cGH * const cctkGH,
                             void (* const function) (cGH * const cctkGH))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        BEGIN_REFLEVEL_LOOP(cctkGH) {
          BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
            function (cctkGH);
          } END_MAP_LOOP;
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    } else if (is_global_mode()) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          function (cctkGH);
        } END_MAP_LOOP;
      } END_REFLEVEL_LOOP;
    } else if (is_level_mode()) {
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        function (cctkGH);
      } END_MAP_LOOP;
    } else {
      BEGIN_SINGLEMAP_MODE(cctkGH) {
        function (cctkGH);
      } END_SINGLEMAP_MODE;
    }
    return 0;
  }
  
  int CallLevelFunction (cGH * const cctkGH,
                         void (* const function) (cGH * const cctkGH))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        BEGIN_REFLEVEL_LOOP(cctkGH) {
          function (cctkGH);
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    } else if (is_global_mode()) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        function (cctkGH);
      } END_REFLEVEL_LOOP;
    } else {
      BEGIN_LEVEL_MODE(cctkGH) {
        function (cctkGH);
      } END_LEVEL_MODE;
    }
    return 0;
  }
  
  int CallGlobalFunction (cGH * const cctkGH,
                          void (* const function) (cGH * const cctkGH))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        function (cctkGH);
      } END_MGLEVEL_LOOP;
    } else {
      BEGIN_GLOBAL_MODE(cctkGH) {
        function (cctkGH);
      } END_GLOBAL_MODE;
    }
    return 0;
  }
  
  int CallMetaFunction (cGH * const cctkGH,
                        void (* const function) (cGH * const cctkGH))
  {
    BEGIN_META_MODE(cctkGH) {
      function (cctkGH);
    } END_META_MODE;
    return 0;
  }
  
  
  
  //
  // Mode setters
  //
  
  // mglevel setter
  
  mglevel_setter::mglevel_setter (cGH const * const cctkGH_, int const ml)
    : cctkGH(const_cast<cGH*>(cctkGH_))
  {
    assert (is_meta_mode());
    enter_global_mode (cctkGH, ml);
  }
  
  mglevel_setter::~mglevel_setter ()
  {
    leave_global_mode (cctkGH);
  }
  
  // reflevel setter
  
  reflevel_setter::reflevel_setter (cGH const * const cctkGH_, int const rl)
    : cctkGH(const_cast<cGH*>(cctkGH_))
  {
    assert (is_global_mode());
    enter_level_mode (cctkGH, rl);
  }
  
  reflevel_setter::~reflevel_setter ()
  {
    leave_level_mode (cctkGH);
  }
  
  // map setter
  
  map_setter::map_setter (cGH const * const cctkGH_,
                          int const m, int const grouptype)
    : cctkGH(const_cast<cGH*>(cctkGH_))
  {
    assert (is_level_mode());
    enter_singlemap_mode (cctkGH, m, grouptype);
  }
  
  map_setter::~map_setter ()
  {
    leave_singlemap_mode (cctkGH);
  }
  
  // component setter
  
  component_setter::component_setter (cGH const * const cctkGH_,
                                      int const c, int const grouptype)
    : cctkGH(const_cast<cGH*>(cctkGH_))
  {
    assert (is_singlemap_mode());
    enter_local_mode (cctkGH, c, grouptype);
  }
  
  component_setter::~component_setter ()
  {
    leave_local_mode (cctkGH);
  }
  
  
  
  //
  // Call a scheduling group
  //
  
  int CallScheduleGroup (cGH * const cctkGH, const char * const group)
  {
    CCTK_ScheduleTraverse (group, cctkGH, CallFunction);
    return 0;
  }
  
  extern "C" void CCTK_FCALL CCTK_FNAME(CallScheduleGroup)
    (int * const ierr, cGH * const * const cctkGH, ONE_FORTSTRING_ARG)
  {
    ONE_FORTSTRING_CREATE (group);
    *ierr = CallScheduleGroup (*cctkGH, group);
    free (group);
  }
  
} // namespace Carpet
