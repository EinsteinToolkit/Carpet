#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "defs.hh"
#include "ggf.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/modes.cc,v 1.9 2004/08/07 20:07:27 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_modes_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  //
  // Mode indicators
  //
  
  bool is_meta_mode ()
  {
    if (mglevel==-1) assert (reflevel==-1 && map==-1 && component==-1);
    return mglevel==-1;
  }
  
  bool is_global_mode ()
  {
    if (mglevel==-1 && reflevel!=-1) assert (map==-1 && component==-1);
    return mglevel!=-1 && reflevel==-1;
  }
  
  bool is_level_mode ()
  {
    if (mglevel!=-1 && reflevel!=-1 && map==-1) assert (component==-1);
    return mglevel!=-1 && reflevel!=-1 && map==-1;
  }
  
  bool is_singlemap_mode ()
  {
    return mglevel!=-1 && reflevel!=-1 && map!=-1 && component==-1;
  }
  
  bool is_local_mode ()
  {
    return mglevel!=-1 && reflevel!=-1 && map!=-1 && component!=-1;
//       assert (mglevel>=0 && mglevel<mglevels);
//       assert (reflevel>=0 && reflevel<reflevels);
//       assert (map>=0 && map<maps);
//       assert (vhh.at(map)->local_components(reflevel)==1 || component==-1);
  }
  
  
  
  //
  // Mode setting
  //
  
  // Set global mode
  
  void enter_global_mode (cGH * const cgh, int const ml)
  {
    assert (is_meta_mode());
    assert (ml>=0 && ml<mglevels);
    Checkpoint ("Entering global mode");
    
    mglevel = ml;
    mglevelfact = ipow(mgfact, mglevel);
    // TODO: this could also just be "mglevel" instead
    cgh->cctk_convlevel = basemglevel + mglevel;
    
    // Set time and space delta
    cgh->cctk_delta_time = delta_time * mglevelfact;
    for (int d=0; d<dim; ++d) {
      cgh->cctk_origin_space[d] = origin_space.at(mglevel)[d];
      cgh->cctk_delta_space[d] = delta_space[d] * mglevelfact;
    }
    
    // Set array information
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) != CCTK_GF) {
        
        const int rl = 0;
        const int m = 0;
        const int c = CCTK_MyProc(cgh);
        
        const ibbox& base = arrdata.at(group).at(m).hh->bases.at(rl).at(ml);
        const bbvect& obnds = arrdata.at(group).at(m).hh->outer_boundaries.at(rl).at(c);
	const ibbox& ext = arrdata.at(group).at(m).dd->boxes.at(rl).at(c).at(ml).exterior;
        
        ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
          = arrdata.at(group).at(m).dd->lghosts;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
          = base.shape() / base.stride();
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
          = ext.shape() / ext.stride();
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
          = (ext.lower() - base.lower()) / ext.stride();
        ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
          = (ext.upper() - base.lower()) / ext.stride();
        for (int d=0; d<dim; ++d) {
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ] = obnds[d][0];
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1] = obnds[d][1];
        }
        
        for (int d=0; d<dim; ++d) {
          assert (groupdata.at(group).info.lsh[d]>=0);
          assert (groupdata.at(group).info.lsh[d]<=groupdata.at(group).info.gsh[d]);
          if (d>=groupdata.at(group).info.dim) {
            assert (groupdata.at(group).info.lsh[d]==1);
          }
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
          const int num_tl = CCTK_NumTimeLevelsFromVarI (firstvar);
          
          assert (arrdata.at(group).at(m).hh->is_local(rl,c));
          
          for (int var=0; var<numvars; ++var) {
            assert (firstvar+var<CCTK_NumVars());
            for (int tl=0; tl<num_tl; ++tl) {
              ggf<dim> * const ff = arrdata.at(group).at(m).data.at(var);
              if (ff) {
                gdata<dim> * const data = (*ff) (-tl, rl, c, ml);
                assert (data);
                cgh->data[firstvar+var][tl] = data->storage();
              } else {
                cgh->data[firstvar+var][tl] = 0;
              }
            }
          }
        }
        
      } // if grouptype
    } // for group
    
    assert (is_global_mode());
  }
  
  void leave_global_mode (cGH * const cgh)
  {
    assert (is_global_mode() || is_meta_mode());
    Checkpoint ("Leaving global mode");
    
    if (mglevel == -1) return;  // early return
    
    // Save and unset time and space delta
    delta_time = cgh->cctk_delta_time / mglevelfact;
    cgh->cctk_delta_time = 0.0;
    for (int d=0; d<dim; ++d) {
      origin_space.at(mglevel)[d] = cgh->cctk_origin_space[d];
      delta_space[d] = cgh->cctk_delta_space[d] / mglevelfact;
      cgh->cctk_origin_space[d] = -424242.0;
      cgh->cctk_delta_space[d] = 0.0;
    }
   
    // Set array information
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) != CCTK_GF) {
        
        const int m = 0;
        
//         ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
//           = 0xdead;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
          = arrdata.at(group).at(m).dd->lghosts;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
          = 0xdead;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
          = 0xdead;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
          = -0xdead;
        ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
          = 0xdead;
        for (int d=0; d<dim; ++d) {
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ] = 0xdead;
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1] = 0xdead;
        }
        
        const int numvars = CCTK_NumVarsInGroupI (group);
        if (numvars>0) {
          const int firstvar = CCTK_FirstVarIndexI (group);
          assert (firstvar>=0);
          const int num_tl = CCTK_NumTimeLevelsFromVarI (firstvar);
          
          assert (group<(int)arrdata.size());
          for (int var=0; var<numvars; ++var) {
            assert (firstvar+var<CCTK_NumVars());
            for (int tl=0; tl<num_tl; ++tl) {
              cgh->data[firstvar+var][tl] = 0;
            }
          }
        }
        
      } // if grouptype
    } // for group
    
    mglevel = -1;
    mglevelfact = -0xdead;
    cgh->cctk_convlevel = -0xdead;
    
    assert (is_meta_mode());
  }
  
  
  
  // Set level mode
  
  void enter_level_mode (cGH * const cgh, int const rl)
  {
    assert (is_global_mode());
    assert (rl>=0 && rl<reflevels);
    Checkpoint ("Entering level mode");
    
    reflevel = rl;
    reflevelfact = ipow(reffact, reflevel);
    ivect::ref(cgh->cctk_levfac) = reflevelfact;
    cgh->cctk_timefac = reflevelfact;
    
    // Set current time
    assert (mglevel>=0 && mglevel<(int)leveltimes.size());
    assert (reflevel>=0 && reflevel<(int)leveltimes.at(mglevel).size());
    cgh->cctk_time = leveltimes.at(mglevel).at(reflevel);
    
    assert (is_level_mode());
  }
  
  void leave_level_mode (cGH * const cgh)
  {
    assert (is_level_mode() || is_global_mode());
    Checkpoint ("Leaving level mode");
    
    if (reflevel == -1) return; // early return
    
    // Save and unset current time
    assert (mglevel>=0 && mglevel<(int)leveltimes.size());
    assert (reflevel>=0 && reflevel<(int)leveltimes.at(mglevel).size());
    leveltimes.at(mglevel).at(reflevel) = cgh->cctk_time;
    cgh->cctk_time = global_time;
    
    reflevel = -1;
    reflevelfact = -0xdead;
    ivect::ref(cgh->cctk_levfac) = -0xdead;
    cgh->cctk_timefac = -0xdead;
    
    assert (is_global_mode());
  }
  
  
  
  // Set singlemap mode
  
  void enter_singlemap_mode (cGH * const cgh, int const m)
  {
    assert (is_level_mode());
    assert (m>=0 && m<maps);
    Checkpoint ("Entering singlemap mode");
    
    carpetGH.map = map = m;
    
    // Set grid shape
    const ibbox& coarseext = vdd.at(map)->bases.at(0       ).at(mglevel).exterior;
    const ibbox& baseext   = vdd.at(map)->bases.at(reflevel).at(mglevel).exterior;
    assert (all (baseext.lower() % baseext.stride() == 0));
    assert (all ((baseext.lower() - coarseext.lower()) % baseext.stride() == 0));
    ivect::ref(cgh->cctk_levoff) = (baseext.lower() - coarseext.lower()) / baseext.stride();
    ivect::ref(cgh->cctk_levoffdenom) = 1;
    ivect::ref(cgh->cctk_gsh) = baseext.shape() / baseext.stride();
    assert (all (vdd.at(map)->lghosts == vdd.at(map)->ughosts));
    ivect::ref(cgh->cctk_nghostzones) = vdd.at(map)->lghosts;
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
          = ivect::ref(cgh->cctk_gsh);
        ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
          = ivect::ref(cgh->cctk_nghostzones);
      }
    }
    
    assert (is_singlemap_mode());
  }
  
  void leave_singlemap_mode (cGH * const cgh)
  {
    assert (is_singlemap_mode() || is_level_mode());
    Checkpoint ("Leaving singlemap mode");
    
    if (map == -1) return;      // early return
    
    // Unset grid shape
    ivect::ref(cgh->cctk_levoff) = 0xdead;
    ivect::ref(cgh->cctk_levoffdenom) = 0;
    ivect::ref(cgh->cctk_gsh) = 0xdead;
//     ivect::ref(cgh->cctk_nghostzones) = 0xdead;
    ivect::ref(cgh->cctk_nghostzones) = vdd.at(map)->lghosts;
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        ivect::ref(const_cast<int*>(groupdata.at(group).info.gsh))
          = ivect::ref(cgh->cctk_gsh);
        ivect::ref(const_cast<int*>(groupdata.at(group).info.nghostzones))
          = ivect::ref(cgh->cctk_nghostzones);
      }
    }
    
    carpetGH.map = map = -1;
    
    assert (is_level_mode());
  }
  
  
  
  // Set local mode
  
  void enter_local_mode (cGH * const cgh, int const c)
  {
    assert (is_singlemap_mode());
    assert (c>=0 && c<vhh.at(map)->components(reflevel));
    Checkpoint ("Entering local mode");
    
    component = c;
    
    // Set cGH fields
    const ibbox& baseext = vdd.at(map)->bases.at(reflevel).at(mglevel).exterior;
    const bbvect& obnds = vhh.at(map)->outer_boundaries.at(reflevel).at(component);
    const ibbox& ext = vdd.at(map)->boxes.at(reflevel).at(component).at(mglevel).exterior;
    
    ivect::ref(cgh->cctk_lsh) = ext.shape() / ext.stride();
    ivect::ref(cgh->cctk_lbnd)
      = (ext.lower() - baseext.lower()) / ext.stride();
    ivect::ref(cgh->cctk_ubnd)
      = (ext.upper() - baseext.lower()) / ext.stride();
    
    for (int d=0; d<dim; ++d) {
      cgh->cctk_bbox[2*d  ] = obnds[d][0];
      cgh->cctk_bbox[2*d+1] = obnds[d][1];
    }
      
    for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
      for (int d=0; d<dim; ++d) {
        // TODO: support staggering
        cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
      }
    }
    
    for (int d=0; d<dim; ++d) {
      assert (cgh->cctk_lsh[d] >= 0);
      assert (cgh->cctk_lsh[d] <= cgh->cctk_gsh[d]);
      assert (cgh->cctk_lbnd[d] >= 0);
      assert (cgh->cctk_lbnd[d] <= cgh->cctk_ubnd[d] + 1);
      assert (cgh->cctk_ubnd[d] < cgh->cctk_gsh[d]);
      assert (cgh->cctk_lbnd[d] + cgh->cctk_lsh[d] - 1 == cgh->cctk_ubnd[d]);
      assert (cgh->cctk_lbnd[d] <= cgh->cctk_ubnd[d]+1);
    }
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
          = ivect::ref(cgh->cctk_lsh);
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
          = ivect::ref(cgh->cctk_lbnd);
        ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
          = ivect::ref(cgh->cctk_ubnd);
        
        for (int d=0; d<dim; ++d) {
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ]
            = cgh->cctk_bbox[2*d  ];
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1]
            = cgh->cctk_bbox[2*d+1];
        }
        
        const int numvars = CCTK_NumVarsInGroupI (group);
        if (numvars>0) {
          const int firstvar = CCTK_FirstVarIndexI (group);
          assert (firstvar>=0);
          const int num_tl = CCTK_NumTimeLevelsFromVarI (firstvar);
          
//           assert (vhh.at(map)->is_local(reflevel,component));
          
          assert (group<(int)arrdata.size());
          for (int var=0; var<numvars; ++var) {
            assert (firstvar+var<CCTK_NumVars());
            for (int tl=0; tl<num_tl; ++tl) {
              ggf<dim> * const ff = arrdata.at(group).at(map).data.at(var);
              if (ff) {
                gdata<dim> * const data
                  = (*ff) (-tl, reflevel, component, mglevel);
                assert (data);
                cgh->data[firstvar+var][tl] = data->storage();
              } else {
                cgh->data[firstvar+var][tl] = 0;
              }
            }
          }
        }
        
      } // if grouptype
    } // for group
    
    assert (is_local_mode());
  }
  
  void leave_local_mode (cGH * const cgh)
  {
    assert (is_local_mode() || is_singlemap_mode());
    Checkpoint ("Leaving local mode");
    
    if (component == -1) return; // early return
    
    // Unset cGH fields
    ivect::ref(cgh->cctk_lsh) = 0xdead;
    ivect::ref(cgh->cctk_lbnd) = -0xdead;
    ivect::ref(cgh->cctk_ubnd) = 0xdead;
    
    for (int d=0; d<dim; ++d) {
      cgh->cctk_bbox[2*d  ] = 0xdead;
      cgh->cctk_bbox[2*d+1] = 0xdead;
    }
      
    for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
      for (int d=0; d<dim; ++d) {
        // TODO: support staggering
        cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
      }
    }
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lsh))
          = ivect::ref(cgh->cctk_lsh);
        ivect::ref(const_cast<int*>(groupdata.at(group).info.lbnd))
          = ivect::ref(cgh->cctk_lbnd);
        ivect::ref(const_cast<int*>(groupdata.at(group).info.ubnd))
          = ivect::ref(cgh->cctk_ubnd);
        
        for (int d=0; d<dim; ++d) {
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d  ]
            = cgh->cctk_bbox[2*d  ];
          const_cast<int*>(groupdata.at(group).info.bbox)[2*d+1]
            = cgh->cctk_bbox[2*d+1];
        }
        
        const int numvars = CCTK_NumVarsInGroupI (group);
        if (numvars>0) {
          const int firstvar = CCTK_FirstVarIndexI (group);
          assert (firstvar>=0);
          const int num_tl = CCTK_NumTimeLevelsFromVarI (firstvar);
          
          assert (group<(int)arrdata.size());
          for (int var=0; var<numvars; ++var) {
            assert (firstvar+var<CCTK_NumVars());
            for (int tl=0; tl<num_tl; ++tl) {
              cgh->data[firstvar+var][tl] = 0;
            }
          }
        }
        
      } // if grouptype
    } // for group
    
    component = -1;
    
    assert (is_singlemap_mode());
  }
  
  
  
  //
  // Mode iterators
  //
  
  // mglevel iterator
  
  mglevel_iterator::mglevel_iterator (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), ml(mglevels-1)
  {
    enter_global_mode (cgh, ml);
  }
  
  mglevel_iterator::~mglevel_iterator ()
  {
    leave_global_mode (cgh);
  }
  
  bool mglevel_iterator::done () const
  {
    return ml < 0;
  }
  
  void mglevel_iterator::step ()
  {
    -- ml;
    if (! done()) {
      leave_global_mode (cgh);
      enter_global_mode (cgh, ml);
    }
  }
  
  
  
  // reverse mglevel iterator
  
  reverse_mglevel_iterator::reverse_mglevel_iterator (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), ml(0)
  {
    enter_global_mode (cgh, ml);
  }
  
  reverse_mglevel_iterator::~reverse_mglevel_iterator ()
  {
    leave_global_mode (cgh);
  }
  
  bool reverse_mglevel_iterator::done () const
  {
    return ml >= mglevels;
  }
  
  void reverse_mglevel_iterator::step ()
  {
    ++ ml;
    if (! done()) {
      leave_global_mode (cgh);
      enter_global_mode (cgh, ml);
    }
  }
  
  
  
  // reflevel iterator
  
  reflevel_iterator::reflevel_iterator (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), rl(0)
  {
    enter_level_mode (cgh, rl);
  }
  
  reflevel_iterator::~reflevel_iterator ()
  {
    leave_level_mode (cgh);
  }
  
  bool reflevel_iterator::done () const
  {
    return rl >= reflevels;
  }
  
  void reflevel_iterator::step ()
  {
    ++ rl;
    if (! done()) {
      leave_level_mode (cgh);
      enter_level_mode (cgh, rl);
    }
  }
  
  
  
  // reverse reflevel iterator
  
  reverse_reflevel_iterator::reverse_reflevel_iterator (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), rl(reflevels-1)
  {
    enter_level_mode (cgh, rl);
  }
  
  reverse_reflevel_iterator::~reverse_reflevel_iterator ()
  {
    leave_level_mode (cgh);
  }
  
  bool reverse_reflevel_iterator::done () const
  {
    return rl < 0;
  }
  
  void reverse_reflevel_iterator::step ()
  {
    -- rl;
    if (! done()) {
      leave_level_mode (cgh);
      enter_level_mode (cgh, rl);
    }
  }
  
  
  
  // map iterator
  
  map_iterator::map_iterator (cGH const * const cgh_,
                              int const grouptype_)
    : cgh(const_cast<cGH*>(cgh_)), grouptype(grouptype_), m(0)
  {
    assert (grouptype == CCTK_GF
            || grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
    enter_singlemap_mode (cgh, m);
  }
  
  map_iterator::~map_iterator ()
  {
    leave_singlemap_mode (cgh);
  }
  
  bool map_iterator::done () const
  {
    return grouptype == CCTK_GF ? m >= maps : m >= 1;
  }
  
  void map_iterator::step ()
  {
    ++ m;
    if (! done()) {
      leave_singlemap_mode (cgh);
      enter_singlemap_mode (cgh, m);
    }
  }
  
  
  
  // Local mode iterator
  
  component_iterator::component_iterator (cGH const * const cgh_,
                                          int const grouptype_)
    : cgh(const_cast<cGH*>(cgh_)), grouptype(grouptype_), c(0)
  {
    assert (grouptype == CCTK_GF
            || grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
    enter_local_mode (cgh, c);
  }
  
  component_iterator::~component_iterator ()
  {
    leave_local_mode (cgh);
  }
  
  bool component_iterator::done () const
  {
    return (grouptype == CCTK_GF
            ? c >= vhh.at(map)->components(reflevel)
            : c >= CCTK_nProcs(cgh));
  }
  
  void component_iterator::step ()
  {
    ++ c;
    if (! done()) {
      leave_local_mode (cgh);
      enter_local_mode (cgh, c);
    }
  }
  
  
  
  // Processor local mode iterator
  
  local_component_iterator::local_component_iterator (cGH const * const cgh_,
                                                      int const grouptype_)
    : cgh(const_cast<cGH*>(cgh_)), grouptype(grouptype_), c(-1)
  {
    assert (grouptype == CCTK_GF
            || grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
    assert (is_singlemap_mode());
    step ();
  }
  
  local_component_iterator::~local_component_iterator ()
  {
    leave_local_mode (cgh);
  }
  
  bool local_component_iterator::done () const
  {
    return c >= (grouptype == CCTK_GF
                 ? vhh.at(map)->components(reflevel)
                 : CCTK_nProcs(cgh));
  }
  
  void local_component_iterator::step ()
  {
    do {
      ++ c;
    } while (! done() && ! (grouptype == CCTK_GF
                            ? vhh.at(map)->is_local(reflevel, c)
                            : c == CCTK_MyProc(cgh)));
    if (! done()) {
      leave_local_mode (cgh);
      enter_local_mode (cgh, c);
    }
  }
  
  
  
  //
  // Mode escapes
  //
  
  // Singlemap escape
  
  singlemap_escape::singlemap_escape (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), c(component)
  {
    assert (! is_meta_mode());
    assert (! is_global_mode());
    assert (! is_level_mode());
    if (! is_singlemap_mode()) {
      leave_local_mode (cgh);
    }
  }
  
  singlemap_escape::~singlemap_escape ()
  {
    assert (is_singlemap_mode());
    if (c != -1) {
      enter_local_mode (cgh, c);
    }
  }
  
  
  
  // Level escape
  
  level_escape::level_escape (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), m(map), c(component)
  {
    assert (! is_meta_mode());
    assert (! is_global_mode());
    if (! is_level_mode()) {
      if (! is_singlemap_mode()) {
        leave_local_mode (cgh);
      }
      leave_singlemap_mode (cgh);
    }
  }
  
  level_escape::~level_escape ()
  {
    assert (is_level_mode());
    if (m != -1) {
      enter_singlemap_mode (cgh, m);
      if (c != -1) {
        enter_local_mode (cgh, c);
      }
    }
  }
  
  
  
  // Global escape
  
  global_escape::global_escape (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), rl(reflevel), m(map), c(component)
  {
    assert (! is_meta_mode());
    if (! is_global_mode()) {
      if (! is_level_mode()) {
        if (! is_singlemap_mode()) {
          leave_local_mode (cgh);
        }
        leave_singlemap_mode (cgh);
      }
      leave_level_mode (cgh);
    }
  }
  
  global_escape::~global_escape ()
  {
    assert (is_global_mode());
    if (rl != -1) {
      enter_level_mode (cgh, rl);
      if (m != -1) {
        enter_singlemap_mode (cgh, m);
        if (c != -1) {
          enter_local_mode (cgh, c);
        }
      }
    }
  }
  
  
  
  // Meta escape
  
  meta_escape::meta_escape (cGH const * const cgh_)
    : cgh(const_cast<cGH*>(cgh_)), ml(mglevel), rl(reflevel), m(map), c(component)
  {
    if (! is_meta_mode()) {
      if (! is_global_mode()) {
        if (! is_level_mode()) {
          if (! is_singlemap_mode()) {
            leave_local_mode (cgh);
          }
          leave_singlemap_mode (cgh);
        }
        leave_level_mode (cgh);
      }
      leave_global_mode (cgh);
    }
  }
  
  meta_escape::~meta_escape ()
  {
    assert (is_meta_mode());
    if (ml != -1) {
      enter_global_mode (cgh, ml);
      if (rl != -1) {
        enter_level_mode (cgh, rl);
        if (m != -1) {
          enter_singlemap_mode (cgh, m);
          if (c != -1) {
            enter_local_mode (cgh, c);
          }
        }
      }
    }
  }
  
  
  
  //
  // Call functions in specific modes
  //
  
  int CallLocalFunction (cGH * const cgh,
                         void (* const function) (cGH * const cgh))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        BEGIN_REFLEVEL_LOOP(cgh) {
          BEGIN_MAP_LOOP(cgh, CCTK_GF) {
            BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
              function (cgh);
            } END_LOCAL_COMPONENT_LOOP;
          } END_MAP_LOOP;
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    } else if (is_global_mode()) {
      BEGIN_REFLEVEL_LOOP(cgh) {
        BEGIN_MAP_LOOP(cgh, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
            function (cgh);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } END_REFLEVEL_LOOP;
    } else if (is_level_mode()) {
      BEGIN_MAP_LOOP(cgh, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
          function (cgh);
        } END_LOCAL_COMPONENT_LOOP;
      } END_MAP_LOOP;
    } else if (is_singlemap_mode()) {
      BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
        function (cgh);
      } END_LOCAL_COMPONENT_LOOP;
    } else {
      function (cgh);
    }
    return 0;
  }
  
  int CallSinglemapFunction (cGH * const cgh,
                             void (* const function) (cGH * const cgh))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        BEGIN_REFLEVEL_LOOP(cgh) {
          BEGIN_MAP_LOOP(cgh, CCTK_GF) {
            function (cgh);
          } END_MAP_LOOP;
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    } else if (is_global_mode()) {
      BEGIN_REFLEVEL_LOOP(cgh) {
        BEGIN_MAP_LOOP(cgh, CCTK_GF) {
          function (cgh);
        } END_MAP_LOOP;
      } END_REFLEVEL_LOOP;
    } else if (is_level_mode()) {
      BEGIN_MAP_LOOP(cgh, CCTK_GF) {
        function (cgh);
      } END_MAP_LOOP;
    } else {
      BEGIN_SINGLEMAP_MODE(cgh) {
        function (cgh);
      } END_SINGLEMAP_MODE;
    }
    return 0;
  }
  
  int CallLevelFunction (cGH * const cgh,
                         void (* const function) (cGH * const cgh))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        BEGIN_REFLEVEL_LOOP(cgh) {
          function (cgh);
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    } else if (is_global_mode()) {
      BEGIN_REFLEVEL_LOOP(cgh) {
        function (cgh);
      } END_REFLEVEL_LOOP;
    } else {
      BEGIN_LEVEL_MODE(cgh) {
        function (cgh);
      } END_LEVEL_MODE;
    }
    return 0;
  }
  
  int CallGlobalFunction (cGH * const cgh,
                          void (* const function) (cGH * const cgh))
  {
    if (is_meta_mode()) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        function (cgh);
      } END_MGLEVEL_LOOP;
    } else {
      BEGIN_GLOBAL_MODE(cgh) {
        function (cgh);
      } END_GLOBAL_MODE;
    }
    return 0;
  }
  
  int CallMetaFunction (cGH * const cgh,
                        void (* const function) (cGH * const cgh))
  {
    BEGIN_META_MODE(cgh) {
      function (cgh);
    } END_META_MODE;
    return 0;
  }
  
  
  
  //
  // Call a scheduling group
  //
  
  int CallScheduleGroup (cGH * const cgh, const char * const group)
  {
    CCTK_ScheduleTraverse (group, cgh, CallFunction);
    return 0;
  }
  
  extern "C" void CCTK_FCALL CCTK_FNAME(CallScheduleGroup)
    (int * const ierr, cGH * const * const cgh, ONE_FORTSTRING_ARG)
  {
    ONE_FORTSTRING_CREATE (group);
    *ierr = CallScheduleGroup (*cgh, group);
    free (group);
  }
  
} // namespace Carpet
