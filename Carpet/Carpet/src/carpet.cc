// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Attic/carpet.cc,v 1.5 2001/03/10 20:55:03 eschnett Exp $

/* It is assumed that the number of components of all arrays is equal
   to the number of components of the grid functions, and that their
   distribution onto the processors is the same.  */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_GHExtensions.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_WarnLevel.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/dist.hh"
#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/th.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Attic/carpet.cc,v 1.5 2001/03/10 20:55:03 eschnett Exp $";



namespace Carpet {
  
  static void RecursiveInitialise (cGH* cgh);
  static void RecursiveEvolve (cGH* cgh);
  static void RecursiveShutdown (cGH* cgh);
  
  static void Recompose (cGH* cgh);
  static void CycleTimeLevels (cGH* cgh);
  static void Restrict (cGH* cgh);
  
  // Debugging output
  static void Checkpoint (const char* fmt, ...);
  
  
  
  // handle from CCTK_RegisterGHExtension
  int GHExtension;
  
  // data for scalars
  vector<vector<vector<void*> > > scdata;// [group][var][tl]
  
  // data for arrays
  vector<arrdesc> arrdata;	// [group]
  
  // the grid hierarchy
  gh<dim>* hh;
  th<dim>* tt;
  dh<dim>* dd;
  int gfsize[dim];
  
  // data for grid functions
  vector<gfdesc> gfdata;	// [group]
  
  // active time level
  int activetimelevel;		// 0 for current, 1 for next
  
  // current position on the grid hierarchy
  int mglevel;
  int reflevel;
  int component;
  
  
  
  int CarpetStartup()
  {
    CCTK_RegisterBanner ("AMR driver provided by Carpet");
    
    GHExtension = CCTK_RegisterGHExtension("Carpet");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    CCTK_OverloadInitialise (Initialise);
    CCTK_OverloadEvolve (Evolve);
    CCTK_OverloadShutdown (Shutdown);
    
    CCTK_OverloadSyncGroup (SyncGroup);
    CCTK_OverloadEnableGroupStorage (EnableGroupStorage);
    CCTK_OverloadDisableGroupStorage (DisableGroupStorage); 
    CCTK_OverloadEnableGroupComm (EnableGroupComm);
    CCTK_OverloadDisableGroupComm (DisableGroupComm);
    CCTK_OverloadBarrier (Barrier);
    CCTK_OverloadExit (Exit);
    CCTK_OverloadAbort (Abort);
    CCTK_OverloadMyProc (myProc);
    CCTK_OverloadnProcs (nProcs);
    CCTK_OverloadArrayGroupSizeB (ArrayGroupSizeB);
    CCTK_OverloadQueryGroupStorageB (QueryGroupStorageB);
    
    return 0;
  }
  
  
  
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // not sure what to do with that
    assert (convLevel==0);
    
    dist::pseudoinit();
    Checkpoint ("starting SetupGH...");
    
    // ghost zones
    vect<int,dim> lghosts, ughosts;
    if (ghost_size == -1) {
      lghosts = vect<int,dim>(ghost_size_x, ghost_size_y, ghost_size_z);
      ughosts = vect<int,dim>(ghost_size_x, ghost_size_y, ghost_size_z);
    } else {
      lghosts = vect<int,dim>(ghost_size, ghost_size, ghost_size);
      ughosts = vect<int,dim>(ghost_size, ghost_size, ghost_size);
    }
    
    // grid size
    const int stride
      = (int)floor(pow(refinement_factor, max_refinement_levels) + 0.5);
    vect<int,dim> npoints;
    if (global_nsize == -1) {
      npoints =  vect<int,dim>(global_nx, global_ny, global_nz);
    } else {
      npoints =  vect<int,dim>(global_nsize, global_nsize, global_nsize);
    }
    
    const vect<int,dim> str(stride);
    const vect<int,dim> lb(lghosts * str);
    const vect<int,dim> ub = (npoints - ughosts - 1) * str;
    
    const bbox<int,dim> baseext(lb, ub, str);
    
    // allocate grid hierarchy
    hh = new gh<dim>(refinement_factor, vertex_centered,
		     multigrid_factor, vertex_centered,
		     baseext);
    
    // allocate time hierarchy
    tt = new th<dim>
      (*hh, (int)floor(pow(refinement_factor, max_refinement_levels) + 0.5));
    
    // allocate data hierarchy
    dd = new dh<dim>(*hh, lghosts, ughosts);
    
    // allocate space for groups
    scdata.resize(CCTK_NumGroups());
    arrdata.resize(CCTK_NumGroups());
    gfdata.resize(CCTK_NumGroups());
    
    // allocate space for variables in group
    // (don't enable storage yet)
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      switch (CCTK_GroupTypeI(group)) {
	
      case CCTK_SCALAR: {
	scdata[group].resize(CCTK_NumVarsInGroupI(group));
	for (int var=0; var<(int)scdata[group].size(); ++var) {
	  const int n = (CCTK_FirstVarIndexI(group) + var);
	  scdata[group][var].resize(CCTK_NumTimeLevelsFromVarI(n));
	  for (int tl=0; tl<(int)scdata[group][var].size(); ++tl) {
	    scdata[group][var][tl] = 0;
	  }
	}
	break;
      }
	
      case CCTK_ARRAY: {
	vect<int,dim> alb, aub;
	for (int d=0; d<dim; ++d) {
	  alb[d] = 0;
	  aub[d] = (*CCTK_GroupSizesI(group))[d];
	}
	const bbox<int,dim> arrext(alb, aub-str, str);
	
	arrdata[group].hh = new gh<dim>
	  (refinement_factor, vertex_centered,
	   multigrid_factor, vertex_centered,
	   arrext);
	
	arrdata[group].tt = new th<dim>
	  (*hh,
	   (int)floor(pow(refinement_factor, max_refinement_levels) + 0.5));
	
	vect<int,dim> alghosts, aughosts;
	for (int d=0; d<dim; ++d) {
	  alghosts[d] = (*CCTK_GroupGhostsizesI(group))[d];
	  aughosts[d] = (*CCTK_GroupGhostsizesI(group))[d];
	}
	
	arrdata[group].dd = new dh<dim>(*hh, alghosts, aughosts);
	
	arrdata[group].data.resize(CCTK_NumVarsInGroupI(group));
	for (int var=0; var<(int)scdata[group].size(); ++var) {
	  arrdata[group].data[var] = 0;
	}
	break;
      }
	
      case CCTK_GF: {
	gfdata[group].data.resize(CCTK_NumVarsInGroupI(group));
	for (int var=0; var<(int)scdata[group].size(); ++var) {
	  gfdata[group].data[var] = 0;
	}
	break;
      }
	
      default:
	abort();
      }
    }
    
    // active time level
    activetimelevel = 0;
    
    // current position
    reflevel  = 0;
    mglevel   = 0;
    component = -1;
    
    // initialise cgh
    cgh->cctk_convlevel = mglevel;
    for (int d=0; d<dim; ++d) {
      cgh->cctk_levfac[d] = (int)floor(pow(hh->reffact, reflevel)+0.5);
    }
    
    // Recompose grid hierarchy
    Recompose (cgh); 
    
    Checkpoint ("done with SetupGH.");
    
    // We register only once, ergo we get only one handle, ergo there
    // is only one grid hierarchy for us.  We store that statically,
    // so there is no need to pass it to Cactus.
    return 0;
  }
  
  
  
  int Initialise (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Initialise stuff
    const int convlev = 0;
    cGH* const cgh = CCTK_SetupGH (fc, convlev);
    CCTKi_AddGH (fc, convlev, cgh);
    
    // Delay checkpoint until MPI has been initialised
    Checkpoint ("starting Initialise...");
    
    // Initialise stuff
    cgh->cctk_iteration = 0;
    cgh->cctk_time = cctk_initial_time;
    
    // Enable storage and communtication
    CCTKi_ScheduleGHInit (cgh);
    
    // Initialise stuff
    CCTKi_InitGHExtensions (cgh);
    
    // Check parameters
    CCTK_ScheduleTraverse ("CCTK_PARAMCHECK", cgh, CallFunction);
    CCTKi_FinaliseParamWarn();
    
    // Initialise all levels recursively
    RecursiveInitialise (cgh);
    
    // Output
    CCTK_OutputGH (cgh);
    
    Checkpoint ("done with Initialise.");
    
    return 0;
  }
  
  static void RecursiveInitialise (cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("%*sstarting RecursiveInitialise on level %d...",
		2*reflevel, "", reflevel);
    
    // Set up the grid
    CCTK_ScheduleTraverse ("CCTK_BASEGRID", cgh, CallFunction);
    
    // Set up the initial data
    CCTK_ScheduleTraverse ("CCTK_INITIAL", cgh, CallFunction);
    CCTK_ScheduleTraverse ("CCTK_POSTINITIAL", cgh, CallFunction);
    
    // Recover
    CCTK_ScheduleTraverse ("CCTK_RECOVER_VARIABLES", cgh, CallFunction);
    CCTK_ScheduleTraverse ("CCTK_CPINITIAL", cgh, CallFunction);
    
    // Recurse
    if (reflevel < hh->reflevels()-1) {
      reflevel_up (cgh);
      RecursiveInitialise (cgh);
      reflevel_down (cgh);
    }
    
    // Poststep
    CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
    
    // Analysis
    CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
    
    Checkpoint ("%*sdone with RecursiveInitialise on level %d...",
		2*reflevel, "", reflevel);
  }
  
  
  
  int Evolve (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("starting Evolve...");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    // Main loop
    while (cgh->cctk_iteration < cctk_itlast
	   || (cctk_final_time >= cctk_initial_time
	       && cgh->cctk_time < cctk_final_time)) {
      
      // Next iteration
      ++cgh->cctk_iteration;
      
      Checkpoint ("Evolving iteration %d...", cgh->cctk_iteration);
      
      RecursiveEvolve (cgh);
      
      // Output
      CCTK_OutputGH (cgh);
    }
    
    Checkpoint ("done with Evolve.");
    
    return 0;
  }
  
  static void RecursiveEvolve (cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("%*sstarting RecursiveEvolve on level %d...",
		2*reflevel, "", reflevel);
    
    // Move activity to next time level
    assert (activetimelevel == 0);
    activetimelevel = 1;
    
    // Prestep
    CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
    
    // Evolve
    CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
    
    // Recurse
    if (reflevel < hh->reflevels()-1) {
      reflevel_up (cgh);
      for (int i=0; i<hh->reffact; ++i) {
	RecursiveEvolve (cgh);
      }
      reflevel_down (cgh);
    } else {
      cgh->cctk_time += cgh->cctk_delta_time;
    }
    
    // Restrict
    Restrict (cgh);
    
    // Cycle time levels
    CycleTimeLevels (cgh);
    
    // Move activity back to current time level
    assert (activetimelevel == 1);
    activetimelevel = 0;
    
    // Poststep
    CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
    
    // Checkpoint
    CCTK_ScheduleTraverse ("CCTK_CHECKPOINT", cgh, CallFunction);
    
    // Analysis
    CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
    
    Checkpoint ("%*sdone with RecursiveEvolve on level %d...",
		2*reflevel, "", reflevel);
  }
  
  
  
  int Shutdown (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("starting Shutdown...");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    RecursiveShutdown (cgh);
    
    CCTK_PRINTSEPARATOR;
    printf ("Done.\n");
    
    // earlier checkpoint before calling finalising MPI
    Checkpoint ("done with Shutdown.");
    
    dist::finalize();
    
    return 0;
  }
  
  static void RecursiveShutdown (cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("%*sstarting RecursiveShutdown on level %d...",
		2*reflevel, "", reflevel);
    
    // Recurse
    if (reflevel < hh->reflevels()-1) {
      reflevel_up (cgh);
      RecursiveShutdown (cgh);
      reflevel_down (cgh);
    }
    
    // Terminate
    CCTK_ScheduleTraverse ("CCTK_TERMINATE", cgh, CallFunction);
    
    // Shutdown
    CCTK_ScheduleTraverse ("CCTK_SHUTDOWN", cgh, CallFunction);
    
    Checkpoint ("%*sdone with RecursiveShutdown on level %d...",
		2*reflevel, "", reflevel);
  }
  
  
  
  int ScheduleTraverse (cGH* cgh, const char* rfrPoint)
  {
    // traverse all functions on all refinement levels on all
    // multigrid levels
    
    abort();
    
    assert (mglevel==-1);
    for (mglevel=0; mglevel<hh->mglevels(0,0); ++mglevel) {
      assert (reflevel==-1);
      for (reflevel=0; reflevel<hh->reflevels(); ++reflevel) {
	CCTK_ScheduleTraverse (rfrPoint, cgh, CallFunction);
      }
      reflevel = -1;
    }
    mglevel = -1;
    
    return 0;
  }
  
  
  
  int CallFunction (void* function, cFunctionData* attribute, void* data)
  {
    // traverse one function on all components of one refinement level
    // of one multigrid level
    
    assert (mglevel>=0);
    assert (reflevel>=0);
    
//     Checkpoint ("%*sstarting CallFunction...", 2*reflevel, "");
    
    cGH* cgh = (cGH*)data;
    
    // set Cactus parameters
    for (int d=0; d<dim; ++d) {
      cgh->cctk_nghostzones[d] = dd->lghosts[d];
    }
    
    if (attribute->global) {
      // global operation: call once per refinement level
      
      assert (component==-1);
      
      // set Cactus parameters to pseudo values
      for (int d=0; d<dim; ++d) {
	cgh->cctk_lsh[d]      = 0xdeadbeef;
	cgh->cctk_gsh[d]      = 0xdeadbeef;
	cgh->cctk_bbox[2*d  ] = 0xdeadbeef;
	cgh->cctk_bbox[2*d+1] = 0xdeadbeef;
	cgh->cctk_lbnd[d]     = 0xdeadbeef;
	cgh->cctk_ubnd[d]     = 0xdeadbeef;
	for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	  cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = 0xdeadbeef;
	}
      }
      
      // set local grid function and array sizes
      for (int d=0; d<dim; ++d) {
	gfsize[d] = 0xdeadbeef;
      }
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	const int n0 = CCTK_FirstVarIndexI(group);
	switch (CCTK_GroupTypeFromVarI(n0)) {
	case CCTK_SCALAR:
	  break;
	case CCTK_ARRAY:
	  for (int d=0; d<dim; ++d) {
	    arrdata[group].size[d] = 0xdeadbeef;
	  }
	  break;
	case CCTK_GF:
	  break;
	default:
	  abort();
	}
      }
      
      // set Cactus pointers to data
      for (int n=0; n<CCTK_NumVars(); ++n) {
	for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n); ++tl) {
	  
	  const int group = CCTK_GroupIndexFromVarI(n);
	  assert (group>=0);
	  const int var   = n - CCTK_FirstVarIndexI(group);
	  assert (var>=0);
	  
	  if (CCTK_QueryGroupStorageI(cgh, group)) {
	    // group has storage
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR:
	      // scalar variables can be accessed
	      assert (group<(int)scdata.size());
	      assert (var<(int)scdata[group].size());
	      assert (tl<(int)scdata[group][var].size());
	      cgh->data[n][tl] = scdata[group][var][tl];
	      break;
	      
	    case CCTK_ARRAY:
	    case CCTK_GF:
	      // arrays and grid functions cannot be accessed
	      cgh->data[n][tl] = 0;
	      break;
	      
	    default:
	      abort();
	    }
	    
	  } else {
	    
	    // group has no storage
	    cgh->data[n][tl] = 0;
	    
	  }
	  
	} // for tl
      }	// for n
      
      // traverse
      const int res = CCTK_CallFunction (function, attribute, data);
      assert (res==0);
      
    } else {
      // local operation: call once per component
      
      assert (component==-1);
      for (component=0; component<hh->components(reflevel); ++component) {
	
	// maybe: traverse only if the component is local to the processor
	// maybe not, because arrays might have different distribution
	// than grid functions
	
	// this requires that all processors have the same number of
	// local components
 	if (hh->is_local(reflevel, component)) {
	  
	  // set Cactus parameters to pseudo values
	  for (int d=0; d<dim; ++d) {
	    typedef bbox<int,dim> ibbox;
	    const ibbox ext = dd->boxes[reflevel][component][mglevel].exterior;
	    const ibbox base = hh->baseextent;
	    cgh->cctk_lsh[d] = (ext.shape() / ext.stride())[d];
	    cgh->cctk_gsh[d]
	      = ((base.shape() / base.stride() + dd->lghosts + dd->ughosts)[d]
		 * cgh->cctk_levfac[d]);
	    cgh->cctk_lbnd[d]     = (ext.lower() / ext.stride())[d];
	    cgh->cctk_ubnd[d]     = (ext.upper() / ext.stride())[d];
	    cgh->cctk_bbox[2*d  ] = cgh->cctk_lbnd[d] == 0;
	    cgh->cctk_bbox[2*d+1] = cgh->cctk_ubnd[d] == cgh->cctk_gsh[d]-1;
	    for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	      // TODO: support staggering
	      cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
	    }
	  }
	  
	  // set local grid function and array sizes
#if 0
	  if (hh->is_local(reflevel, component)) {
#endif
	    const bbox<int,dim> ext =
	      dd->boxes[reflevel][component][mglevel].exterior;
	    for (int d=0; d<dim; ++d) {
	      gfsize[d] = ext.shape()[d] / ext.stride()[d];
	    }
#if 0
	  } else {
	    for (int d=0; d<dim; ++d) {
	      gfsize[d] = 0;
	    }
	  }
#endif
	  
	  for (int group=0; group<CCTK_NumGroups(); ++group) {
	    const int n0 = CCTK_FirstVarIndexI(group);
	    switch (CCTK_GroupTypeFromVarI(n0)) {
	    case CCTK_SCALAR:
	      break;
	    case CCTK_ARRAY:
	      if (arrdata[group].hh->is_local(reflevel, component)) {
		const bbox<int,dim> ext =
		  arrdata[group].dd->
		  boxes[reflevel][component][mglevel].exterior;
		for (int d=0; d<dim; ++d) {
		  arrdata[group].size[d] = ext.shape()[d] / ext.stride()[d];
		}
	      } else {
		for (int d=0; d<dim; ++d) {
		  arrdata[group].size[d] = 0;
		}
	      }
	      break;
	    case CCTK_GF:
	      break;
	    default:
	      abort();
	    }
	  }
	  
	  // set Cactus pointers to data
	  for (int n=0; n<CCTK_NumVars(); ++n) {
	    for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n); ++tl) {
	      
	      const int group = CCTK_GroupIndexFromVarI(n);
	      assert (group>=0);
	      const int var   = n - CCTK_FirstVarIndexI(group);
	      assert (var>=0);
	      
	      if (CCTK_QueryGroupStorageI(cgh, group)) {
		// group has storage
		
		switch (CCTK_GroupTypeFromVarI(n)) {
		  
		case CCTK_SCALAR:
		  assert (group<(int)scdata.size());
		  assert (var<(int)scdata[group].size());
		  assert (tl<(int)scdata[group][var].size());
		  cgh->data[n][tl] = scdata[group][var][tl];
		  break;
		  
		case CCTK_ARRAY:
		  assert (group<(int)arrdata.size());
		  assert (var<(int)arrdata[group].data.size());
		  cgh->data[n][tl] =
		    ((*arrdata[group].data[var])(tl,reflevel,component,mglevel)
		     ->storage());
		  break;
		  
		case CCTK_GF:
		  assert (group<(int)gfdata.size());
		  assert (var<(int)gfdata[group].data.size());
		  cgh->data[n][tl] =
		    ((*gfdata[group].data[var])(tl,reflevel,component,mglevel)
		     ->storage());
		  break;
		  
		default:
		  abort();
		}
		assert (cgh->data[n][tl]);
		
	      } else {
		
		// group has no storage
		cgh->data[n][tl] = 0;
		
	      }	// if ! has storage
	      
	    } // for tl
	  } // for n
	  
	  const int res = CCTK_CallFunction (function, attribute, data);
	  assert (res==0);
	  
	} // if is_local
	
      }	// for component
      component = -1;
      
    } // if local operation
    
//     Checkpoint ("%*sdone with CallFunction.", 2*reflevel, "");
    
#if 0
    // return 0: let the flesh do the synchronisation, if necessary
    return 0;
#endif
    
    // synchronise, because our bbox information was wrong
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      SyncGroup (cgh, CCTK_GroupName(group));
    }
    
    // return 1: we did synchronise
    return 1;
  }
  
  
  
  void reflevel_up (cGH* cgh)
  {
    // Check
    assert (reflevel < hh->reflevels()-1);
    assert (mglevel == 0);
    
    // Change
    cgh->cctk_delta_time /= hh->reffact;
    ++reflevel;
  }
  
  void reflevel_down (cGH* cgh)
  {
    // Check
    assert (reflevel > 0);
    assert (mglevel == 0);
    
    // Change
    cgh->cctk_delta_time *= hh->reffact;
    --reflevel;
  }
  
  
  
  int SyncGroup (cGH* cgh, const char* groupname)
  {
    assert (component == -1);
    
    Checkpoint ("SyncGroup %s", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot synchronise group \"%s\" because it has no storage",
		  groupname);
      return -1;
    }
    
    const int n0 = CCTK_FirstVarIndexI(group);
    const int tl = max(0, (CCTK_NumTimeLevelsFromVarI(n0)
			   + activetimelevel - 2));
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      break;
      
    case CCTK_ARRAY:
      assert (group<(int)arrdata.size());
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	if (reflevel>0) {
	  for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	    arrdata[group].data[var]->ref_bnd_prolongate
	      (tl, reflevel, c, mglevel);
	  }
	}
	for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	  arrdata[group].data[var]->sync (tl, reflevel, c, mglevel);
	}
      }
      break;
      
    case CCTK_GF:
      assert (group<(int)gfdata.size());
      for (int var=0; var<(int)gfdata[group].data.size(); ++var) {
	if (reflevel>0) {
	  for (int c=0; c<hh->components(reflevel); ++c) {
	    gfdata[group].data[var]->ref_bnd_prolongate
	      (tl, reflevel, c, mglevel);
	  }
	}
	for (int c=0; c<hh->components(reflevel); ++c) {
	  gfdata[group].data[var]->sync (tl, reflevel, c, mglevel);
	}
      }
      break;
      
    default:
      abort();
    }
    
    return 0;
  }
  
  
  
  int EnableGroupStorage (cGH* cgh, const char* groupname)
  {
    Checkpoint ("EnableGroupStorage %s", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    // The return values seems to be whether storage was enabled
    // previously.
    const int retval = CCTK_QueryGroupStorageI (cgh, group);
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      if (scdata[group].size()==0 || scdata[group][0].size()==0
	  || scdata[group][0][0] != 0) {
	// group already has storage
	break;
      }
      for (int var=0; var<(int)scdata[group].size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	for (int tl=0; tl<(int)scdata[group][tl].size(); ++tl) {
	  switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)				\
	  case N:				\
	    scdata[group][var][tl] = new T;	\
	    break;
#include "typecase"
#undef TYPECASE
	  default:
	    abort();
	  }
	}
      }
      break;
      
    case CCTK_ARRAY:
      if (arrdata[group].data.size()==0
	  || arrdata[group].data[0] != 0) {
	// group already has storage
	break;
      }
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)							\
	case N:								\
	  arrdata[group].data[var] = new gf<T,dim>			\
	    (CCTK_VarName(n), *arrdata[group].tt, *arrdata[group].dd,	\
	     0, CCTK_NumTimeLevelsFromVarI(n)-1);			\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  abort();
	}
      }
      break;
      
    case CCTK_GF:
      if (gfdata[group].data.size()==0
	  || gfdata[group].data[0] != 0) {
	// group already has storage
	break;
      }
      for (int var=0; var<(int)gfdata[group].data.size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	case N:						\
	  gfdata[group].data[var] = new gf<T,dim>	\
	    (CCTK_VarName(n), *tt, *dd,			\
	     0, CCTK_NumTimeLevelsFromVarI(n)-1);	\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  abort();
	}
      }
      break;
      
    default:
      abort();
    }
    
    // The return value seems to be 1 (success) no matter whether
    // storage has actually been disabled.
    return retval;
  }
  
  
  
  int DisableGroupStorage (cGH* cgh, const char* groupname)
  {
    Checkpoint ("DisableGroupStorage %s", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      if (scdata[group].size()==0 || scdata[group][0].size()==0
	  || scdata[group][0][0] == 0) {
	// group already has no storage
	break;
      }
      for (int var=0; var<(int)scdata[group].size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	for (int tl=0; tl<(int)scdata[group][tl].size(); ++tl) {
	  switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)				\
	  case N:				\
	    delete (T*)scdata[group][var][tl];	\
	    break;
#include "typecase"
#undef TYPECASE
	  default:
	    abort();
	  }
	  scdata[group][var][tl] = 0;
	}
      }
      break;
      
    case CCTK_ARRAY:
      if (arrdata[group].data.size()==0 || arrdata[group].data[0] == 0) {
	// group already has no storage
	break;
      }
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	case N:						\
	  delete (gf<T,dim>*)arrdata[group].data[var];	\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  abort();
	}
	arrdata[group].data[var] = 0;
      }
      break;
      
    case CCTK_GF:
      if (gfdata[group].data.size()==0
	  || gfdata[group].data[0] == 0) {
	// group already has no storage
	break;
      }
      for (int var=0; var<(int)gfdata[group].data.size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	case N:						\
	  delete (gf<T,dim>*)gfdata[group].data[var];	\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  abort();
	}
	gfdata[group].data[var] = 0;
      }
      break;
      
    default:
      abort();
    }
    
    // The return value seems to be 1 (success) no matter whether
    // storage has actually been disabled.
    return 1;
  }
  
  
  
  int EnableGroupComm (cGH* cgh, const char* groupname)
  {
    // Communication is always enabled
    return 0;
  }
  
  int DisableGroupComm (cGH* cgh, const char* groupname)
  {
    // Communication is always enabled
    return -1;
  }
  
  
  
  static void Recompose (cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (component == -1);
    Checkpoint ("%*sRecompose", 2*reflevel, "");
    
    const int nprocs    = CCTK_nProcs(cgh);
    const int reflevels = max_refinement_levels; // arbitrary value
    const int mglevels  = 1;	// arbitrary value
    vector<vector<bbox<int,dim> > > bbss(reflevels);
    // note: what this routine calls "ub" is "ub+str" elsewhere
    vect<int,dim> rstr = hh->baseextent.stride();
    vect<int,dim> rlb  = hh->baseextent.lower();
    vect<int,dim> rub  = hh->baseextent.upper() + rstr;
    for (int rl=0; rl<reflevels; ++rl) {
      if (rl>0) {
	// refined boxes have smaller stride
	assert (all(rstr%hh->reffact == 0));
	rstr /= hh->reffact;
	// refine (arbitrarily) the center only
	rlb /= hh->reffact;
	rub /= hh->reffact;
      }
      vector<bbox<int,dim> > bbs(nprocs);
      for (int c=0; c<nprocs; ++c) {
	vect<int,dim> cstr = rstr;
	vect<int,dim> clb = rlb;
	vect<int,dim> cub = rub;
	const int glonpz = (rub[dim-1] - rlb[dim-1]) / cstr[dim-1];
	const int locnpz = (glonpz + nprocs - 1) / nprocs;
	const int zstep = locnpz * cstr[dim-1];
	clb[dim-1] = rlb[dim-1] + zstep *  c;
	cub[dim-1] = rlb[dim-1] + zstep * (c+1);
	if (c == nprocs-1) cub[dim-1] = rub[dim-1];
	bbs[c] = bbox<int,dim>(clb, cub-cstr, cstr);
      }
      bbss[rl] = bbs;
    }
    vector<vector<vector<bbox<int,dim> > > > bbsss
      = hh->make_multigrid_boxes(bbss, mglevels);
    vector<vector<int> > pss(bbss.size());
    for (int rl=0; rl<reflevels; ++rl) {
      pss[rl] = vector<int>(bbss[rl].size());
      // make sure all processors have the same number of components
      assert (bbss[rl].size() % nprocs == 0);
      for (int c=0; c<(int)bbss[rl].size(); ++c) {
	pss[rl][c] = c % nprocs; // distribute among processors
      }
    }
    hh->recompose(bbsss, pss);
  }
  
  
  
  static void CycleTimeLevels (cGH* cgh)
  {
    Checkpoint ("%*sCycleTimeLevels", 2*reflevel, "");
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  const int n = CCTK_FirstVarIndexI(group) + var;
	  switch (CCTK_GroupTypeFromVarI(n)) {
	  case CCTK_SCALAR: {
	    assert (group<(int)scdata.size());
	    assert (var<(int)scdata[group].size());
	    const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n));
	    for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n)-1; ++tl) {
	      memcpy (scdata[group][var][tl], scdata[group][var][tl+1], sz);
	    }
	    break;
	  }
	  case CCTK_ARRAY: {
	    assert (group<(int)arrdata.size());
	    assert (var<(int)arrdata[group].data.size());
	    for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	      for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n)-1; ++tl) {
		arrdata[group].data[var]->copy(tl, reflevel, c, mglevel);
	      }
	    }
	    break;
	  }
	  case CCTK_GF: {
	    assert (group<(int)gfdata.size());
	    assert (var<(int)gfdata[group].data.size());
	    for (int c=0; c<hh->components(reflevel); ++c) {
	      for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n)-1; ++tl) {
		gfdata[group].data[var]->copy(tl, reflevel, c, mglevel);
	      }
	    }
	    break;
	  }
	  default:
	    abort();
	  }
	}
      }
    }
  }
  
  
  
  static void Restrict (cGH* cgh)
  {
    assert (component == -1);
    
    Checkpoint ("%*sRestrict", 2*reflevel, "");
    
    if (reflevel == hh->reflevels()-1) return;
    
    for (component=0; component<hh->components(reflevel); ++component) {
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	
	const int n0 = CCTK_FirstVarIndexI(group);
	const int tl = max(0, (CCTK_NumTimeLevelsFromVarI(n0)
			       + activetimelevel - 2));
	
	switch (CCTK_GroupTypeI(group)) {
	  
	case CCTK_SCALAR:
	  break;
	  
	case CCTK_ARRAY:
	  for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	    arrdata[group].data[var]->ref_restrict
	      (tl, reflevel, component, mglevel);
	  }
	  break;
	  
	case CCTK_GF:
	  for (int var=0; var<(int)gfdata[group].data.size(); ++var) {
	    gfdata[group].data[var]->ref_restrict
	      (tl, reflevel, component, mglevel);
	  }
	  break;
	  
	default:
	  abort();
	}
      }
    }
    component = -1;
  }
  
  
  
  int Barrier (cGH* cgh)
  {
    MPI_Barrier (dist::comm);
    return 0;
  }
  
  
  
  int Exit (cGH* cgh, int retval)
  {
    CCTK_Barrier (cgh);
    dist::finalize();
    exit (retval);
  }
  
  int Abort (cGH* cgh, int retval)
  {
    MPI_Abort (dist::comm, retval);
    abort ();
  }
  
  
  
  int myProc (cGH* cgh)
  {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    return rank;
  }
  
  int nProcs (cGH* cgh)
  {
    int size;
    MPI_Comm_size (dist::comm, &size);
    return size;
  }
  
  
  
  const int* ArrayGroupSizeB (cGH* cgh, int dir, int group,
			      const char* groupname)
  {
    static const int zero = 0, one = 1;
    
    if (component == -1) {
      // global routine
      return &zero;
    }
    
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    assert (dir>=0 && dir<dim);
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      
      const int var = CCTK_FirstVarIndexI(group);
      assert (var>=0 && var<CCTK_NumVars());
      
      switch (CCTK_GroupTypeFromVarI(var)) {
	
      case CCTK_SCALAR:
	assert (group<(int)scdata.size());
	return &one;
	
      case CCTK_ARRAY:
	assert (group<(int)arrdata.size());
	return &arrdata[group].size[dir];
	
      case CCTK_GF:
	assert (group<(int)gfdata.size());
	return &gfsize[dir];
	
      default:
	abort();
      }
      
    } else {
      
      // no storage
      return &zero;
      
    }
  }
  
  
  
  int QueryGroupStorageB (cGH* cgh, int group, const char* groupname)
  {
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    const int n = CCTK_FirstVarIndexI(group);
    assert (n>=0 && n<CCTK_NumVars());
    const int var = 0;
    
    switch (CCTK_GroupTypeFromVarI(n)) {
      
    case CCTK_SCALAR: {
      assert (group<(int)scdata.size());
      assert (var<(int)scdata[group].size());
      const int tl=0;
      assert (tl<(int)scdata[group][var].size());
      return scdata[group][var][tl] != 0;
    }
      
    case CCTK_ARRAY: {
      assert (group<(int)arrdata.size());
      assert (var<(int)arrdata[group].data.size());
      return arrdata[group].data[var] != 0;
    }
      
    case CCTK_GF: {
      assert (group<(int)gfdata.size());
      assert (var<(int)gfdata[group].data.size());
      return gfdata[group].data[var] != 0;
    }
      
    default:
      abort();
    }
  }
  
  
  
  void Checkpoint (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    if (verbose) {
      va_list args;
      char msg[1000];
      va_start (args, fmt);
      vsnprintf (msg, sizeof(msg), fmt, args);
      va_end (args);
      CCTK_INFO (msg);
    }
    if (barriers) {
      MPI_Barrier (dist::comm);
    }
  }
  
  
  
  MPI_Comm CarpetMPICommunicator ()
  {
    return dist::comm;
  }
  
} // namespace Carpet
