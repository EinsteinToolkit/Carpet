// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Attic/carpet.cc,v 1.13 2001/03/16 21:32:17 eschnett Exp $

// It is assumed that the number of components of all arrays is equal
// to the number of components of the grid functions, and that their
// distribution onto the processors is the same, and that all
// processors own the same number of components.

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

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Attic/carpet.cc,v 1.13 2001/03/16 21:32:17 eschnett Exp $";



namespace Carpet {
  
  static void Recompose (cGH* cgh);
  static void CycleTimeLevels (cGH* cgh);
  static void Restrict (cGH* cgh);
  
  // Debugging output
  static void Checkpoint (const char* fmt, ...);
  
  
  
  // Handle from CCTK_RegisterGHExtension
  int GHExtension;
  
  // Maximum number of refinement levels
  int maxreflevels;
  
  // Refinement factor on finest grid
  int maxreflevelfact;
  
  // Active time level
  int activetimelevel;		// 0 for current, 1 for next
  
  // Current iteration per refinement level
  vector<int> iteration;
  
  // Current position on the grid hierarchy
  int mglevel;
  int reflevel;
  int component;
  
  // refinement factor of current level: pow(refinement_factor, reflevel)
  int reflevelfact;
  
  // Time step on base grid
  CCTK_REAL base_delta_time;
  
  
  
  // Data for scalars
  vector<vector<vector<void*> > > scdata; // [group][var][ti]
  
  // Data for arrays
  // TODO: have replicated arrays
  vector<arrdesc> arrdata;	// [group]
  
  // The grid hierarchy
  gh<dim>* hh;
  th<dim>* tt;
  dh<dim>* dd;
  int gfsize[dim];
  
  // Data for grid functions
  vector<gfdesc> gfdata;	// [group]
  
  
  
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
    CCTK_OverloadMyProc (MyProc);
    CCTK_OverloadnProcs (nProcs);
    CCTK_OverloadArrayGroupSizeB (ArrayGroupSizeB);
    CCTK_OverloadQueryGroupStorageB (QueryGroupStorageB);
    
    return 0;
  }
  
  
  
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Not sure what to do with that
    assert (convLevel==0);
    
    dist::pseudoinit();
    Checkpoint ("starting SetupGH...");
    
    // Refinement information
    maxreflevels = max_refinement_levels;
    maxreflevelfact = (int)floor(pow(refinement_factor, maxreflevels-1) + 0.5);
    
    // Ghost zones
    vect<int,dim> lghosts, ughosts;
    if (ghost_size == -1) {
      lghosts = vect<int,dim>(ghost_size_x, ghost_size_y, ghost_size_z);
      ughosts = vect<int,dim>(ghost_size_x, ghost_size_y, ghost_size_z);
    } else {
      lghosts = vect<int,dim>(ghost_size, ghost_size, ghost_size);
      ughosts = vect<int,dim>(ghost_size, ghost_size, ghost_size);
    }
    
    // Grid size
    const int stride = maxreflevelfact;
    vect<int,dim> npoints;
    if (global_nsize == -1) {
      npoints = vect<int,dim>(global_nx, global_ny, global_nz);
    } else {
      npoints = vect<int,dim>(global_nsize, global_nsize, global_nsize);
    }
    
    const vect<int,dim> str(stride);
    const vect<int,dim> lb(lghosts * str);
    const vect<int,dim> ub = (npoints - ughosts - 1) * str;
    
    const bbox<int,dim> baseext(lb, ub, str);
    
    // Allocate grid hierarchy
    hh = new gh<dim>(refinement_factor, vertex_centered,
		     multigrid_factor, vertex_centered,
		     baseext);
    
    // Allocate time hierarchy
    tt = new th<dim>(*hh, maxreflevelfact);
    
    // Allocate data hierarchy
    dd = new dh<dim>(*hh, lghosts, ughosts);
    
    // Allocate space for groups
    scdata.resize(CCTK_NumGroups());
    arrdata.resize(CCTK_NumGroups());
    gfdata.resize(CCTK_NumGroups());
    
    // Allocate space for variables in group (but don't enable storage
    // yet)
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      switch (CCTK_GroupTypeI(group)) {
	
      case CCTK_SCALAR: {
	scdata[group].resize(CCTK_NumVarsInGroupI(group));
	for (int var=0; var<(int)scdata[group].size(); ++var) {
	  const int n = (CCTK_FirstVarIndexI(group) + var);
	  scdata[group][var].resize(CCTK_NumTimeLevelsFromVarI(n));
	  for (int ti=0; ti<(int)scdata[group][var].size(); ++ti) {
	    scdata[group][var][ti] = 0;
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
	   (int)floor(pow(refinement_factor, max_refinement_levels-1) + 0.5));
	
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
    
    // Initialise cgh
    for (int d=0; d<dim; ++d) {
      cgh->cctk_nghostzones[d] = dd->lghosts[d];
    }
    
    // Initialise current position
    reflevel  = 0;
    mglevel   = 0;
    component = -1;
    
    // Recompose grid hierarchy
    Recompose (cgh); 
    
    // Initialise time step on coarse grid
    base_delta_time = 0;
    
    // Active time level
    activetimelevel = 0;
    
    // Current iteration
    iteration.resize(maxreflevels, 0);
    
    // Set current position (this time for real)
    set_reflevel  (cgh, 0);
    set_mglevel   (cgh, 0);
    set_component (cgh, -1);
    
    Checkpoint ("done with SetupGH.");
    
    // We register only once, ergo we get only one handle, ergo there
    // is only one grid hierarchy for us.  We store that statically,
    // so there is no need to pass anything to Cactus.
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
    
    BEGIN_REFLEVEL_LOOP(cgh) {
      
      // Set up the grid
      CCTK_ScheduleTraverse ("CCTK_BASEGRID", cgh, CallFunction);
      if (reflevel==0) {
	base_delta_time = cgh->cctk_delta_time;
      } else {
// 	assert (abs(cgh->cctk_delta_time - base_delta_time / reflevelfactor)
// 		< 1e-6 * base_delta_time);
	// This fixes a bug in CactusBase/Time
	cgh->cctk_delta_time = base_delta_time / reflevelfact;
      }
      
      // Set up the initial data
      CCTK_ScheduleTraverse ("CCTK_INITIAL", cgh, CallFunction);
      CCTK_ScheduleTraverse ("CCTK_POSTINITIAL", cgh, CallFunction);
      
      // Recover
      CCTK_ScheduleTraverse ("CCTK_RECOVER_VARIABLES", cgh, CallFunction);
      CCTK_ScheduleTraverse ("CCTK_CPINITIAL", cgh, CallFunction);
      
    } END_REFLEVEL_LOOP(cgh);
    
    BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
      
      // Restrict
      Restrict (cgh);
      
    } END_REVERSE_REFLEVEL_LOOP(cgh);
    
    BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
      
      // Poststep
      CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
      
      // Analysis
      CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
      
      // Output
      CCTK_OutputGH (cgh);
      
    } END_REVERSE_REFLEVEL_LOOP(cgh);
    
    Checkpoint ("done with Initialise.");
    
    return 0;
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
      
      BEGIN_REFLEVEL_LOOP(cgh) {
	if ((cgh->cctk_iteration-1) % (maxreflevelfact/reflevelfact) == 0) {
	  
	  // Prestep
	  CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	  
	  // Move activity to next time level
	  assert (activetimelevel == 0);
	  ++activetimelevel;
	  
	  // Evolve
	  CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	  
	  // Move activity back to current time level
	  --activetimelevel;
	  assert (activetimelevel == 0);
	  
	}
      } END_REFLEVEL_LOOP(cgh);
      
      // Advance time
      cgh->cctk_time += base_delta_time / maxreflevelfact;
      
      BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
	if (cgh->cctk_iteration % (maxreflevelfact/reflevelfact) == 0) {
	  
	  // Cycle time levels
	  CycleTimeLevels (cgh);
	  
	  // Advance level times
	  tt->advance_time (reflevel, mglevel);
	  for (int group=0; group<CCTK_NumGroups(); ++group) {
	    switch (CCTK_GroupTypeI(group)) {
	    case CCTK_SCALAR:
	      break;
	    case CCTK_ARRAY:
	      arrdata[group].tt->advance_time (reflevel, mglevel);
	      break;
	    case CCTK_GF:
	      break;
	    default:
	      abort();
	    }
	  }
	  
	  // Restrict
	  Restrict (cgh);
	  
	}
      } END_REVERSE_REFLEVEL_LOOP(cgh);
      
      BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
	if (cgh->cctk_iteration % (maxreflevelfact/reflevelfact) == 0) {
	  
	  // Poststep
	  CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	  
	  // Checkpoint
	  CCTK_ScheduleTraverse ("CCTK_CHECKPOINT", cgh, CallFunction);
	  
	  // Analysis
	  CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
	  
	  // Output
	  CCTK_OutputGH (cgh);
	  
	}
      } END_REVERSE_REFLEVEL_LOOP(cgh);
      
    } // main loop
    
    Checkpoint ("done with Evolve.");
    
    return 0;
  }
  
  
  
  int Shutdown (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("starting Shutdown...");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    // Terminate
    BEGIN_REFLEVEL_LOOP(cgh) {
      CCTK_ScheduleTraverse ("CCTK_TERMINATE", cgh, CallFunction);
    } END_REFLEVEL_LOOP(cgh);
    
    // Shutdown
    BEGIN_REFLEVEL_LOOP(cgh) {
      CCTK_ScheduleTraverse ("CCTK_SHUTDOWN", cgh, CallFunction);
    } END_REFLEVEL_LOOP(cgh);
    
    CCTK_PRINTSEPARATOR;
    printf ("Done.\n");
    
    // earlier checkpoint before calling finalising MPI
    Checkpoint ("done with Shutdown.");
    
    dist::finalize();
    
    return 0;
  }
  
  
  
  int CallFunction (void* function, cFunctionData* attribute, void* data)
  {
    // Traverse one function on all components of one refinement level
    // of one multigrid level
    
    assert (mglevel>=0);
    assert (reflevel>=0);
    
//     Checkpoint ("%*sstarting CallFunction...", 2*reflevel, "");
    
    cGH* cgh = (cGH*)data;
    
    if (attribute->global) {
      // Global operation: call once per refinement level
      
      const int res = CCTK_CallFunction (function, attribute, data);
      assert (res==0);
      
    } else {
      // Local operation: call once per component
      
      BEGIN_COMPONENT_LOOP(cgh) {
	// This requires that all processors have the same number of
	// local components
 	if (hh->is_local(reflevel, component)) {
	  
	  const int res = CCTK_CallFunction (function, attribute, data);
	  assert (res==0);
	  
	} // if is_local
	
      }	END_COMPONENT_LOOP(cgh);
      
    }
    
//     Checkpoint ("%*sdone with CallFunction.", 2*reflevel, "");
    
    // return 0: let the flesh do the synchronisation, if necessary
    return 0;
  }
  
  
  
  int SyncGroup (cGH* cgh, const char* groupname)
  {
    assert (component == -1);
    
    Checkpoint ("%*sSyncGroup %s", 2*reflevel, "", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot synchronise group \"%s\" because it has no storage",
		  groupname);
      return -1;
    }
    
    const int n0 = CCTK_FirstVarIndexI(group);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    const int tl = activetimelevel;
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      break;
      
    case CCTK_ARRAY:
      assert (group<(int)arrdata.size());
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	if (num_tl>1 && reflevel>0) {
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
	if (num_tl>1 && reflevel>0) {
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
    Checkpoint ("%*sEnableGroupStorage %s", 2*reflevel, "", groupname);
    
    // TODO: Enabling storage for one refinement level has to enable
    // it for all other refinement levels as well.  Disabling must
    // wait until all refinement levels have been disabled.
    
    // TODO: Invent a mode "reflevel==-1" that is global, i. e. has
    // effect for all refinement levels.  This mode is used during
    // INITIAL, and en-/disabling storage in it is also global.
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    // The return values seems to be whether storage was enabled
    // previously.
    const int retval = CCTK_QueryGroupStorageI (cgh, group);
    
    // There is a difference between the Cactus time levels and the
    // Carpet time levels.  If there are n time levels, then the
    // Cactus time levels are numbered 0 ... n-1, with the current
    // time level being max(0,n-2).  In Carpet, the time levels are
    // numbered 1-max(1,n-1) ... min(1,n-1), where the current time
    // level is always 0.
    const int n0 = CCTK_FirstVarIndexI(group);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    const int tmin = min(0, 2 - num_tl);
    const int tmax = tmin + num_tl - 1;
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      if (scdata[group].size()==0 || scdata[group][0].size()==0
	  || scdata[group][0][0] != 0) {
	// group already has storage
	break;
      }
      for (int var=0; var<(int)scdata[group].size(); ++var) {
	for (int ti=0; ti<(int)scdata[group][var].size(); ++ti) {
	  const int n = n0 + var;
	  switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)				\
	  case N:				\
	    scdata[group][var][ti] = new T;	\
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
	const int n = n0 + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)							\
	case N:								\
	  arrdata[group].data[var] = new gf<T,dim>			\
	    (CCTK_VarName(n), *arrdata[group].tt, *arrdata[group].dd,	\
	     tmin, tmax);						\
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
	const int n = n0 + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	case N:						\
	  gfdata[group].data[var] = new gf<T,dim>	\
	    (CCTK_VarName(n), *tt, *dd, tmin, tmax);	\
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
    
    // The return values seems to be whether storage was enabled
    // previously, and not whether storage is enabled now.
    return retval;
  }
  
  
  
  int DisableGroupStorage (cGH* cgh, const char* groupname)
  {
    Checkpoint ("%*sDisableGroupStorage %s", 2*reflevel, "", groupname);
    
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
	for (int ti=0; ti<(int)scdata[group][var].size(); ++ti) {
	  switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)				\
	  case N:				\
	    delete (T*)scdata[group][var][ti];	\
	    break;
#include "typecase"
#undef TYPECASE
	  default:
	    abort();
	  }
	  scdata[group][var][ti] = 0;
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
	// save old values
	const vect<int,dim> oldrlb = rlb;
	const vect<int,dim> oldrub = rub;
	// calculate extent and centre
	const vect<int,dim> rextent = rub - rlb + rstr;
	const vect<int,dim> rcentre = rlb + (rextent / 2 / rstr) * rstr;
	// calculate new extent
	assert (all(rextent % hh->reffact == 0));
	const vect<int,dim> newrextent = rextent / hh->reffact;
	// refined boxes have smaller stride
	assert (all(rstr%hh->reffact == 0));
	rstr /= hh->reffact;
	// refine (arbitrarily) around the center only
	rlb = rcentre - (newrextent/2 / rstr) * rstr;
	rub = rlb + newrextent - rstr;
// 	// refine (arbitrarily) around the lower boundary only
// 	rlb = rlb;
// 	rub = rlb + newrextent - rstr;
	assert (all(rlb >= oldrlb && rub <= oldrub));
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
    // TODO: recompose arrays too
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
	    for (int ti=0; ti<CCTK_NumTimeLevelsFromVarI(n)-1; ++ti) {
	      // TODO: Which refinement level to use?
	      memcpy (scdata[group][var][ti], scdata[group][var][ti+1], sz);
	    }
	    break;
	  }
	  case CCTK_ARRAY: {
	    assert (group<(int)arrdata.size());
	    assert (var<(int)arrdata[group].data.size());
	    for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	      for (int ti=0; ti<CCTK_NumTimeLevelsFromVarI(n)-1; ++ti) {
		const int tmin = min(0, 2 - CCTK_NumTimeLevelsFromVarI(n));
		const int tl   = tmin + ti;
		arrdata[group].data[var]->copy(tl, reflevel, c, mglevel);
	      }
	    }
	    break;
	  }
	  case CCTK_GF: {
	    assert (group<(int)gfdata.size());
	    assert (var<(int)gfdata[group].data.size());
	    for (int c=0; c<hh->components(reflevel); ++c) {
	      for (int ti=0; ti<CCTK_NumTimeLevelsFromVarI(n)-1; ++ti) {
		const int tmin = min(0, 2 - CCTK_NumTimeLevelsFromVarI(n));
		const int tl   = tmin + ti;
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
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      // Restrict only groups with storage
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	
	const int tl = activetimelevel;
	
	switch (CCTK_GroupTypeI(group)) {
	  
	case CCTK_SCALAR:
	  break;
	  
	case CCTK_ARRAY:
	  for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	    for (int c=0; c<hh->components(reflevel); ++c) {
	      arrdata[group].data[var]->ref_restrict
		(tl, reflevel, c, mglevel);
	    }
	    for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	      arrdata[group].data[var]->sync (tl, reflevel, c, mglevel);
	    }
	  }
	  break;
	  
	case CCTK_GF:
	  for (int var=0; var<(int)gfdata[group].data.size(); ++var) {
	    for (int c=0; c<hh->components(reflevel); ++c) {
	      gfdata[group].data[var]->ref_restrict
		(tl, reflevel, c, mglevel);
	    }
	    for (int c=0; c<hh->components(reflevel); ++c) {
	      gfdata[group].data[var]->sync (tl, reflevel, c, mglevel);
	    }
	  }
	  break;
	  
	default:
	  abort();
	}
	
      }	// if group has storage
      
    } // loop over groups
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
  
  
  
  int MyProc (cGH* cgh)
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
      const int ti=0;
      assert (ti<(int)scdata[group][var].size());
      return scdata[group][var][ti] != 0;
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
  
  
  
  void set_reflevel (cGH* cgh, const int rl)
  {
    // Check
    assert (rl>=0 && rl<hh->reflevels());
    assert (component == -1);
    
    // Change
    reflevel = rl;
    const bbox<int,dim>& base = hh->baseextent;
    reflevelfact = (int)floor(pow(hh->reffact, reflevel)+0.5);
    cgh->cctk_delta_time = base_delta_time / reflevelfact;
    for (int d=0; d<dim; ++d) {
      cgh->cctk_gsh[d]
	= ((base.shape() / base.stride()
	    + dd->lghosts + dd->ughosts)[d] - 1) * reflevelfact + 1;
      cgh->cctk_levfac[d] = reflevelfact;
    }
  }
  
  
  
  void set_mglevel (cGH* cgh, const int ml)
  {
    assert (ml==0);
    assert (component==-1);
    mglevel = ml;
    cgh->cctk_convlevel = mglevel;
  }
  
  
  
  void set_component (cGH* cgh, const int c)
  {
    assert (c==-1 || (c>=0 && c<hh->components(reflevel)));
    component = c;
    
    if (component == -1) {
      // Global mode -- no component is active
      
      // Set Cactus parameters to pseudo values
      for (int d=0; d<dim; ++d) {
	cgh->cctk_lsh[d]      = 0xdeadbeef;
	cgh->cctk_bbox[2*d  ] = 0xdeadbeef;
	cgh->cctk_bbox[2*d+1] = 0xdeadbeef;
	cgh->cctk_lbnd[d]     = 0xdeadbeef;
	cgh->cctk_ubnd[d]     = 0xdeadbeef;
	for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	  cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = 0xdeadbeef;
	}
      }
      
      // Set local grid function and array sizes
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
      
      // Set Cactus pointers to data
      for (int n=0; n<CCTK_NumVars(); ++n) {
	for (int ti=0; ti<CCTK_NumTimeLevelsFromVarI(n); ++ti) {
	  
	  const int group = CCTK_GroupIndexFromVarI(n);
	  assert (group>=0);
	  const int var   = n - CCTK_FirstVarIndexI(group);
	  assert (var>=0);
	  
	  if (CCTK_QueryGroupStorageI(cgh, group)) {
	    // Group has storage
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR:
	      // Scalar variables can be accessed
	      assert (group<(int)scdata.size());
	      assert (var<(int)scdata[group].size());
	      assert (ti<(int)scdata[group][var].size());
	      cgh->data[n][ti] = scdata[group][var][ti];
	      break;
	      
	    case CCTK_ARRAY:
	    case CCTK_GF:
	      // Arrays and grid functions cannot be accessed
	      cgh->data[n][ti] = 0;
	      break;
	      
	    default:
	      abort();
	    }
	    
	  } else {
	    // Group has no storage
	    
	    cgh->data[n][ti] = 0;
	    
	  }
	  
	} // for ti
      }	// for n
      
    } else {
      // Local mode -- a component is active
      
      // Set Cactus parameters
      for (int d=0; d<dim; ++d) {
	const bbox<int,dim>& ext
	  = dd->boxes[reflevel][component][mglevel].exterior;
	cgh->cctk_lsh[d] = (ext.shape() / ext.stride())[d];
	cgh->cctk_lbnd[d] = (ext.lower() / ext.stride())[d];
	cgh->cctk_ubnd[d] = (ext.upper() / ext.stride())[d];
	assert (cgh->cctk_lsh[d]>=0 && cgh->cctk_lsh[d]<=cgh->cctk_gsh[d]);
	assert (cgh->cctk_lbnd[d]>=0 && cgh->cctk_ubnd[d]<cgh->cctk_gsh[d]);
	assert (cgh->cctk_lbnd[d]<=cgh->cctk_ubnd[d]+1);
	// No outer boundaries on the finer grids
	cgh->cctk_bbox[2*d  ]
	  = reflevel==0 && cgh->cctk_lbnd[d] == 0;
	cgh->cctk_bbox[2*d+1]
	  = reflevel==0 && cgh->cctk_ubnd[d] == cgh->cctk_gsh[d]-1;
#if 0
	// Do allow outer boundaries on the finer grids (but this is
	// generally inconsistent -- c. f. periodicity)
	const bbox<int,dim>& base = hh->baseextent;
	cgh->cctk_bbox[2*d  ] = (ext.lower() < base.lower())[d];
	cgh->cctk_bbox[2*d+1] = (ext.upper() > base.upper())[d];
#endif
	for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	  // TODO: support staggering
	  cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
	}
      }
      
      // Set local grid function and array sizes
      const bbox<int,dim> ext =
	dd->boxes[reflevel][component][mglevel].exterior;
      for (int d=0; d<dim; ++d) {
	gfsize[d] = ext.shape()[d] / ext.stride()[d];
      }
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	const int n0 = CCTK_FirstVarIndexI(group);
	switch (CCTK_GroupTypeFromVarI(n0)) {
	case CCTK_SCALAR:
	  break;
	case CCTK_ARRAY: {
	  const bbox<int,dim> ext =
	    arrdata[group].dd->
	    boxes[reflevel][component][mglevel].exterior;
	  for (int d=0; d<dim; ++d) {
	    arrdata[group].size[d] = ext.shape()[d] / ext.stride()[d];
	  }
	  break;
	}
	case CCTK_GF:
	  break;
	default:
	  abort();
	}
      }
      
      // Set Cactus pointers to data
      for (int n=0; n<CCTK_NumVars(); ++n) {
	for (int ti=0; ti<CCTK_NumTimeLevelsFromVarI(n); ++ti) {
	  
	  const int group = CCTK_GroupIndexFromVarI(n);
	  assert (group>=0);
	  const int var   = n - CCTK_FirstVarIndexI(group);
	  assert (var>=0);
	  const int tmin  = min(0, 2 - CCTK_NumTimeLevelsFromVarI(n));
	  const int tl    = tmin + ti;
	  
	  if (CCTK_QueryGroupStorageI(cgh, group)) {
	    // Group has storage
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR:
	      assert (group<(int)scdata.size());
	      assert (var<(int)scdata[group].size());
	      assert (ti<(int)scdata[group][var].size());
	      cgh->data[n][ti] = scdata[group][var][ti];
	      break;
	      
	    case CCTK_ARRAY:
	      assert (group<(int)arrdata.size());
	      assert (var<(int)arrdata[group].data.size());
	      cgh->data[n][ti]
		= ((*arrdata[group].data[var])
		   (tl, reflevel, component, mglevel)->storage());
	      break;
	      
	    case CCTK_GF:
	      assert (group<(int)gfdata.size());
	      assert (var<(int)gfdata[group].data.size());
	      cgh->data[n][ti]
		= ((*gfdata[group].data[var])
		   (tl, reflevel, component, mglevel)->storage());
	      break;
	      
	    default:
	      abort();
	    }
	    if (hh->is_local(reflevel,component)) {
	      assert (cgh->data[n][ti]);
	    } else {
	      assert (! cgh->data[n][ti]);
	    }
	    
	  } else {
	    // Group has no storage
	    
	    cgh->data[n][ti] = 0;
	    
	  } // if ! has storage
	  
	} // for ti
      } // for n
      
    }
  }
  
  
  
  MPI_Comm CarpetMPICommunicator ()
  {
    return dist::comm;
  }
  
} // namespace Carpet
