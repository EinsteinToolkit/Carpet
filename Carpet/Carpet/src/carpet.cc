// It is assumed that the number of components of all arrays is equal
// to the number of components of the grid functions, and that their
// distribution onto the processors is the same, and that all
// processors own the same number of components.

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <complex>

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

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Attic/carpet.cc,v 1.27 2001/06/12 15:14:16 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  static void Recompose (cGH* cgh);
  static void CycleTimeLevels (cGH* cgh);
  static void Restrict (cGH* cgh);
  
  enum checktimes { currenttime,
		    currenttimebutnotifonly,
		    allbutlasttime,
		    allbutcurrenttime,
		    alltimes };
  
  static void Poison (cGH* cgh, checktimes where);
  static void PoisonGroup (cGH* cgh, int group, checktimes where);
  static void PoisonCheck (cGH* cgh, checktimes where);
  
  static void CalculateChecksums (cGH* cgh, checktimes where);
  static void CheckChecksums (cGH* cgh, checktimes where);
  
  static int mintl (checktimes where, int num_tl);
  static int maxtl (checktimes where, int num_tl);
  
  // Debugging output
  static void Checkpoint (const char* fmt, ...);
  
  // Error output
  static void UnsupportedVarType (int vindex);
  
  
  
  // Handle from CCTK_RegisterGHExtension
  int GHExtension;
  
  // Maximum number of refinement levels
  int maxreflevels;
  
  // Refinement factor on finest grid
  int maxreflevelfact;
  
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
  vector<scdesc> scdata;	// [group]
  
  // Data for arrays
  // TODO: have replicated arrays
  vector<arrdesc> arrdata;	// [group]
  
  // The grid hierarchy
  gh<dim>* hh;
  th* tt;
  dh<dim>* dd;
  int gfsize[dim];
  
  // Data for grid functions
  vector<gfdesc> gfdata;	// [group]
  
  // Checksums
  vector<vector<vector<vector<ckdesc> > > > checksums; // [n][rl][tl][c]
  
  
  
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
    
    assert (cgh->cctk_dim == dim);
    
    // Not sure what to do with that
    assert (convLevel==0);
    
    dist::pseudoinit();
    Checkpoint ("starting SetupGH...");
    
    // Refinement information
    maxreflevels = max_refinement_levels;
    maxreflevelfact
      = (int)floor(pow((double)refinement_factor, maxreflevels-1) + 0.5);
    
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
    tt = new th(hh, maxreflevelfact);
    
    // Allocate data hierarchy
    dd = new dh<dim>(*hh, lghosts, ughosts, prolongation_order_space);
    
    if (max_refinement_levels > 1) {
      const int prolongation_stencil_size = dd->prolongation_stencil_size();
      const int min_nghosts
	= ((prolongation_stencil_size + refinement_factor - 1)
	   / (refinement_factor-1));
      if (any(min(lghosts,ughosts) < min_nghosts)) {
	CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		    "There are not enough ghost zones for the desired spatial prolongation order.  With Carpet::prolongation_order_space=%d, you need at least %d ghost zones.",
		    prolongation_order_space, min_nghosts);
      }
    }
    
    // Allocate space for groups
    scdata.resize(CCTK_NumGroups());
    arrdata.resize(CCTK_NumGroups());
    gfdata.resize(CCTK_NumGroups());
    
    // Allocate space for variables in group (but don't enable storage
    // yet)
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      switch (CCTK_GroupTypeI(group)) {
	
      case CCTK_SCALAR: {
	scdata[group].data.resize(CCTK_NumVarsInGroupI(group));
	for (int var=0; var<(int)scdata[group].data.size(); ++var) {
	  const int n = (CCTK_FirstVarIndexI(group) + var);
	  scdata[group].data[var].resize(maxreflevels);
	  for (int rl=0; rl<maxreflevels; ++rl) {
	    scdata[group].data[var][rl].resize(CCTK_NumTimeLevelsFromVarI(n));
	    for (int tl=0; tl<(int)scdata[group].data[var].size(); ++tl) {
	      scdata[group].data[var][rl][tl] = 0;
	    }
	  }
	}
	break;
      }
	
      case CCTK_ARRAY: {
	cGroup gp;
	CCTK_GroupData (group, &gp);
	
	// TODO
	assert (CCTK_GroupDimI(group) == dim);
	
	// TODO
// 	assert (gp.disttype == CCTK_DISTRIB_CONSTANT);
	assert (gp.disttype == CCTK_DISTRIB_DEFAULT);
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
	
	arrdata[group].tt = new th (arrdata[group].hh, maxreflevelfact);
	
	vect<int,dim> alghosts, aughosts;
	for (int d=0; d<dim; ++d) {
	  alghosts[d] = (*CCTK_GroupGhostsizesI(group))[d];
	  aughosts[d] = (*CCTK_GroupGhostsizesI(group))[d];
	}
	
	arrdata[group].dd = new dh<dim>(*arrdata[group].hh, alghosts, aughosts,
					prolongation_order_space);
	
	if (max_refinement_levels > 1) {
	  const int prolongation_stencil_size
	    = arrdata[group].dd->prolongation_stencil_size();
	  const int min_nghosts
	    = ((prolongation_stencil_size + refinement_factor - 2)
	       / (refinement_factor-1));
	  if (any(min(alghosts,aughosts) < min_nghosts)) {
	    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
			"There are not enough ghost zones for the desired spatial prolongation order for the grid function group \"%s\".  With Carpet::prolongation_order_space=%d, you need at least %d ghost zones.",
			CCTK_GroupName(group),
			prolongation_order_space, min_nghosts);
	  }
	}
	
	arrdata[group].data.resize(CCTK_NumVarsInGroupI(group));
	for (int var=0; var<(int)scdata[group].data.size(); ++var) {
	  arrdata[group].data[var] = 0;
	}
	break;
      }
	
      case CCTK_GF: {
	/* TODO */
	assert (CCTK_GroupDimI(group) == dim);
	
	gfdata[group].data.resize(CCTK_NumVarsInGroupI(group));
	for (int var=0; var<(int)scdata[group].data.size(); ++var) {
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
    
    // Current iteration
    iteration.resize(maxreflevels, 0);
    
    // Set current position (this time for real)
    set_reflevel  (cgh, 0);
    set_mglevel   (cgh, 0);
    set_component (cgh, -1);
    
    // Enable storage for all groups if desired
    // XXX
    if (true || enable_all_storage) {
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	EnableGroupStorage (cgh, CCTK_GroupName(group));
      }
    }
    
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
    Checkpoint ("PARAMCHECK");
    CCTK_ScheduleTraverse ("CCTK_PARAMCHECK", cgh, CallFunction);
    CCTKi_FinaliseParamWarn();
    
    BEGIN_REFLEVEL_LOOP(cgh) {
      
      // Checking
      Poison (cgh, allbutlasttime);
      
      // Set up the grid
      Checkpoint ("%*sScheduling BASEGRID", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_BASEGRID", cgh, CallFunction);
      if (reflevel==0) {
	base_delta_time = cgh->cctk_delta_time;
      } else {
// 	assert (abs(cgh->cctk_delta_time - base_delta_time / reflevelfactor)
// 		< 1e-6 * base_delta_time);
	// This circumvents a bug in CactusBase/Time
	cgh->cctk_delta_time = base_delta_time / reflevelfact;
      }
      
      // Set up the initial data
      Checkpoint ("%*sScheduling INITIAL", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_INITIAL", cgh, CallFunction);
      Checkpoint ("%*sScheduling POSTINITIAL", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_POSTINITIAL", cgh, CallFunction);
      
      // Recover
      Checkpoint ("%*sScheduling RECOVER_VARIABLES", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_RECOVER_VARIABLES", cgh, CallFunction);
      Checkpoint ("%*sScheduling CPINITIAL", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_CPINITIAL", cgh, CallFunction);
      
      // Poststep
      Checkpoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
      
      // Checking
      PoisonCheck (cgh, allbutlasttime);
      
    } END_REFLEVEL_LOOP(cgh);
    
    BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
      
      // Restrict
      Restrict (cgh);
      
      // Checking
      CalculateChecksums (cgh, allbutlasttime);
      
      // Analysis
      Checkpoint ("%*sScheduling ANALYSIS", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
      
      // Output
      Checkpoint ("%*sOutputGH", 2*reflevel, "");
      CCTK_OutputGH (cgh);
      
      // Checking
      CheckChecksums (cgh, allbutlasttime);
      
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
      
      // Advance time
      ++cgh->cctk_iteration;
      cgh->cctk_time += base_delta_time / maxreflevelfact;
      
      Checkpoint ("Evolving iteration %d...", cgh->cctk_iteration);
      
      BEGIN_REFLEVEL_LOOP(cgh) {
	if ((cgh->cctk_iteration-1) % (maxreflevelfact/reflevelfact) == 0) {
	  
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
	  
	  // Checking
	  CalculateChecksums (cgh, allbutcurrenttime);
	  Poison (cgh, currenttimebutnotifonly);
	  
	  // Evolve
	  Checkpoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	  Checkpoint ("%*sScheduling EVOL", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	  Checkpoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	  
	  // Checking
	  PoisonCheck (cgh, currenttimebutnotifonly);
	  
	}
      } END_REFLEVEL_LOOP(cgh);
      
      BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
	if (cgh->cctk_iteration % (maxreflevelfact/reflevelfact) == 0) {
	  
	  // Restrict
	  Restrict (cgh);
	  
	  // Checking
	  CalculateChecksums (cgh, currenttime);
	  
	  // Checkpoint
	  Checkpoint ("%*sScheduling CHECKPOINT", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_CHECKPOINT", cgh, CallFunction);
	  
	  // Analysis
	  Checkpoint ("%*sScheduling ANALYSIS", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
	  
	  // Output
	  Checkpoint ("%*sOutputGH", 2*reflevel, "");
	  CCTK_OutputGH (cgh);
	  
	  // Checking
	  CheckChecksums (cgh, alltimes);
	  
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
      Checkpoint ("%*sScheduling TERMINATE", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_TERMINATE", cgh, CallFunction);
    } END_REFLEVEL_LOOP(cgh);
    
    // Shutdown
    BEGIN_REFLEVEL_LOOP(cgh) {
      Checkpoint ("%*sScheduling SHUTDOWN", 2*reflevel, "");
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
    
//     Checkpoint ("%*sStarting CallFunction...", 2*reflevel, "");
    
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
    
    // The return value indicates whether the grid functions have been
    // synchronised.
    // return 0: let the flesh do the synchronisation, if necessary
    return 0;
  }
  
  
  
  int SyncGroup (cGH* cgh, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (hh->components(reflevel) > 1) assert (component == -1);
    
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
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tl = 0;
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      // TODO: Check whether the local values are consistent
      break;
      
    case CCTK_ARRAY:
      assert (group<(int)arrdata.size());
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	if (num_tl>1 && reflevel>0) {
	  for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	    arrdata[group].data[var]->ref_bnd_prolongate
	      (tl, reflevel, c, mglevel,
	       prolongation_order_space, prolongation_order_time);
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
	      (tl, reflevel, c, mglevel,
	       prolongation_order_space, prolongation_order_time);
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
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      // storage was enabled previously
      return 1;
    }
    
    // There is a difference between the Cactus time levels and the
    // Carpet time levels.  If there are n time levels, then the
    // Cactus time levels are numbered 0 ... n-1, with the current
    // time level being n-1.  In Carpet, the time levels are numbered
    // -(n-1) ... 0, where the current time level is always 0.
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tmin = 1 - num_tl;
    const int tmax = 0;
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      assert (scdata[group].data.size()==0
	      || scdata[group].data[0].size()==0
	      || scdata[group].data[0][0].size()==0
	      || scdata[group].data[0][0][0] == 0);
      for (int var=0; var<(int)scdata[group].data.size(); ++var) {
	for (int rl=0; rl<(int)scdata[group].data[var].size(); ++rl) {
	  for (int tl=0; tl<(int)scdata[group].data[var][rl].size(); ++tl) {
	    const int n = n0 + var;
	    switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	    case N:					\
	      scdata[group].data[var][rl][tl] = new T;	\
	      break;
#include "typecase"
#undef TYPECASE
	    default:
	      UnsupportedVarType(n);
	    }
	  }
	}
      }
      break;
      
    case CCTK_ARRAY:
      assert (arrdata[group].data.size()==0
	      || arrdata[group].data[0] == 0);
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
	  UnsupportedVarType(n);
	}
      }
      break;
      
    case CCTK_GF:
      assert (gfdata[group].data.size()==0
	      || gfdata[group].data[0] == 0);
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
	  UnsupportedVarType(n);
	}
      }
      break;
      
    default:
      abort();
    }
    
    // Reinitialise Cactus variables
    set_component (cgh, component);
    PoisonGroup (cgh, group, alltimes);
    
    // storage was not enabled previously
    return 0;
  }
  
  
  
  int DisableGroupStorage (cGH* cgh, const char* groupname)
  {
    Checkpoint ("%*sDisableGroupStorage %s", 2*reflevel, "", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      // storage was disabled previously
      return 0;
    }
    
    // XXX
    CCTK_WARN (1, "Cannot disable storage -- storage management is not yet consistent for FMR");
    return 1;
    
    const int n0 = CCTK_FirstVarIndexI(group);
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      if (scdata[group].data.size()==0
	  || scdata[group].data[0].size()==0
	  || scdata[group].data[0][0].size()==0
	  || scdata[group].data[0][0][0] == 0) {
	// group already has no storage
	break;
      }
      for (int var=0; var<(int)scdata[group].data.size(); ++var) {
	const int n = n0 + var;
	for (int rl=0; rl<(int)scdata[group].data[var].size(); ++rl) {
	  for (int tl=0; tl<(int)scdata[group].data[var][rl].size(); ++tl) {
	    switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)						\
	    case N:						\
	      delete (T*)scdata[group].data[var][rl][tl];	\
	      break;
#include "typecase"
#undef TYPECASE
	    default:
	      UnsupportedVarType(n);
	    }
	    scdata[group].data[var][rl][tl] = 0;
	  }
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
	  UnsupportedVarType(n);
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
	  UnsupportedVarType(n);
	}
	gfdata[group].data[var] = 0;
      }
      break;
      
    default:
      abort();
    }
    
    // Reinitialise Cactus variables
    set_component (cgh, component);
    
    // storage was not disabled previously
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
	    assert (var<(int)scdata[group].data.size());
	    const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	    assert (num_tl>0);
	    void* tmpdata = scdata[group].data[var][reflevel][0];
	    for (int tl=1; tl<num_tl; ++tl) {
	      scdata[group].data[var][reflevel][tl]
		= scdata[group].data[var][reflevel][tl-1];
	    }
	    scdata[group].data[var][reflevel][0] = tmpdata;
	    tmpdata = 0;
	    break;
	  }
	  case CCTK_ARRAY: {
	    assert (group<(int)arrdata.size());
	    assert (var<(int)arrdata[group].data.size());
	    for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	      arrdata[group].data[var]->cycle (reflevel, c, mglevel);
	    }
	    break;
	  }
	  case CCTK_GF: {
	    assert (group<(int)gfdata.size());
	    assert (var<(int)gfdata[group].data.size());
	    for (int c=0; c<hh->components(reflevel); ++c) {
	      gfdata[group].data[var]->cycle (reflevel, c, mglevel);
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
	
	const int tl = 0;
	
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
      assert (var<(int)scdata[group].data.size());
      assert (reflevel<(int)scdata[group].data[var].size());
      const int tl=0;
      assert (tl<(int)scdata[group].data[var][reflevel].size());
      return scdata[group].data[var][reflevel][tl] != 0;
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
  
  
  
  void Poison (cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! poison_new_timelevels) return;
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	PoisonGroup (cgh, group, where);
      } // if has storage
    } // for group
  }
  
  
  
  void PoisonGroup (cGH* cgh, const int group, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (! poison_new_timelevels) return;
    
    Checkpoint ("%*sPoisonGroup %s", 2*reflevel, "", CCTK_GroupName(group));
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot poison group \"%s\" because it has no storage",
		  CCTK_GroupName(group));
      return;
    }
    
    for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
      
      const int n = CCTK_FirstVarIndexI(group) + var;
      const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n));
      assert (sz>0);
      
      const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
      assert (num_tl>0);
      const int min_tl = mintl(where, num_tl);
      const int max_tl = maxtl(where, num_tl);
      
      for (int tl=min_tl; tl<=max_tl; ++tl) {
	
	switch (CCTK_GroupTypeFromVarI(n)) {
	case CCTK_SCALAR: {
	  assert (group<(int)scdata.size());
	  assert (var<(int)scdata[group].data.size());
	  memset (cgh->data[n][tl], poison_value, sz);
	  break;
	}
	case CCTK_ARRAY:
	case CCTK_GF: {
	  BEGIN_COMPONENT_LOOP(cgh) {
	    if (hh->is_local(reflevel,component)) {
	      int np = 1;
	      for (int d=0; d<dim; ++d) {
		np *= *CCTK_ArrayGroupSizeI(cgh, d, group);
	      }
	      memset (cgh->data[n][tl], poison_value, np*sz);
	    }
	  } END_COMPONENT_LOOP(cgh);
	  break;
	}
	default:
	  abort();
	}
	
      } // for tl
      
    } // for var
  }
  
  
  
  void PoisonCheck (cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! check_for_poison) return;
    
    Checkpoint ("%*sPoisonCheck", 2*reflevel, "");
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  
	  const int n = CCTK_FirstVarIndexI(group) + var;
	  
	  const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	  assert (num_tl>0);
	  const int min_tl = mintl(where, num_tl);
	  const int max_tl = maxtl(where, num_tl);
	  
	  for (int tl=min_tl; tl<=max_tl; ++tl) {
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR: {
	      bool poisoned=false;
	      switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)						\
	      case N: {						\
		T worm;						\
		memset (&worm, poison_value, sizeof(worm));	\
		poisoned = *(const T*)cgh->data[n][tl] == worm;	\
		break;						\
	      }
#include "typecase"
#undef TYPECASE
	      default:
		UnsupportedVarType(n);
	      }
	      if (poisoned) {
		char* fullname = CCTK_FullName(n);
		CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
			    "The variable \"%s\" contains poison in timelevel %d",
			    fullname, tl);
		free (fullname);
	      }
	      break;
	    }
	      
	    case CCTK_ARRAY:
	    case CCTK_GF: {
	      BEGIN_COMPONENT_LOOP(cgh) {
		if (hh->is_local(reflevel,component)) {
		  int size[dim];
		  for (int d=0; d<dim; ++d) {
		    size[d] = *CCTK_ArrayGroupSizeI(cgh, d, group);
		  }
		  const int tp = CCTK_VarTypeI(n);
		  const void* const data = cgh->data[n][tl];
		  int numpoison=0;
		  for (int k=0; k<size[2]; ++k) {
		    for (int j=0; j<size[1]; ++j) {
		      for (int i=0; i<size[0]; ++i) {
			const int idx = CCTK_GFINDEX3D(cgh,i,j,k);
			bool poisoned=false;
			switch (tp) {
#define TYPECASE(N,T)							\
			case N: {					\
			  T worm;					\
			  memset (&worm, poison_value, sizeof(worm));	\
			  poisoned = ((const T*)data)[idx] == worm;	\
			  break;					\
			}
#include "typecase"
#undef TYPECASE
			default:
			  UnsupportedVarType(n);
			}
			if (poisoned) {
			  ++numpoison;
			  if (numpoison<=10) {
			    char* fullname = CCTK_FullName(n);
			    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
					"The variable \"%s\" contains poison at [%d,%d,%d] in timelevel %d",
					fullname, i,j,k, tl);
			    free (fullname);
			  }
			} // if poisoned
		      } // for i
		    } // for j
		  } // for k
		  if (numpoison>10) {
		    char* fullname = CCTK_FullName(n);
		    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
				"The variable \"%s\" contains poison at %d locations in timelevel %d; not all locations were printed.",
				fullname, numpoison, tl);
		    free (fullname);
		  }
		} // if is local
	      } END_COMPONENT_LOOP(cgh);
	      break;
	    }
	      
	    default:
	      abort();
	    }
	    
	  } // for tl
	  
	} // for var
      } // if has storage
    } // for group
  }
  
  
  
  void CalculateChecksums (cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("%*sCalculateChecksums", 2*reflevel, "");
    
    checksums.resize(CCTK_NumVars());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	
	const int n = CCTK_FirstVarIndexI(group) + var;
	const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	assert (num_tl>0);
	const int min_tl = mintl(where, num_tl);
	const int max_tl = maxtl(where, num_tl);
	
	checksums[n].resize(maxreflevels);
	checksums[n][reflevel].resize(num_tl);
	for (int tl=min_tl; tl<=max_tl; ++tl) {
	  switch (CCTK_GroupTypeFromVarI(n)) {
	    
	  case CCTK_SCALAR: {
	    checksums[n][reflevel][tl].resize(1);
	    const int c=0;
	    checksums[n][reflevel][tl][c].valid = false;
	    break;
	  }
	      
	  case CCTK_ARRAY:
	  case CCTK_GF: {
	    checksums[n][reflevel][tl].resize(hh->components(reflevel));
	    BEGIN_COMPONENT_LOOP(cgh) {
	      checksums[n][reflevel][tl][component].valid = false;
	    } END_COMPONENT_LOOP(cgh);
	    break;
	  }
	    
	  default:
	    abort();
	  }
	}
      }
    }
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  
	  const int n = CCTK_FirstVarIndexI(group) + var;
	  const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n));
	  assert (sz>0);
	  const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	  assert (num_tl>0);
	  const int min_tl = mintl(where, num_tl);
	  const int max_tl = maxtl(where, num_tl);
	  
	  for (int tl=min_tl; tl<=max_tl; ++tl) {
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR: {
	      const int c=0;
	      int chk = 0;
	      const void* data = cgh->data[n][tl];
	      for (int i=0; i<sz/(int)sizeof(chk); ++i) {
		chk += ((const int*)data)[i];
	      }
	      checksums[n][reflevel][tl][c].sum = chk;
	      checksums[n][reflevel][tl][c].valid = true;
	      break;
	    }
	      
	    case CCTK_ARRAY:
	    case CCTK_GF: {
	      BEGIN_COMPONENT_LOOP(cgh) {
		if (hh->is_local(reflevel,component)) {
		  int np = 1;
		  for (int d=0; d<dim; ++d) {
		    np *= *CCTK_ArrayGroupSizeI(cgh, d, group);
		  }
		  const void* data = cgh->data[n][tl];
		  int chk = 0;
		  for (int i=0; i<np*sz/(int)sizeof(chk); ++i) {
		    chk += ((const int*)data)[i];
		  }
		  checksums[n][reflevel][tl][component].sum = chk;
		  checksums[n][reflevel][tl][component].valid = true;
		}
	      } END_COMPONENT_LOOP(cgh);
	      break;
	    }
	      
	    default:
	      abort();
	    }
	    
	  } // for tl
	} // for var
      }	// if has storage
    } // for group
  }
  
  
  
  void CheckChecksums (cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("%*sCheckChecksums", 2*reflevel, "");
    
    assert ((int)checksums.size()==CCTK_NumVars());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  
	  const int n = CCTK_FirstVarIndexI(group) + var;
	  const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n));
	  assert (sz>0);
	  const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	  assert (num_tl>0);
	  const int min_tl = mintl(where, num_tl);
	  const int max_tl = maxtl(where, num_tl);
	  
	  assert ((int)checksums[n].size()==maxreflevels);
	  assert ((int)checksums[n][reflevel].size()==num_tl);
	  
	  for (int tl=min_tl; tl<=max_tl; ++tl) {
	    
	    bool unexpected_change = false;
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR: {
	      assert ((int)checksums[n][reflevel][tl].size()==1);
	      const int c=0;
	      if (checksums[n][reflevel][tl][c].valid) {
		int chk = 0;
		const void* data = cgh->data[n][tl];
		for (int i=0; i<sz/(int)sizeof(chk); ++i) {
		  chk += ((const int*)data)[i];
		}
		unexpected_change
		  = (unexpected_change
		     || chk != checksums[n][reflevel][tl][c].sum);
	      }
	      break;
	    }
	      
	    case CCTK_ARRAY:
	    case CCTK_GF: {
	      assert ((int)checksums[n][reflevel][tl].size()
		      == hh->components(reflevel));
	      BEGIN_COMPONENT_LOOP(cgh) {
		if (checksums[n][reflevel][tl][component].valid) {
		  if (hh->is_local(reflevel,component)) {
		    int np = 1;
		    for (int d=0; d<dim; ++d) {
		      np *= *CCTK_ArrayGroupSizeI(cgh, d, group);
		    }
		    const void* data = cgh->data[n][tl];
		    int chk = 0;
		    for (int i=0; i<np*sz/(int)sizeof(chk); ++i) {
		      chk += ((const int*)data)[i];
		    }
		    unexpected_change
		      = (unexpected_change
			 || chk != checksums[n][reflevel][tl][component].sum);
		  }
		}
	      } END_COMPONENT_LOOP(cgh);
	      break;
	    }
	      
	    default:
	      abort();
	    }
	    
	    if (unexpected_change) {
	      char* fullname = CCTK_FullName(n);
	      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
			  "Timelevel %d of the variable \"%s\" has changed unexpectedly.",
			  tl, fullname);
	      free (fullname);
	    }
	    
	  } // for tl
	  
	} // for var
      }	// if has storage
    } // for group
  }
  
  
  
  int mintl (const checktimes where, const int num_tl)
  {
    assert (num_tl>0);
    switch (where) {
    case currenttime:
      return 0;
    case currenttimebutnotifonly:
      // don't include current time if there is only one time level
      return num_tl>1 ? 0: 1;
    case allbutlasttime:
      // do include current time if there is only one time level
      return 0;
    case allbutcurrenttime:
      return 1;
    case alltimes:
      return 0;
    default:
      abort();
    }
  }
  
  int maxtl (const checktimes where, const int num_tl)
  {
    assert (num_tl>0);
    switch (where) {
    case currenttime:
      return 0;
    case currenttimebutnotifonly:
      return 0;
    case allbutlasttime:
      return num_tl-2;
    case allbutcurrenttime:
      return num_tl-1;
    case alltimes:
      return num_tl-1;
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
  
  
  
  void UnsupportedVarType (const int vindex)
  {
    assert (vindex>=0 && vindex<CCTK_NumVars());
    CCTK_VWarn
      (0, __LINE__, __FILE__, CCTK_THORNSTRING,
       "Carpet does not support the type of the variable \"%s\".\n"
       "Either enable support for this type, "
       "or change the type of this variable.", CCTK_FullName(vindex));
  }
  
  
  
  void set_reflevel (cGH* cgh, const int rl)
  {
    // Check
    assert (rl>=0 && rl<hh->reflevels());
    assert (component == -1);
    
    // Change
    reflevel = rl;
    const bbox<int,dim>& base = hh->baseextent;
    reflevelfact = (int)floor(pow((double)hh->reffact, reflevel)+0.5);
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
	for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n); ++tl) {
	  
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
	      assert (var<(int)scdata[group].data.size());
	      assert (reflevel<(int)scdata[group].data[var].size());
	      assert (tl<(int)scdata[group].data[var][reflevel].size());
	      cgh->data[n][tl] = scdata[group].data[var][reflevel][tl];
	      break;
	      
	    case CCTK_ARRAY:
	    case CCTK_GF:
	      // Arrays and grid functions cannot be accessed
	      cgh->data[n][tl] = 0;
	      break;
	      
	    default:
	      abort();
	    }
	    
	  } else {
	    // Group has no storage
	    
	    cgh->data[n][tl] = 0;
	    
	  }
	  
	} // for tl
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
	
	const int group = CCTK_GroupIndexFromVarI(n);
	assert (group>=0);
	const int var   = n - CCTK_FirstVarIndexI(group);
	assert (var>=0);
	const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	assert (num_tl>0);
	
	for (int tl=0; tl<num_tl; ++tl) {
	  
	  if (CCTK_QueryGroupStorageI(cgh, group)) {
	    // Group has storage
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR:
	      assert (group<(int)scdata.size());
	      assert (var<(int)scdata[group].data.size());
	      assert (reflevel<(int)scdata[group].data[var].size());
	      assert (tl<(int)scdata[group].data[var][reflevel].size());
	      cgh->data[n][tl] = scdata[group].data[var][reflevel][tl];
	      break;
	      
	    case CCTK_ARRAY:
	      assert (group<(int)arrdata.size());
	      assert (var<(int)arrdata[group].data.size());
	      cgh->data[n][tl]
		= ((*arrdata[group].data[var])
		   (-tl, reflevel, component, mglevel)->storage());
	      break;
	      
	    case CCTK_GF:
	      assert (group<(int)gfdata.size());
	      assert (var<(int)gfdata[group].data.size());
	      cgh->data[n][tl]
		= ((*gfdata[group].data[var])
		   (-tl, reflevel, component, mglevel)->storage());
	      break;
	      
	    default:
	      abort();
	    }
	    if (CCTK_GroupTypeFromVarI(n)==CCTK_SCALAR
		|| hh->is_local(reflevel,component)) {
	      assert (cgh->data[n][tl]);
	    } else {
	      assert (! cgh->data[n][tl]);
	    }
	    
	  } else {
	    // Group has no storage
	    
	    cgh->data[n][tl] = 0;
	    
	  } // if ! has storage
	  
	} // for tl
      } // for n
      
    }
  }
  
  
  
  MPI_Comm CarpetMPICommunicator ()
  {
    return dist::comm;
  }
  
} // namespace Carpet
