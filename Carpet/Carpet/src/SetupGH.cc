#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/dist.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/SetupGH.cc,v 1.40 2003/04/30 12:43:21 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_SetupGH_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  static bool CanTransferVariableType (cGH* cgh, const int group)
  {
    // Find out which types correspond to the default types
#if CCTK_INTEGER_PRECISION_1
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT1
#elif CCTK_INTEGER_PRECISION_2
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT2
#elif CCTK_INTEGER_PRECISION_4
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT4
#elif CCTK_INTEGER_PRECISION_8
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT8
#else
#  error "Unsupported default integer type"
#endif
    
#if CCTK_REAL_PRECISION_4
#  define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL4
#  define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX8
#elif CCTK_REAL_PRECISION_8
#  define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL8
#  define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX16
#elif CCTK_REAL_PRECISION_16
#  define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL16
#  define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX32
#else
#  error "Unsupported default real type"
#endif
    
    const int var0 = CCTK_FirstVarIndexI(group);
    const int type0 = CCTK_VarTypeI(var0);
    int type1;
    switch (type0) {
    case CCTK_VARIABLE_INT:
      type1 = CCTK_DEFAULT_INTEGER_TYPE;
      break;
    case CCTK_VARIABLE_REAL:
      type1 = CCTK_DEFAULT_REAL_TYPE;
      break;
    case CCTK_VARIABLE_COMPLEX:
      type1 = CCTK_DEFAULT_COMPLEX_TYPE;
      break;
    default:
      type1 = type0;
    }
    switch (type1) {
      
#ifdef CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      // This type is supported.
      return true;
#endif
      
#ifdef CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
#endif
#ifdef CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
#endif
#ifdef CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#endif
#ifdef CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#endif
#ifdef CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#endif
      // This type is not supported, but could be.
      return false;
      
    case CCTK_VARIABLE_BYTE:
#ifdef CCTK_INT1
    case CCTK_VARIABLE_INT1:
#endif
#ifdef CCTK_INT2
    case CCTK_VARIABLE_INT2:
#endif
#ifdef CCTK_INT4
    case CCTK_VARIABLE_INT4:
#endif
#ifdef CCTK_INT8
    case CCTK_VARIABLE_INT8:
#endif
      // This type is not supported, and cannot be.
      return false;
      
    default:
      assert (0);
    }
    
    // not reached
    return false;
  }
  
  
  
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    int ierr;
    
    assert (cgh->cctk_dim == dim);
    
    // Not sure what to do with that
    assert (convLevel==0);
    
    dist::pseudoinit();
    
    CCTK_VInfo (CCTK_THORNSTRING,
		"Carpet is running on %d processors", CCTK_nProcs(cgh));
    
    Waypoint ("starting SetupGH...");
    
    // Refinement information
    maxreflevels = max_refinement_levels;
    reffact = refinement_factor;
    maxreflevelfact = ipow(reffact, maxreflevels-1);
    
    // Multigrid information
    mglevels = multigrid_levels;
    mgfact = multigrid_factor;
    maxmglevelfact = ipow(mgfact, mglevels-1);
    
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
    const vect<int,dim> lb(0);
    const vect<int,dim> ub((npoints - 1) * str);
    
    const bbox<int,dim> baseext(lb, ub, str);
    
    // Allocate grid hierarchy
    hh = new gh<dim>(refinement_factor, vertex_centered,
		     multigrid_factor, vertex_centered,
		     baseext);
    
    // Allocate time hierarchy
    tt = new th<dim>(hh, 1.0);
    
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
    arrdata.resize(CCTK_NumGroups());
    
    // Allocate space for variables in group (but don't enable storage
    // yet)
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      cGroup gp;
      ierr = CCTK_GroupData (group, &gp);
      assert (!ierr);
      
      switch (gp.grouptype) {
	
      case CCTK_SCALAR:
      case CCTK_ARRAY: {
	
	int disttype;
	vect<int,dim> sizes(1), ghostsizes(0);
	
	switch (gp.grouptype) {
	  
	case CCTK_SCALAR:
	  // treat scalars as DIM=0, DISTRIB=const arrays
	  arrdata[group].info.dim = 0;
	  disttype = CCTK_DISTRIB_CONSTANT;
	  break;
	  
	case CCTK_ARRAY: {
	  assert (gp.dim>=1 || gp.dim<=dim);
	  arrdata[group].info.dim = gp.dim;
	  disttype = gp.disttype;
	  const CCTK_INT * const * const sz  = CCTK_GroupSizesI(group);
	  const CCTK_INT * const * const gsz = CCTK_GroupGhostsizesI(group);
	  for (int d=0; d<gp.dim; ++d) {
	    if (sz) sizes[d] = *sz[d];
	    if (gsz) ghostsizes[d] = *gsz[d];
	  }
	  break;
	}
	  
	default:
	  assert (0);
	}
	
	assert (disttype==CCTK_DISTRIB_CONSTANT
		|| disttype==CCTK_DISTRIB_DEFAULT);
	
	if (disttype==CCTK_DISTRIB_CONSTANT) {
	  sizes[dim-1] = (sizes[dim-1] - 2*ghostsizes[dim-1]) * CCTK_nProcs(cgh) + 2*ghostsizes[dim-1];
	}
	
	const vect<int,dim> alb(0);
	const vect<int,dim> aub(sizes-1);
	const vect<int,dim> astr(1);
	const bbox<int,dim> arrext(alb, aub, astr);
	
	arrdata[group].hh = new gh<dim>(refinement_factor, vertex_centered,
					multigrid_factor, vertex_centered,
					arrext);
	
	arrdata[group].tt = new th<dim>(arrdata[group].hh, 1.0);
	
	vect<int,dim> alghosts(0), aughosts(0);
	for (int d=0; d<gp.dim; ++d) {
	  alghosts[d] = ghostsizes[d];
	  aughosts[d] = ghostsizes[d];
	}
	
	arrdata[group].dd
	  = new dh<dim>(*arrdata[group].hh, alghosts, aughosts, 0);
	
	// Set refinement structure for scalars and arrays:
	
	// Set the basic extent
	vector<bbox<int,dim> >          bbs;
	vector<vect<vect<bool,2>,dim> > obs;
	bbs.push_back (arrdata[group].hh->baseextent);
	obs.push_back (vect<vect<bool,2>,dim>(vect<bool,2>(true)));
	
	// Split it into components, one for each processor
	SplitRegions_AlongZ (cgh, bbs, obs);
	
	// For all refinement levels (but there is only one)
	vector<vector<bbox<int,dim> > >          bbss(1);
	vector<vector<vect<vect<bool,2>,dim> > > obss(1);
	bbss[0] = bbs;
	obss[0] = obs;
	
	// For all multigrid levels
	gh<dim>::rexts bbsss;
	bbsss = hh->make_multigrid_boxes(bbss, mglevels);
	
	// Distribute onto processors
	vector<vector<int> > pss;
	MakeProcessors (cgh, bbsss, pss);
	
	// And recompose.  Done.
        char * groupname = CCTK_GroupName (group);
        assert (groupname);
        Checkpoint ("Recomposing grid array group %s", groupname);
        free (groupname);
	arrdata[group].hh->recompose (bbsss, obss, pss);
	
	break;
      }
	
      case CCTK_GF: {
	assert (gp.dim == dim);
	arrdata[group].info.dim = dim;
	arrdata[group].hh = hh;
	arrdata[group].tt = tt;
	arrdata[group].dd = dd;
	break;
      }
	
      default:
	assert (0);
      }
      
      arrdata[group].info.gsh         = new int [dim];
      arrdata[group].info.lsh         = new int [dim];
      arrdata[group].info.lbnd        = new int [dim];
      arrdata[group].info.ubnd        = new int [dim];
      arrdata[group].info.bbox        = new int [2*dim];
      arrdata[group].info.nghostzones = new int [dim];
      
      arrdata[group].data.resize(CCTK_NumVarsInGroupI(group));
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	arrdata[group].data[var] = 0;
      }
      
      arrdata[group].do_transfer = CanTransferVariableType (cgh, group);
      
    } // for group
    
    
    
    // Initialise cgh
    for (int d=0; d<dim; ++d) {
      cgh->cctk_nghostzones[d] = dd->lghosts[d];
    }
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      for (int d=0; d<dim; ++d) {
	((int*)arrdata[group].info.nghostzones)[d]
	  = arrdata[group].dd->lghosts[d];
      }
    }
    
    // Initialise current position
    reflevel  = 0;
    mglevel   = -1;
    component = -1;
    
    
    
    // Set initial refinement structure
    vector<bbox<int,dim> > bbs;
    vector<vect<vect<bool,2>,dim> > obs;
    if (strcmp(base_extents, "") == 0) {
      // default: one grid component covering everything
      bbs.push_back (hh->baseextent);
      obs.push_back (vect<vect<bool,2>,dim>(vect<bool,2>(true)));
    } else {
      // explicit grid components
      istringstream ext_str(base_extents);
      ext_str >> bbs;
      CCTK_VInfo (CCTK_THORNSTRING, "Using %d grid patches", bbs.size());
      cout << "grid-patches-are " << bbs << endl;
      if (bbs.size()<=0) {
	CCTK_WARN (0, "Cannot evolve with 0 grid patches");
      }
      istringstream ob_str (base_outerbounds);
      ob_str >> obs;
      cout << "outer-boundaries-are " << obs << endl;
      assert (obs.size() == bbs.size());
    }
    
    SplitRegions (cgh, bbs, obs);
    
    vector<vector<bbox<int,dim> > > bbss(1);
    vector<vector<vect<vect<bool,2>,dim> > > obss(1);
    bbss[0] = bbs;
    obss[0] = obs;
    
    gh<dim>::rexts bbsss;
    bbsss = hh->make_multigrid_boxes(bbss, mglevels);
    
    gh<dim>::rprocs pss;
    MakeProcessors (cgh, bbsss, pss);
    
    // Recompose grid hierarchy
    Checkpoint ("Recomposing grid functions");
    Recompose (cgh, bbsss, obss, pss);
    
    
    
    // Initialise time step on coarse grid
    base_delta_time = 1.0;
    
    // Current iteration
    iteration.resize(maxreflevels, 0);
    
    // Set current position (this time for real)
    set_reflevel  (cgh, 0);
    set_mglevel   (cgh, -1);
    set_component (cgh, -1);
    
    // Enable storage for all groups if desired
    // XXX
    if (true || enable_all_storage) {
      BEGIN_REFLEVEL_LOOP(cgh) {
	BEGIN_MGLEVEL_LOOP(cgh) {
	  for (int group=0; group<CCTK_NumGroups(); ++group) {
	    EnableGroupStorage (cgh, CCTK_GroupName(group));
	  }
	} END_MGLEVEL_LOOP(cgh);
      } END_REFLEVEL_LOOP(cgh);
    }
    
    Waypoint ("done with SetupGH.");
    
    // We register only once, ergo we get only one handle, ergo there
    // is only one grid hierarchy for us.  We store that statically,
    // so there is no need to pass anything to Cactus.
    return 0;
  }
  
} // namespace Carpet
