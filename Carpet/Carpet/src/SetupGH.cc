#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/dist.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/SetupGH.cc,v 1.7 2001/08/26 14:43:25 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
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
    
    // Allocate grid hierarchy for scalars
    const vect<int,dim> lb0(0);
    const vect<int,dim> ub0(0);
    const bbox<int,dim> baseext0(lb0, ub0, str);
    hh0 = new gh<dim>(refinement_factor, vertex_centered,
		      multigrid_factor, vertex_centered,
		      baseext0);
    
    // Allocate time hierarchy for scalars
    tt0 = new th(hh0, maxreflevelfact);
    
    // Allocate data hierarchy for scalars
    const vect<int,dim> lghosts0(0);
    const vect<int,dim> ughosts0(0);
    dd0 = new dh<dim>(*hh0, lghosts0, ughosts0, prolongation_order_space);
    
    // Allocate space for groups
    arrdata.resize(CCTK_NumGroups());
    
    // Allocate space for variables in group (but don't enable storage
    // yet)
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      switch (CCTK_GroupTypeI(group)) {
	
      case CCTK_SCALAR: {
	arrdata[group].info.dim = 0;
	arrdata[group].hh = hh0;
	arrdata[group].tt = tt0;
	arrdata[group].dd = dd0;
	break;
      }
	
      case CCTK_ARRAY: {
	cGroup gp;
	CCTK_GroupData (group, &gp);
	
	assert (gp.dim>=1 || gp.dim<=dim);
	arrdata[group].info.dim = gp.dim;
	
	switch (gp.disttype) {
	case CCTK_DISTRIB_CONSTANT:
	case CCTK_DISTRIB_DEFAULT:
	  break;
	default:
	  abort();
	}
	
	const CCTK_INT * const * const sz  = CCTK_GroupSizesI(group);
	const CCTK_INT * const * const gsz = CCTK_GroupGhostsizesI(group);
	vect<int,dim> sizes(1), ghostsizes(0);
	for (int d=0; d<gp.dim; ++d) {
	  if (sz) sizes[d] = *sz[d];
	  if (gsz) ghostsizes[d] = *gsz[d];
	}
	
	vect<int,dim> alb(0), aub(stride), astr(stride);
	for (int d=0; d<gp.dim; ++d) {
	  if (gp.disttype==CCTK_DISTRIB_CONSTANT && d==gp.dim-1) {
	    aub[d] = astr[d] * ((CCTK_nProcs(cgh) * sizes[d]
				 - (CCTK_nProcs(cgh)-1) * ghostsizes[d]) - 1);
	  } else {
	    aub[d] = astr[d] * (sizes[d]-1);
	  }
	}
	const bbox<int,dim> arrext(alb, aub, astr);
	
	arrdata[group].hh = new gh<dim>(refinement_factor, vertex_centered,
					multigrid_factor, vertex_centered,
					arrext);
	
	arrdata[group].tt = new th(arrdata[group].hh, maxreflevelfact);
	
	vect<int,dim> alghosts(0), aughosts(0);
	for (int d=0; d<gp.dim; ++d) {
	  alghosts[d] = ghostsizes[d];
	  aughosts[d] = ghostsizes[d];
	}
	
	arrdata[group].dd
	  = new dh<dim>(*arrdata[group].hh, alghosts, aughosts,
			prolongation_order_space);
	
	if (gp.disttype == CCTK_DISTRIB_DEFAULT) {
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
	}
	break;
      }
	
      case CCTK_GF: {
	assert (CCTK_GroupDimI(group) == dim);
	arrdata[group].info.dim = dim;
	arrdata[group].hh = hh;
	arrdata[group].tt = tt;
	arrdata[group].dd = dd;
	break;
      }
	
      default:
	abort();
      }
      
      arrdata[group].info.gsh         = (int*)malloc(  dim * sizeof(int));
      arrdata[group].info.lsh         = (int*)malloc(  dim * sizeof(int));
      arrdata[group].info.lbnd        = (int*)malloc(  dim * sizeof(int));
      arrdata[group].info.ubnd        = (int*)malloc(  dim * sizeof(int));
      arrdata[group].info.bbox        = (int*)malloc(2*dim * sizeof(int));
      arrdata[group].info.nghostzones = (int*)malloc(  dim * sizeof(int));
      
      arrdata[group].data.resize(CCTK_NumVarsInGroupI(group));
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	arrdata[group].data[var] = 0;
      }
    }
    
    // Initialise cgh
    for (int d=0; d<dim; ++d) {
      cgh->cctk_nghostzones[d] = dd->lghosts[d];
    }
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      for (int d=0; d<dim; ++d) {
	cgh->cctk_nghostzones[d] = arrdata[group].dd->lghosts[d];
      }
    }
    
    // Initialise current position
    reflevel  = 0;
    mglevel   = 0;
    component = -1;
    
    // Invent a refinement structure
    gh<dim>::rexts bbsss;
    gh<dim>::rprocs pss;
    MakeRegions_RefineCentre (cgh, maxreflevels, bbsss, pss);
    RegisterRecomposeRegions (bbsss, pss);
    
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
  
} // namespace Carpet
