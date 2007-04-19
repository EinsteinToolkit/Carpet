#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



static void
enforce (int const varindex, char const * const optstring, void * const arg);



void
enforce_evolution_mask(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_TraverseString (enforce_vars, enforce, cctkGH, CCTK_GROUP_OR_VAR);
}



void
enforce(int const varindex, char const * const optstring, void * const arg)
{

  /* enforce the rhs of a given var to zero where the evolution mask
     says so */

  cGH const * const cctkGH = (cGH const *) arg;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  int rhsindex;
  int vargroup, rhsgroup;
  cGroup vardata, rhsdata;
  CCTK_REAL * rhsptr;

  int index;
  int i,j,k;
  int ierr;

  int nx = cctkGH->cctk_lsh[0];
  int ny = cctkGH->cctk_lsh[1];
  int nz = cctkGH->cctk_lsh[2];


  assert (varindex >= 0);

  /* make sure we do have a rhs for this var */

  rhsindex = MoLQueryEvolvedRHS (varindex);
  if (rhsindex < 0) {
    char * const fullvarname = CCTK_FullName (varindex);
    assert (fullvarname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "There is no RHS variable registered with MoL for the evolved variable \"%s\"",
                fullvarname);
    free (fullvarname);
  }
  

 if (verbose) {
    char * const fullvarname = CCTK_FullName (varindex);
    char * const fullrhsname = CCTK_FullName (rhsindex);
    assert (fullvarname);
    assert (fullrhsname);
    CCTK_VInfo (CCTK_THORNSTRING,
                "Enforcing evolution mask on \"%s\" (RHS \"%s\")",
                fullvarname, fullrhsname);
    free (fullvarname);
    free (fullrhsname);
  }


  assert (rhsindex >= 0);
  
  vargroup = CCTK_GroupIndexFromVarI (varindex);
  assert (vargroup >= 0);
  rhsgroup = CCTK_GroupIndexFromVarI (rhsindex);
  assert (rhsgroup >= 0);

  ierr = CCTK_GroupData (vargroup, &vardata);
  assert (!ierr);
  ierr = CCTK_GroupData (rhsgroup, &rhsdata);
  assert (!ierr);

  assert (vardata.grouptype == CCTK_GF);
  assert (vardata.vartype == CCTK_VARIABLE_REAL);
  assert (vardata.dim == cctk_dim);
  assert (rhsdata.grouptype == CCTK_GF);
  assert (rhsdata.vartype == CCTK_VARIABLE_REAL);
  assert (rhsdata.dim == cctk_dim);

  rhsptr = CCTK_VarDataPtrI (cctkGH, 0, rhsindex);
  assert (rhsptr);



  if(!writeNaNs) {
    
    for (k = 0; k<nz; k++)
      for (j = 0; j<ny; j++)
	for (i = 0; i<nx; i++) {
	  
	  index = CCTK_GFINDEX3D (cctkGH, i, j, k);
	  rhsptr[index] = evolution_mask[index]*rhsptr[index];

	}

  } else {

    for (k = 0; k<nz; k++)
      for (j = 0; j<ny; j++)
	for (i = 0; i<nx; i++) {
	  
	  index = CCTK_GFINDEX3D (cctkGH, i, j, k);
	  if(evolution_mask[index] < 1e-10) {
	    rhsptr[index] = 1.0e+12;
	  }
          
	}
    
  }


}
