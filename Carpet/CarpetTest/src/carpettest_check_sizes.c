/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetTest/src/carpettest_check_sizes.c,v 1.1 2001/07/04 12:29:54 schnetter Exp $ */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void carpettest_check_sizes (CCTK_ARGUMENTS);



void carpettest_check_sizes (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  char msg[1000];
  int group;
  int var;
  int dim;
  CCTK_INT **sizes, **ghostsizes;
  int d;
  int gsh[3], lsh[3], lbnd[3], ubnd[3];
  const int *cip;
  int size[3];
  
  sprintf (msg, "gsh: %d %d %d\n", cctk_gsh[0], cctk_gsh[1], cctk_gsh[2]);
  printf (msg);
  sprintf (msg, "lsh: %d %d %d\n", cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
  printf (msg);
  sprintf (msg, "lbnd: %d %d %d\n", cctk_lbnd[0], cctk_lbnd[1], cctk_lbnd[2]);
  printf (msg);
  sprintf (msg, "ubnd: %d %d %d\n", cctk_ubnd[0], cctk_ubnd[1], cctk_ubnd[2]);
  printf (msg);
  printf ("\n");
  
#if 0
  for (group=0; group<CCTK_NumGroups(); ++group) {
    var = CCTK_FirstVarIndexI(group);
    sprintf (msg, "group: %d %s\n", group, CCTK_GroupNameFromVarI(var));
    printf (msg);
    
    dim = CCTK_GroupDimI(group);
    sprintf (msg, "dim: %d\n", CCTK_GroupDimI(group));
    printf (msg);
    
    sizes = CCTK_GroupSizesI(group);
    sprintf (msg, "sizes:");
    if (sizes) {
      for (d=0; d<dim; ++d) {
	sprintf (msg, "%s %d", msg, (*sizes)[d]);
      }
    } else {
      sprintf (msg, "%s (no sizes)", msg);
    }
    sprintf (msg, "%s\n", msg);
    printf (msg);
    
    ghostsizes = CCTK_GroupGhostsizesI(group);
    sprintf (msg, "ghostsizes:");
    if (ghostsizes) {
      for (d=0; d<dim; ++d) {
	sprintf (msg, "%s %d", msg, (*ghostsizes)[d]);
      }
    } else {
      sprintf (msg, "%s (no ghostsizes)", msg);
    }
    sprintf (msg, "%s\n", msg);
    printf (msg);
    
    CCTK_GroupgshGI (cctkGH, dim, gsh, group);
    sprintf (msg, "gsh:");
    for (d=0; d<dim; ++d) {
      sprintf (msg, "%s %d", msg, gsh[d]);
    }
    sprintf (msg, "%s\n", msg);
    printf (msg);
    
    CCTK_GrouplshGI (cctkGH, dim, lsh, group);
    sprintf (msg, "lsh:");
    for (d=0; d<dim; ++d) {
      sprintf (msg, "%s %d", msg, lsh[d]);
    }
    sprintf (msg, "%s\n", msg);
    printf (msg);
    
    CCTK_GrouplbndGI (cctkGH, dim, lbnd, group);
    sprintf (msg, "lbnd:");
    for (d=0; d<dim; ++d) {
      sprintf (msg, "%s %d", msg, lbnd[d]);
    }
    sprintf (msg, "%s\n", msg);
    printf (msg);
    
    CCTK_GroupubndGI (cctkGH, dim, ubnd, group);
    sprintf (msg, "ubnd:");
    for (d=0; d<dim; ++d) {
      sprintf (msg, "%s %d", msg, ubnd[d]);
    }
    sprintf (msg, "%s\n", msg);
    printf (msg);
    
    printf ("\n");
  }
#endif
  
  for (group=0; group<CCTK_NumGroups(); ++group) {
    var = CCTK_FirstVarIndexI(group);
    sprintf (msg, "group: %d %s\n", group, CCTK_GroupNameFromVarI(var));
    printf (msg);
    
    dim = CCTK_GroupDimI(group);
    sprintf (msg, "dim: %d\n", CCTK_GroupDimI(group));
    printf (msg);
    
    for (d=0; d<dim; ++d) {
      cip = CCTK_ArrayGroupSizeI(cctkGH, d, group);
      assert (cip);
      size[d] = *cip;
    }
    
    sprintf (msg, "size:");
    for (d=0; d<dim; ++d) {
      sprintf (msg, "%s %d", msg, size[d]);
    }
    sprintf (msg, "%s\n", msg);
    printf (msg);
    
    printf ("\n");
  }
  
}
