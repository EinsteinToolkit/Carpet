#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void carpettest_check_sizes (CCTK_ARGUMENTS);

static void print_scalar (const char *name, int sc);
static void print_scalar_descr (const char *name, int sc, const char *descr);
static void print_array (const char *name, int dim, const int *arr);

static const char *grouptype_string (int grouptype);
static const char *disttype_string (int disttype);



void carpettest_check_sizes (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  int group;
  int dim;
  cGroup data;
  cGroupDynamicData dyndata;
  
  dim = cctk_dim;
  print_scalar ("cctk_dim", cctk_dim);
  print_array ("cctk_gsh", dim, cctk_gsh);
  print_array ("cctk_lsh", dim, cctk_lsh);
  print_array ("cctk_lbnd", dim, cctk_lbnd);
  print_array ("cctk_ubnd", dim, cctk_ubnd);
  print_array ("cctk_bbox", 2*dim, cctk_bbox);
  print_array ("cctk_nghostzones", dim, cctk_nghostzones);
  printf ("\n");
  
  for (group=0; group<CCTK_NumGroups(); ++group) {
    CCTK_GroupData (group, &data);
    CCTK_GroupDynamicData (cctkGH, group, &dyndata);
    
    print_scalar_descr ("group", group, CCTK_GroupName(group));
    
    dim = data.dim;
    print_scalar ("dim", data.dim);
    print_scalar_descr ("grouptype", data.grouptype, grouptype_string(data.grouptype));
    print_scalar_descr ("vartype", data.vartype, CCTK_VarTypeName(data.vartype));
    print_scalar_descr ("disttype", data.disttype, disttype_string(data.disttype));
    print_scalar ("stagtype", data.stagtype);
    print_scalar ("numvars", data.numvars);
    print_scalar ("numtimelevels", data.numtimelevels);
    print_array ("gsh", dim, dyndata.gsh);
    print_array ("lsh", dim, dyndata.lsh);
    print_array ("lbnd", dim, dyndata.lbnd);
    print_array ("ubnd", dim, dyndata.ubnd);
    print_array ("bbox", 2*dim, dyndata.bbox);
    print_array ("nghostzones", dim, dyndata.nghostzones);
    printf ("\n");
  }
  
}



static void print_scalar (const char *name, int sc)
{
  printf ("%-15s: %3d\n", name, sc);
}

static void print_scalar_descr (const char *name, int sc, const char *descr)
{
  printf ("%-15s: %3d   %s\n", name, sc, descr);
}

static void print_array (const char *name, int dim, const int *arr)
{
  int d;
  printf ("%-15s:", name);
  for (d=0; d<dim; ++d) {
    printf (" %3d", arr[d]);
  }
  printf ("\n");
}
 
 

static const char *grouptype_string (int grouptype)
{
  switch (grouptype) {
  case CCTK_SCALAR: return "CCTK_SCALAR";
  case CCTK_GF:     return "CCTK_GF";
  case CCTK_ARRAY:  return "CCTK_ARRAY";
  }
  return "[illegal group type]";
}

static const char *disttype_string (int disttype)
{
  switch (disttype) {
  case CCTK_DISTRIB_CONSTANT: return "CCTK_DISTRIB_CONSTANT";
  case CCTK_DISTRIB_DEFAULT:  return "CCTK_DISTRIB_DEFAULT";
  }
  return "[illegal distribution type]";
}
