// It is assumed that each group has at least one map.  All arrays
// have exactly one map.  All maps have the same number of refinement
// levels.

// It is assumed that each group has at least one component.

// It is assumed that the number of components of all arrays is equal
// to the number of components of the grid functions, and that their
// distribution onto the processors is the same, and that all
// processors own the same number of components.

#ifndef VARIABLES_HH
#define VARIABLES_HH

#include <vector>

#include "cctk.h"

#include "data.hh"
#include "dh.hh"
#include "ggf.hh"
#include "gh.hh"
#include "operators.hh"
#include "th.hh"
#include "vect.hh"

#include "carpet_public.h"
#include "defines.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  extern int GHExtension;
  
  // Maximum number of refinement levels
  extern int maxreflevels;
  
  // Refinement levels
  extern int reflevels;
  
  // Refinement factor
  extern int reffact;
  
  // Refinement factor on finest possible grid
  extern int maxreflevelfact;
  
  // Base multigrid level
  extern int basemglevel;
  
  // Multigrid levels
  extern int mglevels;
  
  // Multigrid factor
  extern int mgfact;
  
  // Multigrid factor on coarsest grid
  extern int maxmglevelfact;
  
  // Maps
  extern int maps;
  
  
  
  // Current position on the grid hierarchy
  extern int reflevel;
  extern int mglevel;
  extern int map;
  extern int component;
  
  // Current refinement factor
  extern int reflevelfact;
  
  // Current multigrid factor
  extern int mglevelfact;
  
  
  
  // Carpet's GH
  extern CarpetGH carpetGH;
  
  
  
  // Times and spaces on the refinement levels
  extern CCTK_REAL global_time;
  extern vector<vector<CCTK_REAL> > leveltimes; // [mglevel][reflevel]
  extern CCTK_REAL delta_time;
  
  extern vector<vect<CCTK_REAL,dim> > origin_space; // [mglevel]
  extern vect<CCTK_REAL,dim> delta_space;
  
  
  
  // Is this the time for a global mode call?
  extern bool do_meta_mode;
  extern bool do_global_mode;
  
  // Is prolongation enabled?
  extern bool do_prolongate;
  
  
  
  // Data for grid functions
  
  // The grid hierarchy
  extern vector<gh*> vhh;       // [map]
  extern vector<dh*> vdd;       // [map]
  extern vector<th*> vtt;       // [map]
  
  // Data for the groups
  struct groupdesc {
    cGroupDynamicData info;
    operator_type transport_operator; // prolongation and restriction
  };
  extern vector<groupdesc> groupdata; // [group]
  
  // structure to hold a set of groups which all have the same CCTK vartype
  struct group_set {
    int vartype;                // eg. CCTK_VARIABLE_REAL, etc.
    vector<int> members;        // members of this set
  };                            // (given by their CCTK group indices)

  // Data for everything
  struct arrdesc {
    // points to hh etc. for GF, and is unique for SCALAR and ARRAY
    gh* hh;
    dh* dd;
    th* tt;
    vector<ggf*> data;          // [var]
  };
  extern vector<vector<arrdesc> > arrdata; // [group][map]
  
} // namespace Carpet

#endif // !defined(VARIABLES_HH)
