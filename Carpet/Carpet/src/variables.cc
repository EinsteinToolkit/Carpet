
#include "variables.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  int GHExtension;
  
  // Maximum number of refinement levels
  int maxreflevels;
  
  // Refinement levels
  int reflevels;
  
#if 0
  // Refinement factor
  int reffact;
#endif
  
  // Temporal refinement factors over the coarsest grid
  vector<int> timereffacts;
  
  // Spatial refinement factors over the coarsest grid
  vector<vect<int,dim> > spacereffacts;
  
  // Maximum refinement factors on finest possible grid
  int maxtimereflevelfact;
  vect<int,dim> maxspacereflevelfact;
  
  // Base multigrid level
  int basemglevel;
  
  // Multigrid levels
  int mglevels;
  
  // Multigrid factor
  int mgfact;
  
  // Multigrid factor on coarsest grid
  int maxmglevelfact;
  
  // Maps
  int maps;
  
  
  
  // Current position on the grid hierarchy
  int reflevel;
  int mglevel;
  int map;
  int component;
  
  // Current refinement factors
  int timereflevelfact;
  vect<int,dim> spacereflevelfact;
  
  // Current multigrid factor
  int mglevelfact;
  
  
  
  // Carpet's GH
  CarpetGH carpetGH;
  
  
  
  // Times and spaces on the refinement levels
  CCTK_REAL global_time;
  vector<vector<CCTK_REAL> > leveltimes; // [mglevel][reflevel]
  CCTK_REAL delta_time;
  
  vector<vect<CCTK_REAL,dim> > origin_space; // [mglevel]
  vect<CCTK_REAL,dim> delta_space;
  
  
  
  // Is this the time for a global mode call?
  bool do_meta_mode;
  bool do_global_mode;
  
  // Is prolongation enabled?
  bool do_prolongate;
  
  
  
  // Data for grid functions
  
  // The grid hierarchy
  vector<gh*> vhh;              // [map]
  vector<dh*> vdd;              // [map]
  vector<th*> vtt;              // [map]

  // Data for the groups
  vector<groupdesc> groupdata;  // [group]
  
  // Data for everything
  vector<vector<arrdesc> > arrdata; // [group][map]
  
} // namespace Carpet
