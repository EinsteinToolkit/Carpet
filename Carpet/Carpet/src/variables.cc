#include <vector>
#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/th.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/variables.cc,v 1.5 2002/01/09 17:45:40 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  int GHExtension;
  
  // Multigrid levels
  int mglevels;
  
  // Multigrid factor
  int mgfact;
  
  // Multigrid factor on coarsest grid
  int maxmglevelfact;
  
  // Maximum number of refinement levels
  int maxreflevels;
  
  // Refinement factor
  int reffact;
  
  // Refinement factor on finest grid
  int maxreflevelfact;
  
  // Current iteration per refinement level
  vector<int> iteration;
  
  // Current position on the grid hierarchy
  int mglevel;
  int reflevel;
  int component;
  
  // multigrid factor of current level: ipow(multigrid_factor, mglevel)
  int mglevelfact;
  
  // refinement factor of current level: ipow(refinement_factor, reflevel)
  int reflevelfact;
  
  // Time step on base grid
  CCTK_REAL base_delta_time;
  
  
  
  // Data for grid functions
  
  // The grid hierarchy
  gh<dim>* hh;
  th* tt;
  dh<dim>* dd;
  
  // Data for scalars
  gh<dim>* hh0;
  th* tt0;
  dh<dim>* dd0;

    // Data for everything
  vector<arrdesc> arrdata;	// [group]
  
  // Checksums
  vector<vector<vector<vector<ckdesc> > > > checksums; // [n][rl][tl][c]
  
} // namespace Carpet
