#include <vector>
#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/th.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/variables.cc,v 1.3 2001/12/17 13:34:02 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  int GHExtension;
  
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
  
  // refinement factor of current level: pow(refinement_factor, reflevel)
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
