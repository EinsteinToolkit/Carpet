// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.hh,v 1.4 2001/11/02 16:05:03 schnetter Exp $

// It is assumed that the number of components of all arrays is equal
// to the number of components of the grid functions, and that their
// distribution onto the processors is the same, and that all
// processors own the same number of components.

#include <vector>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/th.hh"

namespace Carpet {
  
  
  
  const int dim = 3;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  extern int GHExtension;
  
  // Maximum number of refinement levels
  extern int maxreflevels;
  
  // Refinement factor on finest grid
  extern int maxreflevelfact;
  
  // Current iteration per refinement level
  extern vector<int> iteration;
  
  // Current position on the grid hierarchy
  extern int reflevel;
  extern int mglevel;
  extern int component;
  
  // Current refinement factor
  extern int reflevelfact;
  
  // Time step on base grid
  extern CCTK_REAL base_delta_time;
  
  
  
  // Data for grid functions
  
  // The grid hierarchy
  extern gh<dim>* hh;
  extern th* tt;
  extern dh<dim>* dd;
  
  // Data for scalars
  extern gh<dim>* hh0;
  extern th* tt0;
  extern dh<dim>* dd0;
  
  // Data for everything
  struct arrdesc {
    // points to hh etc. for GF, and to hh0 etc. for SCALAR
    // is unique for ARRAY
    cGroupDynamicData info;
    gh<dim>* hh;
    th* tt;
    dh<dim>* dd;
    vector<generic_gf<dim>* > data; // [var]
  };
  extern vector<arrdesc> arrdata; // [group]
  
  
  
  // Checksums
  struct ckdesc {
    bool valid;
    int sum;
  };
  extern vector<vector<vector<vector<ckdesc> > > > checksums; // [n][rl][tl][c]
  
  
  
  // Stuff with C linkage
  extern "C" {
#include "carpet_public.h"
  }
  
  
  
  // Registered functions
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh);
  
  int Initialise (tFleshConfig* config);
  int Evolve (const tFleshConfig* config);
  int Shutdown (const tFleshConfig* config);
  int CallFunction (void* function, cFunctionData* attribute, void* data);
  
  int SyncGroup (const cGH* cgh, const char* groupname);
  int EnableGroupStorage (cGH* cgh, const char* groupname);
  int DisableGroupStorage (cGH* cgh, const char* groupname); 
  int EnableGroupComm (const cGH* cgh, const char* groupname);
  int DisableGroupComm (const cGH* cgh, const char* groupname);
  int Barrier (const cGH* cgh);
  int Exit (const cGH* cgh, int retval);
  int Abort (const cGH* cgh, int retval);
  int MyProc (const cGH* cgh);
  int nProcs (const cGH* cgh);
  const int* ArrayGroupSizeB (const cGH* cgh, int dir, int group,
			      const char* groupname);
  int QueryGroupStorageB (const cGH* cgh, int group, const char* groupname);
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data);
  
  
  
  // Functions for recomposing the grid hierarchy
  void RegisterRecomposeRegions (const gh<dim>::rexts& bbsss,
				 const gh<dim>::rprocs& pss);
  
  void MakeRegions_RefineCentre (const cGH* cgh, int reflevels,
				 gh<dim>::rexts& bbsss, gh<dim>::rprocs& pss);
  
  
  
  // Helper functions
  void set_reflevel (cGH* cgh, int rl);
  void set_mglevel (cGH* cgh, int ml);
  void set_component (cGH* cgh, int c);
  
  
  
  // Refinement level iterator
  
#define BEGIN_REFLEVEL_LOOP(cgh)		\
  do {						\
    int _rl;					\
    assert (reflevel==0);			\
    for (;;) {					\
      {
#define END_REFLEVEL_LOOP(cgh)			\
      }						\
      if (reflevel==maxreflevels-1) break;	\
      set_reflevel ((cgh), reflevel+1);		\
    }						\
    set_reflevel ((cgh), 0);			\
    assert (reflevel==0);			\
    _rl = 0;					\
  } while (0)
  
  
  
  // Reverse refinement level iterator
  
#define BEGIN_REVERSE_REFLEVEL_LOOP(cgh)	\
  do {						\
    int _rrl;					\
    assert (reflevel==0);			\
    set_reflevel ((cgh), maxreflevels-1);	\
    for (;;) {					\
      {
#define END_REVERSE_REFLEVEL_LOOP(cgh)		\
      }						\
      if (reflevel==0) break;			\
      set_reflevel ((cgh), reflevel-1);		\
    }						\
    assert (reflevel==0);			\
    _rrl = 0;					\
  } while (0)
  
  
  
  // Component iterator
  
#define BEGIN_COMPONENT_LOOP(cgh)			\
  do {							\
    int _cl;						\
    assert (reflevel>=0 && reflevel<hh->reflevels());	\
    assert (component==-1);				\
    set_component ((cgh), 0);				\
    for (;;) {						\
      {
#define END_COMPONENT_LOOP(cgh)				\
      }							\
      if (component==hh->components(reflevel)-1) break;	\
      set_component ((cgh), component+1);		\
    }							\
    set_component ((cgh), -1);				\
    assert (component==-1);				\
    _cl = 0;						\
  } while (0)
  
} // namespace Carpet
