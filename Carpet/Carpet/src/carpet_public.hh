// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.hh,v 1.20 2002/09/25 15:50:32 schnetter Exp $

// It is assumed that the number of components of all arrays is equal
// to the number of components of the grid functions, and that their
// distribution onto the processors is the same, and that all
// processors own the same number of components.

#ifndef CARPET_PUBLIC_HH
#define CARPET_PUBLIC_HH

#include <vector>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/th.hh"
  
#include "carpet_public.h"

  
  
namespace Carpet {
  
  
  
  const int dim = 3;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  extern int GHExtension;
  
  // Refinement factor
  extern int reffact;
  
  // Refinement factor on finest grid
  extern int maxreflevelfact;
  
  // Multigrid levels
  extern int mglevels;
  
  // Multigrid factor
  extern int mgfact;
  
  // Maximum number of refinement levels
  extern int maxreflevels;
  
  // Multigrid factor on coarsest grid
  extern int maxmglevelfact;
  
  // Current iteration per refinement level
  extern vector<int> iteration;
  
  // Current position on the grid hierarchy
  extern int reflevel;
  extern int mglevel;
  extern int component;
  
  // Current refinement factor
  extern int reflevelfact;
  
  // Current multigrid factor
  extern int mglevelfact;
  
  // Time step on base grid
  extern CCTK_REAL base_delta_time;
  
  // Current time
  extern CCTK_REAL current_time;
  
  
  
  // Data for grid functions
  
  // The grid hierarchy
  extern gh<dim>* hh;
  extern th* tt;
  extern dh<dim>* dd;
  
  // Data for everything
  struct arrdesc {
    // points to hh etc. for GF, and is unique for SCALAR and ARRAY
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
  int Evolve (tFleshConfig* config);
  int Shutdown (tFleshConfig* config);
  int CallFunction (void* function, cFunctionData* attribute, void* data);
  
  int SyncGroup (cGH* cgh, const char* groupname);
  int EnableGroupStorage (cGH* cgh, const char* groupname);
  int DisableGroupStorage (cGH* cgh, const char* groupname); 
  int EnableGroupComm (cGH* cgh, const char* groupname);
  int DisableGroupComm (cGH* cgh, const char* groupname);
  int Barrier (const cGH* cgh);
  int Exit (cGH* cgh, int retval);
  int Abort (cGH* cgh, int retval);
  int MyProc (const cGH* cgh);
  int nProcs (const cGH* cgh);
  const int* ArrayGroupSizeB (const cGH* cgh, int dir, int group,
			      const char* groupname);
  int QueryGroupStorageB (const cGH* cgh, int group, const char* groupname);
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data);
  
  
  
  // Functions for recomposing the grid hierarchy
  void RegisterRegridRoutine (int (*routine)(const cGH * cckgGH,
					     gh<dim>::rexts& bbsss,
					     gh<dim>::rbnds& obss,
					     gh<dim>::rprocs& pss));
  
  void SplitRegions (const cGH* cgh, vector<bbox<int,dim> >& bbs,
		     vector<vect<vect<bool,2>,dim> >& obs);
  void SplitRegions_AlongZ (const cGH* cgh, vector<bbox<int,dim> >& bbs,
			    vector<vect<vect<bool,2>,dim> >& obs);
  
  void MakeProcessors (const cGH* cgh, const gh<dim>::rexts& bbsss,
		       gh<dim>::rprocs& pss);
  
  
  
  // Helper functions
  void set_reflevel (cGH* cgh, int rl);
  void set_mglevel (cGH* cgh, int ml);
  void set_component (cGH* cgh, int c);
  
  
  
  // Refinement level iterator
  
#define BEGIN_REFLEVEL_LOOP(cgh)			\
  do {							\
    int _rll;						\
    assert (reflevel==0);				\
    for (int _rl=0; _rl<hh->reflevels(); ++_rl) {	\
      set_reflevel ((cGH*)(cgh), _rl);			\
      {
#define END_REFLEVEL_LOOP(cgh)			\
      }						\
    }						\
    set_reflevel ((cGH*)(cgh), 0);		\
    assert (reflevel==0);			\
    _rll = 0;					\
  } while (0)
  
  
  
  // Reverse refinement level iterator
  
#define BEGIN_REVERSE_REFLEVEL_LOOP(cgh)		\
  do {							\
    int _rrll;						\
    assert (reflevel==0);				\
    for (int _rl=hh->reflevels()-1; _rl>=0; --_rl) {	\
      set_reflevel ((cGH*)(cgh), _rl);			\
      {
#define END_REVERSE_REFLEVEL_LOOP(cgh)		\
      }						\
    }						\
    assert (reflevel==0);			\
    _rrll = 0;					\
  } while (0)
  
  
  
  // Multigrid level iterator
  
#define BEGIN_MGLEVEL_LOOP(cgh)				\
  do {							\
    int _mgl;						\
    assert (reflevel>=0 && reflevel<hh->reflevels());	\
    assert (mglevel==-1);				\
    for (int _ml=mglevels-1; _ml>=0; --_ml) {		\
      set_mglevel ((cGH*)(cgh), _ml);			\
      {
#define END_MGLEVEL_LOOP(cgh)			\
      }						\
    }						\
    set_mglevel ((cGH*)(cgh), -1);		\
    assert (mglevel==-1);			\
    _mgl = 0;					\
  } while (0)
  
  
  
  // Component iterator
  
#define BEGIN_COMPONENT_LOOP(cgh)			\
  do {							\
    int _cl;						\
    assert (reflevel>=0 && reflevel<hh->reflevels());	\
    assert (mglevel>=0 && mglevel<mglevels);		\
    assert (component==-1);				\
    for (int _c=0; _c<hh->components(reflevel); ++_c) {	\
      set_component ((cGH*)(cgh), _c);			\
      {
#define END_COMPONENT_LOOP(cgh)			\
      }						\
    }						\
    set_component ((cGH*)(cgh), -1);		\
    assert (component==-1);			\
    _cl = 0;					\
  } while (0)


  
#define BEGIN_LOCAL_COMPONENT_LOOP(cgh)			\
  do {							\
    int _lcl;						\
    assert (reflevel>=0 && reflevel<hh->reflevels());	\
    assert (mglevel>=0 && mglevel<mglevels);		\
    assert (component==-1);				\
    for (int _c=0; _c<hh->components(reflevel); ++_c) {	\
      if (hh->is_local(reflevel,_c)) {			\
	set_component ((cGH*)(cgh), _c);		\
	{
#define END_LOCAL_COMPONENT_LOOP(cgh)		\
	}					\
      }						\
    }						\
    set_component ((cGH*)(cgh), -1);		\
    assert (component==-1);			\
    _lcl = 0;					\
  } while (0)
  
} // namespace Carpet

#endif // !defined(CARPET_PUBLIC_HH)
