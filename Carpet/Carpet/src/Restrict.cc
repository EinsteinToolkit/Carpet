#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "cctk.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Restrict.cc,v 1.22 2003/11/05 16:18:37 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Restrict_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  void Restrict (const cGH* cgh)
  {
    assert (reflevel != -1);
    assert (mglevel != -1);
    assert (component == -1);
    
    Checkpoint ("%*sRestrict", 2*reflevel, "");
    
    // Restrict
    if (reflevel < hh->reflevels()-1) {
      for (comm_state<dim> state; !state.done(); state.step()) {
        for (int group=0; group<CCTK_NumGroups(); ++group) {
          if (CCTK_GroupTypeI(group) == CCTK_GF) {
            if (CCTK_QueryGroupStorageI(cgh, group)) {
              if (arrdata[group].do_transfer) {
                for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
                  
                  const int tl = 0;
                  
                  // use background time here (which may not be
                  // modified by the user)
                  const CCTK_REAL time = tt->time (tl, reflevel, mglevel);
                  const CCTK_REAL time1 = tt->time (0, reflevel, mglevel);
                  const CCTK_REAL time2 = cgh->cctk_time / cgh->cctk_delta_time;
                  assert (fabs((time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time))) < 1e-12);
                  
                  
                  for (int c=0; c<hh->components(reflevel); ++c) {
                    arrdata[group].data[var]->ref_restrict
                      (state, tl, reflevel, c, mglevel, time);
                  }
                  
                } // loop over variables
              } else {
                if (state.thestate==state_recv) {
                  char * groupname = CCTK_GroupName(group);
                  Checkpoint ("%*s(no restricting for group %s)",
                              2*reflevel, "", groupname);
                  free (groupname);
                }
              } // if ! do_transfer
            } // if group has storage
          } // if grouptype == CCTK_GF
        } // loop over groups
      } // for state
    } // if not finest refinement level
    
    
    
    // Sync
    if (reflevel < hh->reflevels()-1) {
      for (comm_state<dim> state; !state.done(); state.step()) {
        for (int group=0; group<CCTK_NumGroups(); ++group) {
          if (CCTK_GroupTypeI(group) == CCTK_GF) {
            if (CCTK_QueryGroupStorageI(cgh, group)) {
              if (arrdata[group].do_transfer) {
                for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
                  
                  const int tl = 0;
                  
                  for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
                    arrdata[group].data[var]->sync
                      (state, tl, reflevel, c, mglevel);
                  }
                  
                } // loop over variables
              } // if do_transfer
            } // if group has storage
          } // if grouptype == CCTK_GF
        } // loop over groups
      } // for state
    } // if not finest refinement level
    
  }
  
} // namespace Carpet
