#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "cctk.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Restrict.cc,v 1.19 2003/07/20 21:03:43 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Restrict_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  void Restrict (const cGH* cgh)
  {
    assert (component == -1);
    
    Checkpoint ("%*sRestrict", 2*reflevel, "");
    
    // Loop over grid functions with storage
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF
	  && CCTK_QueryGroupStorageI(cgh, group)) {
        if (arrdata[group].do_transfer) {
          for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
            
            const int tl = 0;
            
            // use background time here (which may not be modified by
            // the user)
            const CCTK_REAL time = tt->time (tl, reflevel, mglevel);
            if (tl==0) {
              const CCTK_REAL time1 = tt->time (tl, reflevel, mglevel);
              const CCTK_REAL time2 = cgh->cctk_time / cgh->cctk_delta_time;
              assert (fabs((time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time))) < 1e-12);
            }
            
            if (reflevel < hh->reflevels()-1) {
              
              for (int c=0; c<hh->components(reflevel); ++c) {
                arrdata[group].data[var]->ref_restrict
                  (tl, reflevel, c, mglevel, time);
              }
              for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
                arrdata[group].data[var]->sync (tl, reflevel, c, mglevel);
              }
              
            } // if not finest refinement level
            
          } // loop over variables
        } else {
          char * groupname = CCTK_GroupName(group);
          Checkpoint ("%*s(no restricting for group %s)",
                      2*reflevel, "", groupname);
          free (groupname);
        }
      } // if group has storage
    } // loop over groups
    
  }
  
} // namespace Carpet
