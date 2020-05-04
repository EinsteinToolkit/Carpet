#ifndef PRESYNC_HH
#define PRESYNC_HH
#include "cctk.h"
#include "cctk_Schedule.h"
#include "carpet.hh"
#include <set>
#include <vector>


namespace Carpet {
extern void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::set<int>& pregroups);
extern void PreCheckValid(cFunctionData *attribute,cGH *cctkGH,std::set<int>& pregroups);
extern void PostCheckValid(cFunctionData *attribute,cGH *cctkGH,const std::vector<int>& sync_groups);

void ApplyPhysicalBCsForGroupI(const cGH *cctkGH, const int group_index);
void ApplyPhysicalBCsForVarI(const cGH *cctkGH, const int var_index);
}

#endif
