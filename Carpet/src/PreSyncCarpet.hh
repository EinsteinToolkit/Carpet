#ifndef PRESYNC_HH
#define PRESYNC_HH
#include "cctk.h"
#include "cctk_Schedule.h"
#include "carpet.hh"
#include <set>
#include <vector>


namespace Carpet {
extern "C" void clear_readwrites();
extern "C" void check_readwrites(cFunctionData const * const attribute);
extern "C" void SetValidRegion(int vi,int tl,int wh);
extern "C" int GetValidRegion(int vi,int tl);

extern void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::set<int>& pregroups);
extern void PreCheckValid(cFunctionData *attribute,cGH *cctkGH,std::set<int>& pregroups);
extern void PostCheckValid(cFunctionData *attribute,cGH *cctkGH,const std::vector<int>& sync_groups);

void ApplyPhysicalBCsForGroupI(const cGH *cctkGH, const int group_index);
void ApplyPhysicalBCsForVarI(const cGH *cctkGH, const int var_index);

extern void invalidate_rdwr(const cGH *cctkGH, int vi, int tl);
}

#endif
