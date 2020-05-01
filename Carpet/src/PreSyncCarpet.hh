#ifndef PRESYNC_HH
#define PRESYNC_HH
#include <cctk.h>
#include <cctk_Arguments.h>
#include <carpet.hh>
#include <set>
#include <vector>


extern "C" void ShowValid();
extern "C" void diagnosticPreValid();

extern "C" void attempt_readwrite(const char *thorn,const char *var, int spec);
extern "C" void diagnosticChanged();

namespace Carpet {
extern "C" void clear_readwrites();
extern "C" void check_readwrites();
extern "C" void SetValidRegion(int vi,int tl,int wh);
extern "C" int GetValidRegion(int vi,int tl);

extern void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::set<int>& pregroups);
extern void PreCheckValid(cFunctionData *attribute,cGH *cctkGH,std::set<int>& pregroups);
extern void PostCheckValid(cFunctionData *attribute,cGH *cctkGH,const std::vector<int>& sync_groups);

extern void cycle_rdwr(const cGH *cctkGH);
extern void uncycle_rdwr(const cGH *cctkGH);
extern void flip_rdwr(const cGH *cctkGH, int vi);
extern void invalidate_rdwr(const cGH *cctkGH, int vi, int tl);
}

#endif
