#ifndef PRESYNC_HH
#define PRESYNC_HH


extern "C" void ShowValid();
extern "C" void diagnosticPreValid();

extern "C" void attempt_readwrite(const char *thorn,const char *var, int spec);
extern "C" void diagnosticChanged();
extern "C" void TraverseReads(const char *func_name,void(*trace_func)(int,int,int));
extern "C" void TraverseReads(const char *func_name,void(*trace_func)(int,int,int));
extern "C" void TraverseWrites(const char *func_name,void(*trace_func)(int,int,int));

namespace Carpet {
extern "C" void clear_readwrites();
extern "C" void check_readwrites();
extern "C" int Carpet_hasAccess(const cGH *cctkGH,int var_index);
extern "C" void Carpet_requestAccess(int var_index,int read_spec,int write_spec);
extern "C" void SetValidRegion(int vi,int tl,int wh);
extern "C" int GetValidRegion(int vi,int tl);
}

#endif
