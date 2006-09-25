#include <assert.h>
#include <string.h>

#include <cctk.h>



void
CCTK_FCALL
CCTK_FNAME(checkindex) (int const * restrict const i,
                        int const * restrict const j,
                        int const * restrict const k,
                        int const * restrict const di,
                        int const * restrict const dj,
                        int const * restrict const dk,
                        int const * restrict const imax,
                        int const * restrict const jmax,
                        int const * restrict const kmax,
                        ONE_FORTSTRING_ARG)
{
  if (*i < 1 || *i+*di-1 > *imax ||
      *j < 1 || *j+*dj-1 > *jmax ||
      *k < 1 || *k+*dk-1 > *kmax)
  {
    ONE_FORTSTRING_CREATE (where);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "%s array index out of bounds: shape is (%d,%d,%d), index is (%d,%d,%d), extent is (%d,%d,%d)",
                where, *imax,*jmax,*kmax, *i,*j,*k, *di,*dj,*dk);
    assert (0);
    free (where);
  }
}
