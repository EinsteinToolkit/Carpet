#ifndef CACHEINFO_HH
#define CACHEINFO_HH

#include <cctk.h>

#include <limits>

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"



static ptrdiff_t next_power_of_2 (ptrdiff_t const x)
  CCTK_ATTRIBUTE_UNUSED;
static ptrdiff_t next_power_of_2 (ptrdiff_t const x)
{
  assert (x > 0);
  ptrdiff_t res = 1;
  while (res < x) {
    res *= 2;
    assert (res > 0);           // try to catch overflows
  }
  return res;
}

static ptrdiff_t previous_power_of_2 (ptrdiff_t const x)
  CCTK_ATTRIBUTE_UNUSED;
static ptrdiff_t previous_power_of_2 (ptrdiff_t const x)
{
  return next_power_of_2(x/2+1);
}

static bool is_power_of_2 (ptrdiff_t const x)
  CCTK_ATTRIBUTE_UNUSED;
static bool is_power_of_2 (ptrdiff_t const x)
{
  return x == next_power_of_2(x);
}

#if 0
static ptrdiff_t lcm (ptrdiff_t const x, ptrdiff_t const y)
  CCTK_ATTRIBUTE_UNUSED;
static ptrdiff_t lcm (ptrdiff_t const x, ptrdiff_t const y)
{
  assert (x > 0 && y > 0);
  ptrdiff_t z = x;
  // TODO: improve LCM algorithm
  while (z % y != 0) z += x;
  assert (z % x == 0 && z % y == 0);
  return z;
}
#endif



// These routines are apparently not pure -- don't know why
template<int D>
vect<int,D> pad_shape (bbox<int,D> const& extent) /*CCTK_ATTRIBUTE_PURE*/;
template<int D>
vect<int,D> pad_shape (vect<int,D> const& shape) /*CCTK_ATTRIBUTE_PURE*/;

#endif  // CACHEINFO_HH