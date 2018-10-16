#ifndef CACHEINFO_HH
#define CACHEINFO_HH

#include <cctk.h>

#include <cassert>
#include <limits>

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"

namespace CarpetLib {

// <http://graphics.stanford.edu/~seander/bithacks.html>
template <typename T> static T next_power_of_2(T x) CCTK_ATTRIBUTE_UNUSED;
template <typename T> static T next_power_of_2(T x) {
  assert(x > 0);
  // T res = 1;
  // while (res < x) {
  //   res *= 2;
  //   assert(res > 0); // try to catch overflows
  // }
  // return res;
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  if (sizeof(T) > 1)
    x |= x >> 8;
  if (sizeof(T) > 2)
    x |= x >> 16;
  if (sizeof(T) > 4)
    x |= x >> 32;
  ++x;
  return x;
}

template <typename T>
static T previous_power_of_2(T const x) CCTK_ATTRIBUTE_UNUSED;
template <typename T> static T previous_power_of_2(T const x) {
  return next_power_of_2((x >> 1) + 1);
}

// <http://graphics.stanford.edu/~seander/bithacks.html>
template <typename T>
static bool is_power_of_2(T const x) CCTK_ATTRIBUTE_UNUSED;
template <typename T> static bool is_power_of_2(T const x) {
  assert(x > 0);
  // return x == next_power_of_2(x);
  return (x & (x - 1)) == 0;
}

template <typename T, int D> struct padding_t {
  vect<T, D> padded_shape;      // padded extent
  vect<T, D> padding_alignment; // padding alignment
  vect<T, D> padding_offset;    // padding offset
  // forall j, k: &var[0,j,k] % padding_alignment == padding_offset
};

// These routines are apparently not pure -- don't know why
template <int D>
padding_t<int, D> pad_shape(bbox<int, D> const &extent,
                            bbox<int, D> const &owned) /*CCTK_ATTRIBUTE_PURE*/;
template <int D>
padding_t<int, D> pad_shape(vect<int, D> const &shape,
                            vect<int, D> const &offset) /*CCTK_ATTRIBUTE_PURE*/;
}

#endif // CACHEINFO_HH
