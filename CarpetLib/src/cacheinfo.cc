#include "cacheinfo.hh"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <utility>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <vectors.h>

namespace CarpetLib {

template <int D>
padding_t<int, D> pad_shape(bbox<int, D> const &extent,
                            bbox<int, D> const &owned) {
  // This is an empty grid group
  // TODO: Use larger alignment and/or make up an offset
  if (extent.empty())
    return {extent.sizes(), vect<int, D>(1), vect<int, D>(0)};
  assert(not extent.empty());
  assert(not owned.empty());
  assert(owned.is_aligned_with(extent));
  assert(owned <= extent);
  return pad_shape(extent.sizes(),
                   (owned.lower() - extent.lower()) / extent.stride());
}

namespace {
struct cache_info_t {
  int type;
  int linesize;
  int stride;
};
bool have_cache_info = false;
vector<cache_info_t> cache_info;
}

template <int D>
padding_t<int, D> pad_shape(vect<int, D> const &shape,
                            vect<int, D> const &offset) {
  DECLARE_CCTK_PARAMETERS;

  assert(all(shape >= 0));
  if (not(all(offset >= 0 and offset <= shape)))
    cerr << "pad_shape: shape=" << shape << " offset=" << offset << "\n";
  assert(all(offset >= 0 and offset <= shape));

  // Don't pad empty arrays; we don't want to handle all the special
  // cases for this below
  if (any(shape == 0))
    return {shape, vect<int, D>(1), vect<int, D>(0)};

  if (CCTK_BUILTIN_EXPECT(not have_cache_info, false)) {
#pragma omp barrier
#pragma omp master
    {
      if (CCTK_IsFunctionAliased("GetCacheInfo1")) {
        int const num_levels =
            GetCacheInfo1(NULL, NULL, NULL, NULL, NULL, NULL, 0);
        vector<CCTK_INT> types(num_levels);
        vector<CCTK_INT> linesizes(num_levels);
        vector<CCTK_INT> strides(num_levels);
        GetCacheInfo1(NULL, &types[0], NULL, &linesizes[0], &strides[0], NULL,
                      num_levels);
        cache_info.resize(num_levels);
        for (int level = 0; level < num_levels; ++level) {
          cache_info[level].type = types[level];
          cache_info[level].linesize = linesizes[level];
          cache_info[level].stride = strides[level];
        }
      }
      have_cache_info = true;
    }
#pragma omp barrier
  }

  vect<int, D> padded_shape;
  vect<int, D> padding_alignment;
  vect<int, D> padding_offset;
  size_type accumulated_npoints = 1;
  for (int d = 0; d < D; ++d) {
    size_type npoints = shape[d];
    size_type bndsize = offset[d];
    size_type alignment = 1;

    if (VECTORISE and VECTORISE_ALIGNED_ARRAYS) {

      if (d == 0) {
        // Pad array to a multiple of the vector size. Note that this is
        // a hard requirement, so that we can emit aligned load/store
        // operations.
        if (VECTORISE_ALIGN_INTERIOR) {
          size_type newbndsize =
              align_up(bndsize, size_type(CCTK_REAL_VEC_SIZE));
          npoints += newbndsize - bndsize;
          bndsize = newbndsize;
        }
        npoints = align_up(npoints, size_type(CCTK_REAL_VEC_SIZE));
        alignment = std::max(alignment, size_type(CCTK_REAL_VEC_SIZE));
      }

      if (VECTORISE_ALIGN_FOR_CACHE) {
        for (size_t cache_level = 0; cache_level < cache_info.size();
             ++cache_level) {
          if (cache_info[cache_level].type == 0) {
            // Pad array in this direction to a multiple of this cache
            // line size
            int const cache_linesize = cache_info[cache_level].linesize;

            assert(cache_linesize % sizeof(CCTK_REAL) == 0);
            int const linesize = cache_linesize / sizeof(CCTK_REAL);
            if (npoints * accumulated_npoints < linesize) {
              // The extent is less than one cache line long: Ensure
              // that the array size divides the cache line size evenly
              // by rounding to the next power of 2
              assert(is_power_of_2(linesize));
              npoints = next_power_of_2(npoints);
              alignment = std::max(alignment, npoints);
            } else {
              // The extent is at least one cache line long: round up to
              // the next full cache line
              if (VECTORISE_ALIGN_INTERIOR) {
                if (d == 0) {
                  size_type newbndsize = align_up(bndsize, size_type(linesize));
                  npoints += newbndsize - bndsize;
                  bndsize = newbndsize;
                }
              }
              size_type total_npoints = npoints * accumulated_npoints;
              total_npoints = align_up(total_npoints, size_type(linesize));
              assert(total_npoints % accumulated_npoints == 0);
              npoints = total_npoints / accumulated_npoints;
              alignment = std::max(alignment, size_type(linesize));
            }

#if 0
        // Avoid multiples of the cache stride
        int const cache_stride = cache_info[cache_level].stride;
        if (cache_stride > 0) {
          assert(cache_stride % sizeof(CCTK_REAL) == 0);
          int const stride = cache_stride / sizeof(CCTK_REAL);
          if (npoints * accumulated_npoints % stride == 0) {
            assert(stride > linesize);
            size_type total_npoints = npoints * accumulated_npoints;
            total_npoints += std::max(size_type(linesize), accumulated_npoints);
            assert(total_npoints % accumulated_npoints == 0);
            npoints = total_npoints / accumulated_npoints;
          }
        }
#endif
          } // if is cache
        }   // for cache_level
      }
    }

    padded_shape[d] = npoints;
    padding_offset[d] = bndsize - offset[d];
    padding_alignment[d] = alignment;
    accumulated_npoints *= npoints;
  }
  assert(prod(vect<size_type, D>(padded_shape)) == accumulated_npoints);

  // self-check
  for (int d = 0; d < D; ++d) {
    assert(padded_shape[d] >= shape[d]);
    assert(padding_offset[d] >= 0);
    assert(padding_alignment[d] > 0);
    assert(padding_offset[d] < padding_alignment[d]);
    assert(padding_offset[d] + shape[d] <= padded_shape[d]);
    if (VECTORISE and VECTORISE_ALIGNED_ARRAYS) {
      if (d == 0) {
        assert(padded_shape[d] % CCTK_REAL_VEC_SIZE == 0);
        if (VECTORISE_ALIGN_INTERIOR)
          assert((padding_offset[d] + offset[d]) % CCTK_REAL_VEC_SIZE == 0);
      }
      if (not VECTORISE_ALIGN_INTERIOR)
        assert(padding_offset[d] == 0);
    } else {
      assert(padded_shape[d] == shape[d]);
    }
  }

  // Safety check
  if (not(prod(vect<size_type, D>(padded_shape)) <=
          2 * prod(vect<size_type, D>(shape)) + 1000)) {
    cerr << "shape=" << shape << "   "
         << "prod(shape)=" << prod(vect<size_type, D>(shape)) << "\n"
         << "padded_shape=" << padded_shape << "   "
         << "prod(padded_shape)=" << prod(vect<size_type, D>(padded_shape))
         << "\n"
         << "offset=" << offset << "\n"
         << "padding_alignment=" << padding_alignment << "\n"
         << "padding_offset=" << padding_offset << "\n";
  }
  assert(prod(vect<size_type, D>(padded_shape)) <=
         2 * prod(vect<size_type, D>(shape)) + 1000);

  if (verbose) {
    ostringstream buf;
    buf << "padding shape " << shape << " to " << padded_shape
        << "; padding offset " << offset << " to " << (padding_offset + offset)
        << "; padding alignment " << padding_alignment;
    CCTK_INFO(buf.str().c_str());
  }

  return {padded_shape, padding_alignment, padding_offset};
}

template padding_t<int, 3> pad_shape(bbox<int, 3> const &extent,
                                     bbox<int, 3> const &owned);
template padding_t<int, 3> pad_shape(vect<int, 3> const &shape,
                                     vect<int, 3> const &offset);

template padding_t<int, 4> pad_shape(bbox<int, 4> const &extent,
                                     bbox<int, 4> const &owned);
template padding_t<int, 4> pad_shape(vect<int, 4> const &shape,
                                     vect<int, 4> const &offset);
}
