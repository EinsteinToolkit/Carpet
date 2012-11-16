#ifndef CACHEINFO_HH
#define CACHEINFO_HH

#include <cctk.h>

#include <limits>

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"



static ptrdiff_t div_down (ptrdiff_t const x, ptrdiff_t const align)
  CCTK_ATTRIBUTE_UNUSED;
static ptrdiff_t div_down (ptrdiff_t const x, ptrdiff_t const align)
{
  assert (x >= 0);
  assert (align > 0);
  return x / align;
}

static ptrdiff_t div_up (ptrdiff_t const x, ptrdiff_t const align)
  CCTK_ATTRIBUTE_UNUSED;
static ptrdiff_t div_up (ptrdiff_t const x, ptrdiff_t const align)
{
  assert (x >= 0);
  assert (align > 0);
  return (x + align - 1) / align;
}

static ptrdiff_t align_down (ptrdiff_t const x, ptrdiff_t const align)
  CCTK_ATTRIBUTE_UNUSED;
static ptrdiff_t align_down (ptrdiff_t const x, ptrdiff_t const align)
{
  assert (x >= 0);
  assert (align > 0);
  return div_down(x, align) * align;
}

static ptrdiff_t align_up (ptrdiff_t const x, ptrdiff_t const align)
  CCTK_ATTRIBUTE_UNUSED;
static ptrdiff_t align_up (ptrdiff_t const x, ptrdiff_t const align)
{
  assert (x >= 0);
  assert (align > 0);
  return div_up(x, align) * align;
}

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



class cacheinfo_t {
  ptrdiff_t m_size;                // bytes
  ptrdiff_t m_linesize;            // bytes (pagesize for TLBs)
  ptrdiff_t m_associativity;
public:
  bool invariant () const
  {
    return
      is_power_of_2(m_size) and
      is_power_of_2(m_linesize) and
      is_power_of_2(m_associativity) and
      m_size % (m_linesize * m_associativity) == 0;
  }
  cacheinfo_t (ptrdiff_t const a_size,
               ptrdiff_t const a_linesize,
               ptrdiff_t const a_associativity)
    : m_size (a_size),
      m_linesize (a_linesize),
      m_associativity (a_associativity)
  {
    assert (invariant());
  }
  cacheinfo_t (ptrdiff_t const a_linesize)
    : m_size (previous_power_of_2(numeric_limits<ptrdiff_t>::max())),
      m_linesize (a_linesize),
      m_associativity (1)
  {
    assert (invariant());
  }
  // size in bytes
  ptrdiff_t size() const
  {
    return m_size;
  }
  // line size in bytes
  ptrdiff_t linesize() const
  {
    return m_linesize;
  }
  // associativity
  ptrdiff_t associativity() const
  {
    return m_associativity;
  }
  // number of cache elements
  ptrdiff_t num_elements() const
  {
    return size() / (linesize() * associativity());
  }
  // stride (between main memory locations that use the same cache
  // element) in bytes
  ptrdiff_t stride() const
  {
    return num_elements() * linesize();
  }
};



// These routines are apparently not pure -- don't know why
template<int D>
vect<int,D> pad_shape (bbox<int,D> const& extent) /*CCTK_ATTRIBUTE_PURE*/;
template<int D>
vect<int,D> pad_shape (vect<int,D> const& shape) /*CCTK_ATTRIBUTE_PURE*/;

#endif  // CACHEINFO_HH
