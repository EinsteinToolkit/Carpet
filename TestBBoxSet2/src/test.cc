#include <cctk_Parameters.h>

#include <cycleclock.h>

#include <bboxset.hh>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

namespace {

const double nano = 1.0e-9;

// Time a region of code
template <typename F>
auto xtime(double &time, const F &f) ->
    typename enable_if<is_void<decltype(f())>::value, void>::type {
  const ticks tb = getticks();
  f();
  const ticks te = getticks();
  time += seconds_per_tick() * elapsed(te, tb);
}

template <typename F>
auto xtime(double &time, const F &f) ->
    typename enable_if<not is_void<decltype(f())>::value, decltype(f())>::type {
  const ticks tb = getticks();
  const decltype(f()) ret = f();
  const ticks te = getticks();
  time += seconds_per_tick() * elapsed(te, tb);
  return ret;
}

// Convert to string
// (Would like to use to_string instead, but this doesn't seem to
// compile on Stampede)
template <typename T> string to_str(const T &x) {
  ostringstream buf;
  buf << x;
  return buf.str();
}

template <typename A1, typename A2, typename T, int D>
void check_equal(const string descr, const A1 &x1, const A2 &x2,
                 const bboxset1::bboxset<T, D> &bset1,
                 const bboxset2::bboxset<T, D> &bset2) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    cout << "  <" << D << "> " << descr << "...";
    cout.flush();
  }

  set<bbox<T, D> > v1, v2;
  bset1.serialise(v1);
  bset2.serialise(v2);

  bool eq = v1 == v2;
  if (not eq) {
    eq = (bset1 == bboxset1::bboxset<T, D>(v2) and
          bset2 == bboxset2::bboxset<T, D>(v1));
  }

  if (eq) {
    if (verbose)
      cout << " [success]\n";
  } else {
    if (verbose)
      cout << " [failure]\n";
    cout.flush();
    cerr.flush();
    fflush(stdout);
    fflush(stderr);
    CCTK_VParamWarn(CCTK_THORNSTRING,
                    "Found difference in operation \"%s\":", descr.c_str());
    fflush(stdout);
    fflush(stderr);
    cout << "\n"
         << "input1=" << x1 << "\n"
         << "input2=" << x2 << "\n"
         << "bset1=" << bset1 << "\n"
         << "bset2=" << bset2 << "\n"
         << "bset1.serialised=" << v1 << "\n"
         << "bset2.serialised=" << v2 << "\n";
    cout.flush();
    cerr.flush();

    const bboxset1::bboxset<T, D> bs(v2);
    cout << "bset1-bset2=" << bset1 - bs << "\n"
         << "bset2-bset1=" << bs - bset1 << "\n";
  }
}

template <typename T, int D>
void check_equal(const string descr, const bboxset1::bboxset<T, D> &bset1,
                 const bboxset2::bboxset<T, D> &bset2) {
  check_equal(descr, "", "", bset1, bset2);
}

#define bvect bvect_
#define b2vect b2vect_
#define ivect ivect_
#define i2vect i2vect_
#define ibbox ibbox_
template <int D> struct full_boxes_t {
  typedef ::vect<bool, D> bvect;
  typedef ::vect<bvect, 2> b2vect;
  typedef ::vect<int, D> ivect;
  typedef ::vect<ivect, 2> i2vect;
  typedef ::bbox<int, D> ibbox;
  typedef bboxset1::bboxset<int, D> ibset1;
  typedef bboxset2::bboxset<int, D> ibset2;
  ibbox interior;
  b2vect is_outer_boundary;
  ibbox exterior;
  ibset1 ghosts1;
  ibset2 ghosts2;
  ibbox communicated;
  ibset1 outer_boundaries1;
  ibset2 outer_boundaries2;
  ibbox owned;
  ibset1 boundaries1;
  ibset2 boundaries2;
  ibset1 buffers1;
  ibset2 buffers2;
  ibset1 overlaps1;
  ibset2 overlaps2;
  ibset1 active1;
  ibset2 active2;
};
template <int D> struct level_boxes_t {
  typedef ::vect<bool, D> bvect;
  typedef ::vect<bvect, 2> b2vect;
  typedef ::vect<int, D> ivect;
  typedef ::vect<ivect, 2> i2vect;
  typedef ::bbox<int, D> ibbox;
  typedef bboxset1::bboxset<int, D> ibset1;
  typedef bboxset2::bboxset<int, D> ibset2;
  ibset1 buffers1;
  ibset2 buffers2;
  ibset1 active1;
  ibset2 active2;
};

template <int D> void test() {
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Testing bboxset2<%d>...", D);

  // Prevent warnings about shadowing declarations
  typedef ::vect<bool, D> bvect;
  typedef ::vect<bvect, 2> b2vect;
  typedef ::vect<int, D> ivect;
  typedef ::vect<ivect, 2> i2vect;
  typedef ::bbox<int, D> ibbox;
  typedef bboxset1::bboxset<int, D> ibset1;
  typedef bboxset2::bboxset<int, D> ibset2;

  const int str = 3; // a non-trivial stride

  double time1, time2, time1a, time2a, time_check;

  // Some warm-up gymnastics with emtpy sets

  ibset1 bs1;
  ibset2 bs2;
  xtime(time_check, [&]() { check_equal("create ibset", bs1, bs2); });
  assert(bs1.empty());
  assert(bs2.empty());

  ibbox b;
  bs1 |= b;
  bs2 |= b;

  xtime(time_check, [&]() { check_equal("union with empty ibbox", bs1, bs2); });
  assert(bs1.empty());
  assert(bs2.empty());

  bs1 &= b;
  bs2 &= b;
  xtime(time_check,
        [&]() { check_equal("intersection with empty ibbox", bs1, bs2); });
  assert(bs1.empty());
  assert(bs2.empty());

  bs1 += b;
  bs2 += b;
  xtime(time_check,
        [&]() { check_equal("symmetric union with empty ibbox", bs1, bs2); });
  assert(bs1.empty());
  assert(bs2.empty());

  bs1 -= b;
  bs2 -= b;
  xtime(time_check,
        [&]() { check_equal("difference with empty ibbox", bs1, bs2); });
  assert(bs1.empty());
  assert(bs2.empty());

  // Some basic tests

  b = ibbox(ivect(0 * str), ivect(9 * str), ivect(str));
  bs1 += b;
  bs2 += b;
  xtime(time_check, [&]() { check_equal("add ibbox", bs1, bs2); });

  bs1 -= b;
  bs2 -= b;
  xtime(time_check, [&]() { check_equal("subtract ibbox", bs1, bs2); });

  bs1 |= b;
  bs2 |= b;
  xtime(time_check, [&]() { check_equal("union", bs1, bs2); });

  bs1 &= b;
  bs2 &= b;
  xtime(time_check, [&]() { check_equal("intersection", bs1, bs2); });

  ibset1 bs1a = bs1.shift(ivect(1));
  ibset2 bs2a = bs2.shift(ivect(1));
  xtime(time_check, [&]() { check_equal("shift", bs1a, bs2a); });

  ibset1 bs1b = bs1.expand(ivect(1), ivect(1));
  ibset2 bs2b = bs2.expand(ivect(1), ivect(1));
  xtime(time_check, [&]() { check_equal("expand", bs1b, bs2b); });

  ibset1 bs1c = bs1 & bs1a;
  ibset2 bs2c = bs2 & bs2a;
  xtime(time_check,
        [&]() { check_equal("non-trivial intersection", bs1c, bs2c); });

  ibset1 bs1d = bs1b - bs1;
  ibset2 bs2d = bs2b - bs2;
  xtime(time_check,
        [&]() { check_equal("non-trivial difference", bs1d, bs2d); });

  // Many tests with random sets

  bs1 = ibset1();
  bs2 = ibset2();
  xtime(time_check,
        [&]() { check_equal("many dense bboxes (init)", bs1, bs2); });

  time1 = time2 = time_check = 0.0;
  for (int n = 0; n < 40; ++n) {
    ivect lo, up;
    for (int d = 0; d < D; ++d)
      lo[d] = random() % 10 * str;
    for (int d = 0; d < D; ++d)
      up[d] = lo[d] + random() % 10 * str;
    b = ibbox(lo, up, ivect(str));
    xtime(time1, [&]() { bs1 |= b; });
    xtime(time2, [&]() { bs2 |= b; });
    xtime(time_check, [&]() {
      check_equal("many dense bboxes (union) #" + to_str(n), bs1, bs2);
    });
  }
  cout << "  "
       << "time1: " << time1 << " s, "
       << "time2: " << time2 << " s, "
       << "time_check: " << time_check << " s\n";

  time1 = time2 = time_check = 0.0;
  for (int n = 0; n < 100; ++n) {
    ivect lo, up;
    for (int d = 0; d < D; ++d)
      lo[d] = random() % 10 * str;
    for (int d = 0; d < D; ++d)
      up[d] = lo[d] + random() % 10 * str;
    b = ibbox(lo, up, ivect(str));
    xtime(time1, [&]() { bs1 -= b; });
    xtime(time2, [&]() { bs2 -= b; });
    xtime(time_check, [&]() {
      check_equal("many dense bboxes (difference) #" + to_str(n), bs1, bs2);
    });
  }
  cout << "  "
       << "time1: " << time1 << " s, "
       << "time2: " << time2 << " s, "
       << "time_check: " << time_check << " s\n";

  bs1 = ibset1();
  bs2 = ibset2();
  xtime(time_check,
        [&]() { check_equal("many sparse bboxes (init)", bs1, bs2); });

  time1 = time2 = time_check = 0.0;
  for (int n = 0; n < 10; ++n) {
    ivect lo, up;
    for (int d = 0; d < D; ++d)
      lo[d] = random() % 1000 * str;
    for (int d = 0; d < D; ++d)
      up[d] = lo[d] + random() % 1000 * str;
    b = ibbox(lo, up, ivect(str));
    xtime(time1, [&]() { bs1 |= b; });
    xtime(time2, [&]() { bs2 |= b; });
    xtime(time_check, [&]() {
      check_equal("many sparse bboxes (union) #" + to_str(n), bs1, bs2);
    });
  }
  cout << "  "
       << "time1: " << time1 << " s, "
       << "time2: " << time2 << " s, "
       << "time_check: " << time_check << " s\n";

  time1 = time2 = time_check = 0.0;
  for (int n = 0; n < 20; ++n) {
    ivect lo, up;
    for (int d = 0; d < D; ++d)
      lo[d] = random() % 1000 * str;
    for (int d = 0; d < D; ++d)
      up[d] = lo[d] + random() % 1000 * str;
    b = ibbox(lo, up, ivect(str));
    xtime(time1, [&]() { bs1 -= b; });
    xtime(time2, [&]() { bs2 -= b; });
    xtime(time_check, [&]() {
      check_equal("many sparse bboxes (difference) #" + to_str(n), bs1, bs2);
    });
  }
  cout << "  "
       << "time1: " << time1 << " s, "
       << "time2: " << time2 << " s, "
       << "time_check: " << time_check << " s\n";

  // Test shifting, expanding, contracting

  int sec = 0;
  for (int ts = 1; ts <= 3; ++ts) {
    const ivect target_stride(1 << ts);
    for (int to = 0; to < target_stride[0]; ++to) {
      const ivect target_offset(to);
      const ibbox target(target_offset, target_offset, target_stride);

      for (int s = 1; s <= 3; ++s) {
        const ivect stride(1 << s);
        for (int o = 0; o < stride[0]; ++o) {
          const ivect offset(o);

          bs1 = ibset1();
          bs2 = ibset2();
          for (int n = 0; n < 10; ++n) {
            ivect lo, up;
            for (int d = 0; d < D; ++d) {
              lo[d] = offset[d] + random() % 10 * stride[d];
              up[d] = lo[d] + random() % 10 * stride[d];
            }
            b = ibbox(lo, up, stride);
            bs1 |= b;
            bs2 |= b;
            xtime(time_check, [&]() {
              check_equal("repeated test #" + to_str(sec) + "." + to_str(n),
                          bs1, bs2);
            });

            bs1a = bs1.expanded_for(target);
            bs2a = bs2.expanded_for(target);
            xtime(time_check, [&]() {
              check_equal("expanded_for #" + to_str(sec) + "." + to_str(n), bs1,
                          target, bs1a, bs2a);
            });
            bs1b = bs1.contracted_for(target);
            bs2b = bs2.contracted_for(target);
            // Omitting check since bboxset1 has known bugs here
            // xtime(time_check, [&]() { check_equal("contracted_for", bs1,
            // target, bs1b, bs2b); });

          } // n
          ++sec;

        } // offset
      }   // stride

    } // target_offset
  }   // target_stride

  // Some real-world tests

  {

    // const auto minus1 =
    //   (ibset1(*)(const ibbox&, const ibbox&))bboxset1::operator-;
    // const auto minus2 =
    //   (ibset2(*)(const ibbox&, const ibbox&))bboxset2::operator-;

    // ibset1(*const minus1)(const ibbox&, const ibbox&)(bboxset1::operator-);
    // ibset2(*const minus2)(const ibbox&, const ibbox&)(bboxset2::operator-);

    typedef ibset1(minus1_t)(const ibbox &, const ibbox &);
    typedef ibset2(minus2_t)(const ibbox &, const ibbox &);
    // typedef ibset2 (&xor2_t)(const ibbox&, const ibbox&);
    minus1_t &minus1(bboxset1::operator-);
    minus2_t &minus2(bboxset2::operator-);
    // const xor2_t xor2(bboxset2::operator^);

    const i2vect ghost_width(3);
    const i2vect buffer_width(9);
    const i2vect overlap_width(0);
    const i2vect boundary_width(3);

    const int components_per_direction = lrint(pow(components_goal, 1.0 / D));
    const int components = lrint(pow(components_per_direction, D));
    const ivect origin(0);
    const ivect stride(1);
    const ivect component_size(30);

    cout << "  components: " << components << "\n";

    // Domain

    time1 = time2 = time1a = time2a = time_check = 0.0;

    const ibbox domain_exterior(
        origin,
        origin + stride * component_size * ivect(components_per_direction) -
            stride,
        stride);

    const ibbox domain_active = domain_exterior.expand(-boundary_width);
    assert(domain_active <= domain_exterior);
    ibset1 domain_boundary1;
    ibset2 domain_boundary2;
    xtime(time1,
          [&]() { domain_boundary1 = minus1(domain_exterior, domain_active); });
    xtime(time2,
          [&]() { domain_boundary2 = minus2(domain_exterior, domain_active); });
    xtime(time_check, [&]() {
      check_equal("domain_boundary", domain_boundary1, domain_boundary2);
    });

    cout << "  Domain:\n"
         << "    time1:        " << time1 << " s\n"
         << "    time2:        " << time2 << " s\n"
         << "    time1_assert: " << time1a << " s\n"
         << "    time2_assert: " << time2a << " s\n"
         << "    time_check:   " << time_check << " s\n";

    // Region

    time1 = time2 = time1a = time2a = time_check = 0.0;

    vector<full_boxes_t<D> > full_boxes(components);
    level_boxes_t<D> level_boxes;
    for (int c = 0; c < components; ++c) {

      // Interior:
      ibbox &intr = full_boxes.AT(c).interior;
      ivect ipos;
      int ic = c;
      for (int d = 0; d < D; ++d) {
        ipos[d] = ic % components_per_direction;
        ic /= components_per_direction;
      }
      assert(ic == 0);
      intr =
          ibbox(origin + stride * component_size * ipos,
                origin + stride * component_size * (ipos + 1) - stride, stride);
      assert(intr <= domain_exterior);
      for (int cc = 0; cc < c; ++cc) {
        assert(not intr.intersects(full_boxes.AT(cc).interior));
      }

      // Outer boundary faces:
      b2vect &is_outer_boundary = full_boxes.AT(c).is_outer_boundary;
      is_outer_boundary[0] = intr.lower() == domain_exterior.lower();
      is_outer_boundary[1] = intr.upper() == domain_exterior.upper();

      // Exterior:
      ibbox &extr = full_boxes.AT(c).exterior;
      extr = intr.expand(i2vect(not is_outer_boundary) * ghost_width);
      assert(extr <= domain_exterior);

      // Cactus ghost zones (which include outer boundaries):
      ibset1 &ghosts1 = full_boxes.AT(c).ghosts1;
      ibset2 &ghosts2 = full_boxes.AT(c).ghosts2;
      xtime(time1, [&]() { ghosts1 = minus1(extr, intr); });
      xtime(time2, [&]() { ghosts2 = minus2(extr, intr); });
      xtime(time_check, [&]() { check_equal("ghosts", ghosts1, ghosts2); });
      xtime(time1a, [&]() { assert(ghosts1 <= domain_exterior); });
      xtime(time2a, [&]() { assert(ghosts2 <= domain_exterior); });

      // Communicated region:
      ibbox &comm = full_boxes.AT(c).communicated;
      comm = extr.expand(i2vect(is_outer_boundary) * (-boundary_width));
      assert(comm <= domain_active);

      // Outer boundary:
      ibset1 &outer_boundaries1 = full_boxes.AT(c).outer_boundaries1;
      ibset2 &outer_boundaries2 = full_boxes.AT(c).outer_boundaries2;
      xtime(time1, [&]() { outer_boundaries1 = minus1(extr, comm); });
      xtime(time2, [&]() { outer_boundaries2 = minus2(extr, comm); });
      xtime(time_check, [&]() {
        check_equal("outer_boundaries", outer_boundaries1, outer_boundaries2);
      });
      xtime(time1a, [&]() { assert(outer_boundaries1 <= domain_boundary1); });
      xtime(time2a, [&]() { assert(outer_boundaries2 <= domain_boundary2); });

      // Owned region:
      ibbox &owned = full_boxes.AT(c).owned;
      owned = intr.expand(i2vect(is_outer_boundary) * (-boundary_width));
      assert(owned <= domain_active);
      for (int cc = 0; cc < c; ++cc) {
        assert(not owned.intersects(full_boxes.AT(cc).owned));
      }

      // Boundary (Carpet ghost zones, which do not include outer
      // boundaries):
      ibset1 &boundaries1 = full_boxes.AT(c).boundaries1;
      ibset2 &boundaries2 = full_boxes.AT(c).boundaries2;
      xtime(time1, [&]() { boundaries1 = minus1(comm, owned); });
      xtime(time2, [&]() { boundaries2 = minus2(comm, owned); });
      xtime(time_check,
            [&]() { check_equal("boundaries", boundaries1, boundaries2); });
      xtime(time1a, [&]() { assert(boundaries1 <= domain_active); });
      xtime(time2a, [&]() { assert(boundaries2 <= domain_active); });

    } // for c

    cout << "  Region:\n"
         << "    time1:        " << time1 << " s\n"
         << "    time2:        " << time2 << " s\n"
         << "    time1_assert: " << time1a << " s\n"
         << "    time2_assert: " << time2a << " s\n"
         << "    time_check:   " << time_check << " s\n";

    // Buffer zones

    time1 = time2 = time1a = time2a = time_check = 0.0;

    // All owned regions
    ibset1 allowned1;
    ibset2 allowned2;
    xtime(time1,
          [&]() { allowned1 = ibset1(full_boxes, &full_boxes_t<D>::owned); });
    xtime(time2,
          [&]() { allowned2 = ibset2(full_boxes, &full_boxes_t<D>::owned); });
    xtime(time_check, [&]() { check_equal("allowned", allowned1, allowned2); });
    xtime(time1a, [&]() { assert(allowned1 <= domain_active); });
    xtime(time2a, [&]() { assert(allowned2 <= domain_active); });

    // All not-owned regions
    ibset1 notowned1;
    ibset2 notowned2;
    xtime(time1, [&]() { notowned1 = domain_active - allowned1; });
    xtime(time2, [&]() { notowned2 = domain_active - allowned2; });
    xtime(time_check, [&]() { check_equal("notowned", notowned1, notowned2); });

    // All not-active points
    ibset1 notactive1;
    ibset2 notactive2;
    xtime(time1, [&]() {
      notactive1 = notowned1.expand(buffer_width + overlap_width);
    });
    xtime(time2, [&]() {
      notactive2 = notowned2.expand(buffer_width + overlap_width);
    });
    xtime(time_check,
          [&]() { check_equal("notactive", notactive1, notactive2); });

    // All not-active points, in stages
    int const num_substeps = any(any(ghost_width == 0))
                                 ? 0
                                 : minval(minval(buffer_width / ghost_width));
    assert(all(all(buffer_width == num_substeps * ghost_width)));
    vector<ibset1> notactive_stepped1(num_substeps + 1);
    vector<ibset2> notactive_stepped2(num_substeps + 1);
    notactive_stepped1.AT(0) = notowned1;
    notactive_stepped2.AT(0) = notowned2;
    for (int substep = 1; substep <= num_substeps; ++substep) {
      xtime(time1, [&]() {
        notactive_stepped1.AT(substep) =
            notactive_stepped1.AT(substep - 1).expand(ghost_width);
      });
      xtime(time2, [&]() {
        notactive_stepped2.AT(substep) =
            notactive_stepped2.AT(substep - 1).expand(ghost_width);
      });
      xtime(time_check, [&]() {
        check_equal("notactive_stepped[" + to_str(substep) + "]",
                    notactive_stepped1.AT(substep),
                    notactive_stepped2.AT(substep));
      });
    }
    ibset1 notactive_overlaps1;
    ibset2 notactive_overlaps2;
    xtime(time1a, [&]() {
      notactive_overlaps1 =
          notactive_stepped1.AT(num_substeps).expand(overlap_width);
    });
    xtime(time2a, [&]() {
      notactive_overlaps2 =
          notactive_stepped2.AT(num_substeps).expand(overlap_width);
    });
    xtime(time_check, [&]() {
      check_equal("notactive_overlaps", notactive_overlaps1,
                  notactive_overlaps2);
    });
    if (all(all(buffer_width == num_substeps * ghost_width))) {
      xtime(time1a, [&]() { assert(notactive_overlaps1 == notactive1); });
      xtime(time2a, [&]() { assert(notactive_overlaps2 == notactive2); });
    }

    // All buffer zones
    ibset1 &allbuffers1 = level_boxes.buffers1;
    ibset2 &allbuffers2 = level_boxes.buffers2;
    xtime(time1,
          [&]() { allbuffers1 = allowned1 & notowned1.expand(buffer_width); });
    xtime(time2,
          [&]() { allbuffers2 = allowned2 & notowned2.expand(buffer_width); });
    xtime(time_check,
          [&]() { check_equal("allbuffers", allbuffers1, allbuffers2); });

    // All overlap zones
    ibset1 alloverlaps1;
    ibset2 alloverlaps2;
    xtime(time1,
          [&]() { alloverlaps1 = (allowned1 & notactive1) - allbuffers1; });
    xtime(time2,
          [&]() { alloverlaps2 = (allowned2 & notactive2) - allbuffers2; });
    xtime(time_check,
          [&]() { check_equal("alloverlaps", alloverlaps1, alloverlaps2); });

    // All active points
    ibset1 &allactive1 = level_boxes.active1;
    ibset2 &allactive2 = level_boxes.active2;
    xtime(time1, [&]() { allactive1 = allowned1 - notactive1; });
    xtime(time2, [&]() { allactive2 = allowned2 - notactive2; });
    xtime(time_check,
          [&]() { check_equal("allactive", allactive1, allactive2); });

    // All stepped buffer zones
    vector<ibset1> allbuffers_stepped1(num_substeps);
    vector<ibset2> allbuffers_stepped2(num_substeps);
    ibset1 allbuffers_stepped_combined1;
    ibset2 allbuffers_stepped_combined2;
    for (int substep = 0; substep < num_substeps; ++substep) {
      xtime(time1, [&]() {
        allbuffers_stepped1.AT(substep) =
            allowned1 & (notactive_stepped1.AT(substep + 1) -
                         notactive_stepped1.AT(substep));
      });
      xtime(time2, [&]() {
        allbuffers_stepped2.AT(substep) =
            allowned2 & (notactive_stepped2.AT(substep + 1) -
                         notactive_stepped2.AT(substep));
      });
      xtime(time_check, [&]() {
        check_equal("notactive_stepped[" + to_str(substep) + "]",
                    notactive_stepped1.AT(substep),
                    notactive_stepped2.AT(substep));
      });
      xtime(time1a, [&]() {
        allbuffers_stepped_combined1 += allbuffers_stepped1.AT(substep);
      });
      xtime(time2a, [&]() {
        allbuffers_stepped_combined2 += allbuffers_stepped2.AT(substep);
      });
      xtime(time_check, [&]() {
        check_equal("allbuffers_stepped_combined", allbuffers_stepped_combined1,
                    allbuffers_stepped_combined2);
      });
    }
    if (all(all(buffer_width == num_substeps * ghost_width))) {
      xtime(time1a,
            [&]() { assert(allbuffers_stepped_combined1 == allbuffers1); });
      xtime(time2a,
            [&]() { assert(allbuffers_stepped_combined2 == allbuffers2); });
    }

    // Overlap zones and buffer zones must be in the active part of
    // the domain
    xtime(time1a, [&]() { assert(allactive1 <= domain_active); });
    xtime(time2a, [&]() { assert(allactive2 <= domain_active); });
    xtime(time1a, [&]() { assert(alloverlaps1 <= domain_active); });
    xtime(time2a, [&]() { assert(alloverlaps2 <= domain_active); });
    xtime(time1a, [&]() { assert(allbuffers1 <= domain_active); });
    xtime(time2a, [&]() { assert(allbuffers2 <= domain_active); });
    xtime(time1a, [&]() { assert((allactive1 & alloverlaps1).empty()); });
    xtime(time2a, [&]() { assert((allactive2 & alloverlaps2).empty()); });
    xtime(time1a, [&]() { assert((allactive1 & allbuffers1).empty()); });
    xtime(time2a, [&]() { assert((allactive2 & allbuffers2).empty()); });
    xtime(time1a, [&]() { assert((alloverlaps1 & allbuffers1).empty()); });
    xtime(time2a, [&]() { assert((alloverlaps2 & allbuffers2).empty()); });
    xtime(time1a, [&]() {
      assert(allactive1 + alloverlaps1 + allbuffers1 == allowned1);
    });
    xtime(time2a, [&]() {
      assert(allactive2 + alloverlaps2 + allbuffers2 == allowned2);
    });

    for (int c = 0; c < components; ++c) {
      auto &box = full_boxes.AT(c);

      // Buffer zones:
      xtime(time1, [&]() { box.buffers1 = box.owned & allbuffers1; });
      xtime(time2, [&]() { box.buffers2 = box.owned & allbuffers2; });
      xtime(time_check,
            [&]() { check_equal("buffers", box.buffers1, box.buffers2); });

      // Overlap zones:
      xtime(time1, [&]() { box.overlaps1 = box.owned & alloverlaps1; });
      xtime(time2, [&]() { box.overlaps2 = box.owned & alloverlaps2; });
      xtime(time_check,
            [&]() { check_equal("overlaps", box.overlaps1, box.overlaps2); });

      // Active region:
      xtime(time1, [&]() { box.active1 = box.owned & allactive1; });
      xtime(time2, [&]() { box.active2 = box.owned & allactive2; });
      xtime(time_check,
            [&]() { check_equal("active", box.active1, box.active2); });
      xtime(time1a, [&]() {
        assert(box.active1 == box.owned - (box.buffers1 + box.overlaps1));
      });
      xtime(time2a, [&]() {
        assert(box.active2 == box.owned - (box.buffers2 + box.overlaps2));
      });
    } // for c

#if 0
      for (int lc = 0; lc < h.local_components(rl); ++ lc) {
        int const c = h.get_component (rl, lc);
        local_dboxes & local_box = local_level.AT(lc);
        full_dboxes const& box = full_boxes.AT(c);
        
        local_box.buffers = box.buffers;
        local_box.overlaps = box.overlaps;
        
        local_box.active = box.active;
      } // for lc
      
      // The conjunction of all buffer zones must equal allbuffers
      ibset const allbuffers1 (full_boxes, & full_dboxes::buffers);
      ASSERT_rl (allbuffers1 == allbuffers,
                 "Buffer zone consistency check");
#endif

    cout << "  Buffer zones:\n"
         << "    time1:        " << time1 << " s\n"
         << "    time2:        " << time2 << " s\n"
         << "    time1_assert: " << time1a << " s\n"
         << "    time2_assert: " << time2a << " s\n"
         << "    time_check:   " << time_check << " s\n";
  }

  CCTK_VInfo(CCTK_THORNSTRING, "Done testing bboxset2<%d>.", D);
}

} // namespace

extern "C" void TestBBoxSet2_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  test<1>();
  test<2>();
  test<3>();
}
