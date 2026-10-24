// Standalone checks for the cubic Hermite dense output in hermite.hh.
// Build: g++ -std=c++11 -Wall -Wextra -Werror -o hermite_test hermite_test.cc

#include "hermite.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

int failures = 0;

void check(bool const ok, char const *const text, int const line) {
  if (ok)
    return;
  std::fprintf(stderr, "FAIL %s:%d: %s\n", "hermite_test.cc", line, text);
  ++failures;
}

#define CHECK(cond) check(static_cast<bool>(cond), #cond, __LINE__)

bool near(double const got, double const want) {
  return std::fabs(got - want) <= 1.0e-12 * (1.0 + std::fabs(want));
}

void check_endpoints() {
  double const y0 = 1.25;
  double const y1 = -0.5;
  double const dy0 = 0.75;
  double const dy1 = -2.0;
  double const dt = 0.125;
  using Carpet::hermite::slope;
  using Carpet::hermite::value;
  CHECK(near(value(0.0, dt, y0, y1, dy0, dy1), y0));
  CHECK(near(value(1.0, dt, y0, y1, dy0, dy1), y1));
  CHECK(near(slope(0.0, dt, y0, y1, dy0, dy1), dy0));
  CHECK(near(slope(1.0, dt, y0, y1, dy0, dy1), dy1));
}

// y = 1 + 2t + 3t^2 + 4t^3 on [0, 1/2] is reproduced exactly.
void check_cubic() {
  double const dt = 0.5;
  double const y0 = 1.0;
  double const y1 = 3.25;
  double const dy0 = 2.0;
  double const dy1 = 8.0;
  double const theta = 0.5;
  using Carpet::hermite::slope;
  using Carpet::hermite::value;
  CHECK(near(value(theta, dt, y0, y1, dy0, dy1), 1.75));
  CHECK(near(slope(theta, dt, y0, y1, dy0, dy1), 4.25));
  // RK4-2(2) evaluates a stage at c = -99/50. The cubic is not clamped.
  double const outside[2] = {-1.98, 1.01};
  for (int i = 0; i < 2; ++i) {
    double const t = outside[i] * dt;
    double const want = 1.0 + 2.0 * t + 3.0 * t * t + 4.0 * t * t * t;
    CHECK(near(value(outside[i], dt, y0, y1, dy0, dy1), want));
  }
}

void check_slope_levels() {
  using Carpet::hermite::slope_levels;
  CHECK(slope_levels(false).start == 1);
  CHECK(slope_levels(false).end == 0);
  CHECK(slope_levels(true).start == 0);
  CHECK(slope_levels(true).end == 2);
}

} // namespace

int main() {
  check_endpoints();
  check_cubic();
  check_slope_levels();
  if (failures != 0) {
    std::fprintf(stderr, "%d Hermite checks failed\n", failures);
    return EXIT_FAILURE;
  }
  std::printf("hermite_test: ok\n");
  return EXIT_SUCCESS;
}
