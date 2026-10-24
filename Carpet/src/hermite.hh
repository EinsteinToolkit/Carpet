#ifndef CARPET_HERMITE_HH
#define CARPET_HERMITE_HH

// Cubic Hermite dense output of one parent Runge-Kutta step.
//
// The unique degree-3 polynomial that matches y and y' at both ends of
// [t_n, t_n+dt]. theta is (t - t_n) / dt. The slope weights multiply dt*y',
// so a caller passes the derivatives themselves and the step size separately.
// k2 and k3 are not part of this polynomial.

namespace Carpet {
namespace hermite {

struct Weights {
  double value_start;
  double value_end;
  double slope_start;
  double slope_end;
};

inline Weights weights(double const theta) {
  double const theta2 = theta * theta;
  double const theta3 = theta2 * theta;
  Weights w;
  w.value_start = 1.0 - 3.0 * theta2 + 2.0 * theta3;
  w.value_end = 3.0 * theta2 - 2.0 * theta3;
  w.slope_start = theta - 2.0 * theta2 + theta3;
  w.slope_end = -theta2 + theta3;
  return w;
}

inline double value(double const theta, double const dt, double const y_start,
                    double const y_end, double const dy_start,
                    double const dy_end) {
  Weights const w = weights(theta);
  return w.value_start * y_start + w.value_end * y_end +
         dt * (w.slope_start * dy_start + w.slope_end * dy_end);
}

// dy/dt of the same polynomial. Used by the unit test to lock the slopes.
inline double slope(double const theta, double const dt, double const y_start,
                    double const y_end, double const dy_start,
                    double const dy_end) {
  double const theta2 = theta * theta;
  double const d_value_start = -6.0 * theta + 6.0 * theta2;
  double const d_value_end = 6.0 * theta - 6.0 * theta2;
  double const d_slope_start = 1.0 - 4.0 * theta + 3.0 * theta2;
  double const d_slope_end = -2.0 * theta + 3.0 * theta2;
  return (d_value_start * y_start + d_value_end * y_end) / dt +
         d_slope_start * dy_start + d_slope_end * dy_end;
}

} // namespace hermite
} // namespace Carpet

#endif
