#pragma once
#include <math.h>

template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
  return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
  return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
  return signum(x, std::is_signed<T>());
}

template <typename T> inline constexpr
T square(T x) {
  return std::pow(x, 2);
}

template <typename T> inline constexpr
T cube(T x) {
  return std::pow(x, 3);
}

inline double erfc(double x) {
  static constexpr double p{0.47047};
  static constexpr double b1{0.1740121};
  static constexpr double b2{-0.0479399};
  static constexpr double b3{0.3739278};

  try
  {
    double x0, cutf, T;
    if (x < 0)
      x0 = -x;
    else
      x0 = x;
    T = 1.0 / (1.0 + p * x0);
    cutf = (b1 * T + b2 * T * T + b3 * T * T * T)
        * std::exp(-(x0 * x0));
    if (x < 0)
      return 2 - 2 * cutf;
    else
      return 2 * cutf;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error(""));
  }
}
