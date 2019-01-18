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
