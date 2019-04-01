#include <core/fitting/parameter/sine_bounded_param.h>
#include <core/util/more_math.h>

//#include <mpreal.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void SineBoundedValue::val(double new_val)
{
  double t = (min() + max() - 2.0 * new_val) / (min() - max());
  if (std::abs(t) <= 1)
    x(std::asin((min() + max() - 2.0 * new_val) / (min() - max())));
  else if (signum(t) < 0)
    x(std::asin(-1));
  else
    x(std::asin(1));
}

double SineBoundedValue::val_at(double at_x) const
{
//  mpfr::mpreal mx = at_x;
//  mpfr::mpreal one = 1.0;
//  mpfr::mpreal two = 2.0;
//  mpfr::mpreal mmax = max_;
//  mpfr::mpreal mmin = min_;
//  auto ret = (one + mpfr::sin(mx)) * (mmax - mmin) / two + mmin;
//  return ret.toDouble();
  return (1.0 + std::sin(at_x)) * (max() - min()) / 2.0 + min();
}

double SineBoundedValue::grad_at(double at_x) const
{
//  mpfr::mpreal mx = at_x;
//  mpfr::mpreal two = 2.0;
//  mpfr::mpreal mmax = max_;
//  mpfr::mpreal mmin = min_;
//  auto ret = mpfr::cos(mx) * (mmax - mmin) / two;
//  return ret.toDouble();
  return std::cos(at_x) * (max() - min()) / 2.0;
}


}
