#include <optimizerBFGS/Value.h>

#include <core/util/custom_logger.h>
#include <optimizerBFGS/more_math.h>

namespace Hypermet
{

double Value::max() const
{
  return max_;
}

void Value::max(double new_max)
{
  //Dim val As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(x_))
  //Value = std::min(NewMax, val)
  max_ = new_max;
  //x_ = std::Asin((_Min + _Max - 2 * val) / (-1 * _Max + _Min))
}

double Value::min()
{
  return min_;
}

void Value::min(double new_min)
{
  //Dim val As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(x_))
  min_ = new_min;
  //temp = std::max(temp, val)
  //x_ = std::Asin((_Min + _Max - 2 * temp) / (-1 * _Max + _Min))
}

void Value::bound(double v1, double v2)
{
  min(std::min(v1, v2));
  max(std::max(v1, v2));
}

double Value::x()
{
  return x_;
}

void Value::x(double new_x)
{
  x_ = new_x;
}

double Value::val() const
{
  return val_at(x_);
}

double Value::grad() const
{
  return grad_at(x_);
}

void Value::val(double new_val)
{
  double t = (min_ + max_ - 2.0 * new_val) / (min_ - max_);
  if (std::abs(t) <= 1)
    x_ = std::asin((min_ + max_ - 2.0 * new_val) / (min_ - max_));
  else if (signum(t) < 0)
    x_ = std::asin(-1);
  else
    x_ = std::asin(1);
}

double Value::val_at(double at_x) const
{
  return min_ + (max_ - min_) / 2.0 * (1.0 + std::sin(at_x));
}

double Value::grad_at(double at_x) const
{
  return (max_ - min_) * std::cos(at_x) / 2.0;
}

double ValueGam::x()
{
  return x_;
}

void ValueGam::x(double new_x)
{
  x_ = new_x;
}

double ValueGam::val() const
{
  return val_at(x_);
}

double ValueGam::grad() const
{
  return grad_at(x_);
}

void ValueGam::val(double new_val)
{
  x_ = std::sqrt(new_val);
}

double ValueGam::val_at(double at_x) const
{
  return square(at_x);
}

double ValueGam::grad_at(double at_x) const
{
  return 2.0 * at_x;
}

double ValueBkgDefault::x() const
{
  return x_;
}

void ValueBkgDefault::x(double new_x)
{
  x_ = new_x;
}

double ValueBkgDefault::val() const
{
  return val_at(x_);
}

void ValueBkgDefault::val(double new_val)
{
  x_ = std::sqrt(new_val);
}

double ValueBkgDefault::val_at(double at_x) const
{
  return square(at_x);
}

double ValueBkgDefault::grad_at(double at_x) const
{
  return 2.0 * at_x;
}

double ValueBkg::x() const
{
  return x_;
}

void ValueBkg::x(double new_x)
{
  x_ = new_x;
}

double ValueBkg::val() const
{
  return x_;
}

void ValueBkg::val(double new_val)
{
  x_ = new_val;
}

double ValueBkg::val_at(double at_x) const
{
  return at_x;
}

double ValueBkg::grad_at(double at_x) const
{
  (void) at_x;
  return 1.0;
}

}
