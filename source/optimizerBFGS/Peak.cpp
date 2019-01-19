#include <optimizerBFGS/Peak.h>

#include <core/util/custom_logger.h>
#include <optimizerBFGS/more_math.h>

namespace Hypermet
{

double ValueDefault::max() const
{
  return max_;
}

void ValueDefault::max(double new_max)
{
  //Dim val As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(x_))
  //Value = std::min(NewMax, val)
  max_ = new_max;
  //x_ = std::Asin((_Min + _Max - 2 * val) / (-1 * _Max + _Min))
}

double ValueDefault::min()
{
  return min_;
}

void ValueDefault::min(double new_min)
{
  //Dim val As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(x_))
  min_ = new_min;
  //temp = std::max(temp, val)
  //x_ = std::Asin((_Min + _Max - 2 * temp) / (-1 * _Max + _Min))
}

double ValueDefault::x()
{
  return x_;
}

void ValueDefault::x(double new_x)
{
  x_ = new_x;
}

double ValueDefault::val() const
{
  return val_at(x_);
}

void ValueDefault::val(double new_val)
{
  double t = (min_ + max_ - 2.0 * new_val) / (min_ - max_);
  if (std::abs(t) <= 1)
    x_ = std::asin((min_ + max_ - 2.0 * new_val) / (min_ - max_));
  else if (signum(t) < 0)
    x_ = std::asin(-1);
  else
    x_ = std::asin(1);
}

double ValueDefault::val_at(double at_x) const
{
  return min_ + (max_ - min_) / 2.0 * (1.0 + std::sin(at_x));
}

double ValueDefault::grad_at(double at_x) const
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

int32_t Peak::step_type() const
{
  return FEP_status_;
}

double Peak::peak_position() const
{
  return position.val();
}

double Peak::peak_position_unc() const
{
  return position.uncert_value;
}

double Peak::peak_energy(const Calibration& cal) const
{
  return cal.channel_to_energy(position.val());
}

double Peak::peak_energy_unc(const Calibration& cal) const
{
  return cal.energy_slope() * position.uncert_value;
}

bool Peak::full_energy_peak() const
{
  if (FEP_status_ == 1)
    return true;
  else if (FEP_status_ == -1)
    return false;
  return false; // false? Ask Laci
}

void Peak::full_energy_peak(bool flag)
{
  if (flag)
    FEP_status_ = 1;
  else
    FEP_status_ = -1;
}

bool Peak::operator<(const Peak& other) const
{
  return position.val() < other.position.val();
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
