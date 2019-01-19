#pragma once

#include <optimizerBFGS/Calibration.h>

namespace Hypermet
{

class ValueDefault
{
 private:
  double x_{0};
  double dx_{0};
  double max_{1};
  double min_{0};
 public:
  double uncert_value{0};
  int32_t x_index{-1};

  //void SetRange(double NewMin, double NewMax) {
  //    if(Not NewMin < NewMax) { return; }
  //    double val = _Min + (_Max - _Min) / 2 * (1 + std::sin(x_));
  //    Value = std::min(NewMin, val);
  //    Value = std::max(NewMax, val);
  //    _Min = NewMin;
  //    _Max = NewMax;
  //    x_ = std::asin((_Min + _Max - 2 * val) / (-1 * _Max + _Min));
  //}

  double max() const;
  void max(double new_max);
  double min();
  void min(double new_min);
  double x();
  void x(double new_x);
  double val() const;
  void val(double new_val);
  double val_at(double at_x) const;
  double grad_at(double at_x) const;
};

class Value : public ValueDefault
{
 public:
  bool to_fit{true};
};

class ValueGam
{
 private:
  double x_{0};
  double dx_{0};
 public:
  double uncert_value{0};
  int32_t x_index{-1};

 public:
  double x();
  void x(double new_x);
  double val() const;
  void val(double new_val);
  double val_at(double at_x) const;
  double grad_at(double at_x) const;
};

class Peak
{
 private:
  int32_t FEP_status_{1};
 public:
  ValueDefault position;
  ValueGam GAM;

  int32_t step_type() const;
  double peak_position() const;
  double peak_position_unc() const;
  double peak_energy(const Calibration& cal) const;
  double peak_energy_unc(const Calibration& cal) const;
  bool full_energy_peak() const;
  void full_energy_peak(bool flag);
  bool operator<(const Peak& other) const;
};

class ValueBkgDefault
{
 private:
  double x_{0};
  double dx_{0};

 public:
  double uncert_value{0};
  double x_index{-1};

  double x() const;
  void x(double new_x);
  double val() const;
  void val(double new_val);
  double val_at(double at_x) const;
  double grad_at(double at_x) const;
};

class ValueBkg
{
 private:
  double x_{0};
  double dx_{0};

 public:
  double uncert_value{0};
  int32_t x_index{-1};
  bool to_fit{true};

  double x() const;
  void x(double new_x);
  double val() const;
  void val(double new_val);
  double val_at(double at_x) const;
  double grad_at(double at_x) const;
};

}
