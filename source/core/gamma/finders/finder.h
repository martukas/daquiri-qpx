#pragma once

#include <core/gamma/fit_settings.h>
//#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

struct DetectedPeak
{
  double center;
  double left;
  double right;
  double highest_y{0};
};

class Finder
{
 public:
  Finder(const std::vector<double>& x,
         const std::vector<double>& y)
      : x_(x), y_(y) {}
  virtual ~Finder() = default;

  DetectedPeak tallest_detected() const;
  double highest_residual(double l, double r) const;
  virtual double find_left(double chan) const = 0;
  virtual double find_right(double chan) const = 0;

  //DATA
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<DetectedPeak> filtered;

 protected:
  int32_t find_index(double chan_val) const;
};

}
