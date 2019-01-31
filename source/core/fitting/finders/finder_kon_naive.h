#pragma once

#include <core/fitting/finders/finder.h>

namespace DAQuiri
{

class NaiveKON : public Finder
{
 public:
  NaiveKON(const std::vector<double>& x,
      const std::vector<double>& y,
      uint16_t kon_width = 4,
      double sigma = 3.0);

  double find_left(double chan) const override;
  double find_right(double chan) const override;

  //DATA
  uint16_t kon_width_{4};
  double sigma_{3.0};
  std::vector<double> y_convolution;

 private:
  void calc_kon();
  void find_peaks();

  size_t left_edge(size_t idx) const;
  size_t right_edge(size_t idx) const;
};

}
