#pragma once

#include <core/gamma/finders/finder.h>
//#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

class CalibratedKON : public Finder
{
 public:
  CalibratedKON(const std::vector<double>& x,
                const std::vector<double>& y,
                const FCalibration& cal,
                double sigma = 3.0,
                double edge_width_factor = 3.5);

  double find_left(double chan) const override;
  double find_right(double chan) const override;

  //DATA
  double sigma_{3.0};
  double edge_width_factor_{3.5};
  std::vector<double> fw_theoretical;
  std::vector<double> y_convolution;

 private:
  void calc_kon();
  void find_peaks();

  size_t left_edge(size_t idx) const;
  size_t right_edge(size_t idx) const;
};

}
