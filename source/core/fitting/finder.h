#pragma once

#include <core/fitting/fit_settings.h>

namespace DAQuiri {

class Finder {

public:
  Finder() {}
  Finder(const std::vector<double> &x, const std::vector<double> &y, const FitSettings &settings);

  void clear();
  void reset();
  bool empty() const;
  
  bool cloneRange(const Finder &other, double l, double r);
  void setFit(const std::vector<double> &x_fit,
              const std::vector<double> &y_fit,
              const std::vector<double> &y_background);
  void find_peaks();

  double find_left(double chan) const;
  double find_right(double chan) const ;
  int32_t find_index(double chan_val) const;

  //DATA

  std::vector<double> x_, y_;
  std::vector<double> y_fit_, y_background_, y_resid_, y_resid_on_background_;
  std::vector<double> fw_theoretical_nrg;
  std::vector<double> fw_theoretical_bin;
  std::vector<double> x_kon, x_conv;

  std::vector<size_t> prelim, filtered, lefts, rights;

  FitSettings settings_;

private:
  void calc_kon();

  size_t left_edge(size_t idx) const;
  size_t right_edge(size_t idx) const;

  void setNewData(const std::vector<double> &x, const std::vector<double> &y);

};

}
