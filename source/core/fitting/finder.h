#pragma once

#include <core/fitting/fit_settings.h>

namespace DAQuiri
{

class Finder
{
 public:
  Finder() = default;
  Finder(const std::vector<double>& x,
         const std::vector<double>& y,
         const FitSettings& settings);

  // \todo constructors for ZDT spectra

  void clear();
  void reset();
  void calc_uncertainties();

  bool empty() const;

  bool cloneRange(const Finder& other, double l, double r);
  void setFit(const std::vector<double>& x_fit,
              const std::vector<double>& y_fit,
              const std::vector<double>& y_background);
  void find_peaks();

  double find_left(double chan) const;
  double find_right(double chan) const;
  int32_t find_index(double chan_val) const;

  //DATA

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> y_weight_true;
  std::vector<double> y_weight_phillips_marlow;
  std::vector<double> y_weight_revay;
  std::vector<double> y_fit_, y_background_, y_resid_, y_resid_on_background_;
  std::vector<double> fw_theoretical_nrg;
  std::vector<double> fw_theoretical_bin;
  std::vector<double> y_kon, y_convolution;

  std::vector<size_t> prelim, filtered, lefts, rights;

  FitSettings settings_;

 private:
  void calc_kon();

  size_t left_edge(size_t idx) const;
  size_t right_edge(size_t idx) const;

  void setNewData(const std::vector<double>& x, const std::vector<double>& y);

  double weight_true(size_t i) const;
  double weight_phillips_marlow(size_t i) const;
  double weight_revay_student(size_t i) const;
};

}
