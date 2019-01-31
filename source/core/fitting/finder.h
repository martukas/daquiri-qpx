#pragma once

#include <core/fitting/fit_settings.h>
#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

struct DetectedPeak
{
  double center;
  double left;
  double right;
  double highest_y {0};
};

class Finder
{
 public:
  Finder() = default;
  Finder(const std::vector<double>& x,
         const std::vector<double>& y,
         const KONSettings& settings);

  void clear();
  void reset();

  bool empty() const;

  bool cloneRange(const Finder& other, double l, double r);
  void setFit(const std::vector<double>& y_fit,
              const std::vector<double>& y_background);
  void find_peaks();
  DetectedPeak tallest_detected() const;

  double find_left(double chan) const;
  double find_right(double chan) const;
  int32_t find_index(double chan_val) const;

  double highest_residual(double l, double r) const;

  //DATA

  std::vector<double> x_;
  std::vector<double> y_;

  SpectrumData weighted_data;

  std::vector<double> y_fit_, y_background_, y_resid_, y_resid_on_background_;
  std::vector<double> fw_theoretical_bin;
  std::vector<double> y_kon, y_convolution;

  std::vector<size_t> prelim;
  std::vector<DetectedPeak> filtered;

  KONSettings settings_;

 private:
  void calc_kon();

  size_t left_edge(size_t idx) const;
  size_t right_edge(size_t idx) const;

  void setNewData(const SpectrumData& d);
};

}
