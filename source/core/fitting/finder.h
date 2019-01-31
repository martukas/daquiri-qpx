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
  double highest_y{0};
};

class KON
{
 public:
  KON(const std::vector<double>& x,
      const std::vector<double>& y,
      bool residuals,
      const KONSettings& settings);

  DetectedPeak tallest_detected() const;
  double find_left(double chan) const;
  double find_right(double chan) const;
  double highest_residual(double l, double r) const;

  //DATA
  KONSettings settings_;

  bool residuals_ {false};

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> fw_theoretical_bin;
  std::vector<double> y_convolution;

  std::vector<DetectedPeak> filtered;

 private:
  void calc_kon();
  void find_peaks();

  size_t left_edge(size_t idx) const;
  size_t right_edge(size_t idx) const;
  int32_t find_index(double chan_val) const;
};


class FitEvaluation
{
 public:
  FitEvaluation() = default;
  FitEvaluation(const WeightedData& data);

  void clear();
  void reset();

  bool empty() const;

  bool cloneRange(const FitEvaluation& other, double l, double r);
  void update_fit(const std::vector<double>& y_fit,
                  const std::vector<double>& y_background);

  //DATA

  std::vector<double> x_;
  std::vector<double> y_;

  WeightedData weighted_data;

  std::vector<double> y_fit_, y_background_, y_resid_, y_resid_on_background_;

 private:

  void setNewData(const WeightedData& d);
};

}
