#pragma once

#include <core/fitting/fit_settings.h>

namespace DAQuiri
{

struct SpectrumDataPoint
{
  double x {0};
  double y {0};
  double weight_true {0};
  double weight_phillips_marlow {0};
  double weight_revay {0};
};

struct SpectrumData
{
  // \todo uncertainty treatment for ZDT spectra
  // \todo empty
  // \todo push_back

  SpectrumData() = default;
  void set(const std::vector<double>& x,
           const std::vector<double>& y);
  SpectrumData subset(double from, double to) const;
  SpectrumData left(size_t size) const;
  SpectrumData right(size_t size) const;
  void clear();

  std::vector<SpectrumDataPoint> data;

  double weight_true(const std::vector<double>& y, size_t i) const;
  double weight_phillips_marlow(const std::vector<double>& y, size_t i) const;
  double weight_revay_student(const std::vector<double>& y, size_t i) const;
};

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
  void setFit(const std::vector<double>& x_fit,
              const std::vector<double>& y_fit,
              const std::vector<double>& y_background);
  void find_peaks();
  DetectedPeak tallest_detected() const;

  double find_left(double chan) const;
  double find_right(double chan) const;
  int32_t find_index(double chan_val) const;

  //DATA

  std::vector<double> x_;
  std::vector<double> y_;

  SpectrumData weighted_data;

  std::vector<double> y_fit_, y_background_, y_resid_, y_resid_on_background_;
  std::vector<double> fw_theoretical_nrg;
  std::vector<double> fw_theoretical_bin;
  std::vector<double> y_kon, y_convolution;

  std::vector<size_t> prelim;
  std::vector<DetectedPeak> filtered;

  KONSettings settings_;

 private:
  void calc_kon();

  size_t left_edge(size_t idx) const;
  size_t right_edge(size_t idx) const;

  void setNewData(const std::vector<double>& x, const std::vector<double>& y);
};

}
