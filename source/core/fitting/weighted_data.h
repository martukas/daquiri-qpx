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

  SpectrumData() = default;
  SpectrumData(const std::vector<double>& x,
               const std::vector<double>& y);
  SpectrumData subset(double b1, double b2) const;
  SpectrumData left(size_t size) const;
  SpectrumData right(size_t size) const;


  void clear();
  bool empty() const;

  std::vector<SpectrumDataPoint> data;

  double weight_true(const std::vector<double>& y, size_t i) const;
  double weight_phillips_marlow(const std::vector<double>& y, size_t i) const;
  double weight_revay_student(const std::vector<double>& y, size_t i) const;
};

}
