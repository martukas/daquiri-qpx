#pragma once

#include <cinttypes>
#include <vector>

namespace DAQuiri
{

struct WeightedDataPoint
{
  double x {0};
  double y {0};
  double weight_true {0};
  double weight_phillips_marlow {0};
  double weight_revay {0};
};

struct WeightedData
{
  // \todo uncertainty treatment for ZDT spectra

  WeightedData() = default;
  WeightedData(const std::vector<double>& x,
               const std::vector<double>& y);
  WeightedData subset(double b1, double b2) const;
  WeightedData left(size_t size) const;
  WeightedData right(size_t size) const;


  void clear();
  bool empty() const;

  std::vector<WeightedDataPoint> data;

  double weight_true(const std::vector<double>& y, size_t i) const;
  double weight_phillips_marlow(const std::vector<double>& y, size_t i) const;
  double weight_revay_student(const std::vector<double>& y, size_t i) const;
};

}
