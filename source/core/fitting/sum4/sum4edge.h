#pragma once

#include <core/fitting/weighted_data.h>
#include <core/fitting/uncertain.h>
#include <core/calibration/polynomial.h>

namespace DAQuiri
{

class SUM4Edge
{
  double Lchan_{std::numeric_limits<double>::quiet_NaN()};
  double Rchan_{std::numeric_limits<double>::quiet_NaN()};
  double min_{std::numeric_limits<double>::quiet_NaN()};
  double max_{std::numeric_limits<double>::quiet_NaN()};
  UncertainDouble dsum_{0.0, 0.0};
  UncertainDouble davg_{0.0, 0.0};

 public:
  SUM4Edge() = default;
  SUM4Edge(const WeightedData& d);

  double left() const { return Lchan_; }
  double right() const { return Rchan_; }
  double width() const;
  UncertainDouble sum() const { return dsum_; }
  UncertainDouble average() const { return davg_; }
  double variance() const;

  double min() const { return min_; }
  double max() const { return max_; }
  double midpoint() const;

  std::string to_string() const;

  friend void to_json(nlohmann::json& j, const SUM4Edge& s);
  friend void from_json(const nlohmann::json& j, SUM4Edge& s);

  static Polynomial sum4_background(const SUM4Edge& L,
                                    const SUM4Edge& R);
};

}
