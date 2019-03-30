/**
 * @file sum4edge.h
 * @brief Construct for simple analysis of a background sample
 *
 * This is an abstraction for background estimation based on:
 * M. Lindstrom, Richard. (1994)
 * Sum and Mean Standard Programs for Activation Analysis.
 * Biological trace element research. 43-45. 597-603.
 * 10.1007/978-1-4757-6025-5_69.
 *
 * Two such samples can be used to generate a polynomial function
 * describing a straight line for background subtraction under a peak.
 *
 * @author Martin Shetty
 */

#pragma once

#include <core/fitting/weighted_data.h>
#include <core/fitting/uncertain.h>
#include <core/calibration/polynomial.h>

namespace DAQuiri
{

struct SUM4Background
{
  double x_offset {0};
  double base {0};
  double slope {0};

  double operator()(double x) const;
};

class SUM4Edge
{
 public:
  SUM4Edge() = default;

  /// \brief constructs from a subset of weighted data, calculating all values
  /// \param spectrum_data subset of spectrum data encompassing selected sample
  SUM4Edge(const WeightedData& spectrum_data);

  /// \returns left channel of edge sample, NaN if uninitialized
  double left() const { return Lchan_; }
  /// \returns right channel of edge sample, NaN if uninitialized
  double right() const { return Rchan_; }
  /// \returns channel width of sample, inclusive of end points
  double width() const;

  /// \returns sum of all counts in sample, with aggregate uncertainty
  UncertainDouble sum() const { return dsum_; }

  /// \returns average number of counts in sample
  UncertainDouble average() const { return davg_; }

  /// \brief convenience function for background estimation
  /// \returns variance of average counts, i.e. square of uncertainty
  double variance() const;

  /// \returns minimum counts in region, NaN if uninitialized
  double min() const { return min_; }

  /// \returns maximum counts in region, NaN if uninitialized
  double max() const { return max_; }

  /// \brief generates a linear function estimating a background spanning two samples
  static SUM4Background sum4_background(const SUM4Edge& LB,
                                        const SUM4Edge& RB);

  std::string to_string() const;

  friend void to_json(nlohmann::json& j, const SUM4Edge& s);
  friend void from_json(const nlohmann::json& j, SUM4Edge& s);

 private:
  double Lchan_{std::numeric_limits<double>::quiet_NaN()};
  double Rchan_{std::numeric_limits<double>::quiet_NaN()};
  double min_{std::numeric_limits<double>::quiet_NaN()};
  double max_{std::numeric_limits<double>::quiet_NaN()};
  UncertainDouble dsum_{0.0, 0.0};
  UncertainDouble davg_{0.0, 0.0};
};

}
