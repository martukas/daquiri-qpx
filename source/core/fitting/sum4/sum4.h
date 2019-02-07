/**
 * @file sum4.h
 * @brief Construct for simple analysis of a background sample
 *
 * This is an abstraction of a peak area based on:
 * M. Lindstrom, Richard. (1994)
 * Sum and Mean Standard Programs for Activation Analysis.
 * Biological trace element research. 43-45. 597-603.
 * 10.1007/978-1-4757-6025-5_69.
 *
 * @author Martin Shetty
 */

#pragma once

#include <core/fitting/sum4/sum4edge.h>
#include <core/calibration/calibration.h>

namespace DAQuiri
{

class SUM4
{
 public:
  SUM4() = default;

  /// \brief constructs from a subset of data and a background definition,
  ///        calculating all values
  /// \param spectrum_data subset of spectrum data encompassing selected sample
  /// \param LB left background sample
  /// \param RB Right background sample
  SUM4(const WeightedData& spectrum_data, const SUM4Edge& LB, const SUM4Edge& RB);

  /// \returns left channel of peak sample, NaN if uninitialized
  double left() const { return Lchan_; }
  /// \returns right channel of peak sample, NaN if uninitialized
  double right() const { return Rchan_; }
  /// \returns channel width of peak sample, inclusive of end points
  double peak_width() const;

  /// \returns Currie quality indicator; see publication
  // \todo use unsigned type
  int quality() const;

  /// \returns sum of all counts in peak sample, with aggregate uncertainty
  UncertainDouble gross_area() const { return gross_area_; }
  /// \returns estimated background counts in peak sample
  UncertainDouble background_area() const { return background_area_; }
  /// \returns total counts in peak
  UncertainDouble peak_area() const { return peak_area_; }
  /// \returns the channel center for peak
  UncertainDouble centroid() const { return centroid_; }
  /// \returns estimation of peak width
  UncertainDouble fwhm() const { return fwhm_; }

  /// \returns energy calibrated peak center
  UncertainDouble peak_energy(const Calibration& cal) const;
  /// \returns energy calibrated peak width
  UncertainDouble fwhm_energy(const Calibration& cal) const;

  // \todo move this out of the class
  static int get_currie_quality_indicator(double peak_net_area, double background_variance);

  std::string to_string() const;

  friend void to_json(nlohmann::json& j, const SUM4& s);
  friend void from_json(const nlohmann::json& j, SUM4& s);

 private:
  double Lchan_{std::numeric_limits<double>::quiet_NaN()};
  double Rchan_{std::numeric_limits<double>::quiet_NaN()};

  UncertainDouble gross_area_{0.0, 0.0};
  UncertainDouble background_area_{0.0, 0.0};
  UncertainDouble peak_area_{0.0, 0.0};
  UncertainDouble centroid_;
  UncertainDouble fwhm_;
};

}
