/**
 * @file weighted_data.h
 * @brief
 *
 * This is a class for storing spectrum counts and associated weights for fitting.
 *
 * @author Martin Shetty
 */

#pragma once

#include <cinttypes>
#include <vector>

namespace DAQuiri
{

/// \returns true statistical weight for histogram channel, i.e. sqrt(count)
double weight_true(double count);

/// \returns low-count compensated weight for histogram channel using:
///          G. W. Phillips and K. W. Marlow,
///          "Peak Search and Analysis of Gamma-Ray Spectra with Very Low Statistics,"
///          IEEE Transactions on Nuclear Science, vol. 24, no. 1, pp. 154-157, Feb. 1977.
///          doi: 10.1109/TNS.1977.4328659
/// \param y vector of channel counts
/// \param i channel index
double weight_phillips_marlow(const std::vector<double>& y, size_t i);

/// \returns low-count compensated weight for histogram channel using: NEED CITATION!!!
// \todo NEED CITATION!!!
double weight_revay_student(double count);

/// \struct WeightedDataPoint weighted_data.h <core/fitting/weighted_data.h>
/// \brief provides all relevant information from the results of an optimization attempt.
struct WeightedDataPoint
{
  double x {0}; /// < channel
  double y {0}; /// < counts
  double weight_true {0};            /// < true statistical weight
  double weight_phillips_marlow {0}; /// < low-count compensated weight (need ref)
  double weight_revay {0};           /// < low-count compensated weight (need ref)
};

// \todo uncertainty treatment for ZDT spectra

/// \struct WeightedData weighted_data.h <core/fitting/weighted_data.h>
/// \brief stores histogram data for a region with true and low-count compensated weights
///         and provides interface for extracting sub-regions.
struct WeightedData
{
  WeightedData() = default;

  /// \brief constructs and calculates weights
  /// \param x channels/bins
  /// \param y counts
  WeightedData(const std::vector<double>& x,
               const std::vector<double>& y);

  /// \returns a subset of data between two provided channel bounds
  /// \param b1 first bound (min or max)
  /// \param b2 second bound (min or max)
  WeightedData subset(double b1, double b2) const;

  /// \returns a subset of data from the front end
  /// \param size number of bins from the left to include
  WeightedData left(size_t size) const;

  /// \returns a subset of data from the back end
  /// \param size number of bins from the right to include
  WeightedData right(size_t size) const;

  /// \brief clears data
  void clear();

  /// \returns if data is empty
  bool empty() const;

  std::vector<WeightedDataPoint> data;
};

}
