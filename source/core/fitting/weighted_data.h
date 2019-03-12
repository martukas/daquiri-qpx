/**
 * @file weighted_data.h
 * @brief
 *
 * This is a class for storing spectrum counts and associated weights for fitting.
 *
 * @author Martin Shetty
 */

#pragma once

#include <vector>
#include <limits>
#include <cstddef>

namespace DAQuiri
{

/// \returns true statistical weight for histogram channel, i.e. sqrt(count)
double weight_true(double count);

/// \returns low-count compensated weight for histogram channel using:
///          G. W. Phillips and K. W. Marlow,
///          "Peak Search and Analysis of Gamma-Ray Spectra with Very Low Statistics,"
///          IEEE Transactions on Nuclear Science, vol. 24, no. 1, pp. 154-157, Feb. 1977.
///          doi: 10.1109/TNS.1977.4328659
/// \param counts vector of channel counts
/// \param index channel index
double weight_phillips_marlow(const std::vector<double>& counts, size_t index);

/// \returns low-count compensated weight for histogram channel using: NEED CITATION!!!
// \todo NEED CITATION!!!
double weight_revay_student(double count);

/// \struct WeightedDataPoint weighted_data.h <core/fitting/weighted_data.h>
/// \brief provides all relevant information from the results of an optimization attempt.
struct WeightedDataPoint
{
  double channel {0};
  double count {0};
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
  /// \param channels channels/bins
  /// \param counts counts
  WeightedData(const std::vector<double>& channels,
               const std::vector<double>& counts);

  /// \returns a subset of data between two provided channel bounds
  /// \param bound1 first bound (min or max)
  /// \param bound2 second bound (min or max)
  WeightedData subset(double bound1, double bound2) const;

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
  double count_min {std::numeric_limits<double>::quiet_NaN()};
  double count_max {std::numeric_limits<double>::quiet_NaN()};
};

}
