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

  /// \returns if data is nonempty and array sizes match
  bool valid() const;

  std::vector<double> chan;
  std::vector<double> count;
  std::vector<double> count_weight;

  double count_min() const;
  double count_max() const;
};

}
