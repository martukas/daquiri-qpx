/**
 * @file weighted_data.h
 * @brief
 *
 *
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
std::vector<double> weight_true(const std::vector<double>& counts);

/// \returns low-count compensated weight for histogram channel using:
///          G. W. Phillips and K. W. Marlow,
///          "Peak Search and Analysis of Gamma-Ray Spectra with Very Low Statistics,"
///          IEEE Transactions on Nuclear Science, vol. 24, no. 1, pp. 154-157, Feb. 1977.
///          doi: 10.1109/TNS.1977.4328659
/// \param counts vector of channel counts
/// \param index channel index
double weight_phillips_marlow(const std::vector<double>& counts, size_t index);
std::vector<double> weight_phillips_marlow(const std::vector<double>& counts);

/// \returns low-count compensated weight for histogram channel using: NEED CITATION!!!
// \todo NEED CITATION!!!
double weight_revay_student(double count);
std::vector<double> weight_revay_student(const std::vector<double>& counts);


}
