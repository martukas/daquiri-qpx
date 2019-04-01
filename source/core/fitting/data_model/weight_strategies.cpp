/**
 * @file weighted_data.h
 * @brief
 *
 *
 *
 * @author Martin Shetty
 */

#include <core/fitting/data_model/weighted_data.h>
#include <stdexcept>
#include <cmath>
#include <range/v3/all.hpp>

namespace DAQuiri
{

double weight_true(double count)
{
  return std::sqrt(count);
}

std::vector<double> weight_true(const std::vector<double>& counts)
{
  return counts | ranges::view::transform([](double i) { return weight_true(i); });
}


double weight_phillips_marlow(const std::vector<double>& counts, size_t index)
{
  double k0 = counts[index];

  if (k0 >= 25.)
    return std::sqrt(k0);
  else
  {
    k0 = 1.;
    if ((index > 0) && ((index + 1) < counts.size()))
      k0 = (counts[index - 1] + counts[index] + counts[index + 1]) / 3.;
    // \todo what if on edges?
    return std::max(std::sqrt(k0), 1.);
  }
}

std::vector<double> weight_phillips_marlow(const std::vector<double>& counts)
{
  std::vector<double> ret;
  ret.resize(counts.size());
  for (size_t i = 0; i < counts.size(); ++i)
    ret[i] = weight_phillips_marlow(counts, i);
  return ret;
}


double weight_revay_student(double count)
{
  return std::sqrt(count + 1.);
}

std::vector<double> weight_revay_student(const std::vector<double>& counts)
{
  return counts | ranges::view::transform([](double i) { return weight_revay_student(i); });
}

}
