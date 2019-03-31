/**
 * @file weighted_data.h
 * @brief
 *
 * This is a class for storing spectrum counts and associated weights for fitting.
 *
 * @author Martin Shetty
 */

#include <core/fitting/weighted_data.h>
#include <stdexcept>
#include <cmath>

namespace DAQuiri
{

double weight_true(double count)
{
  return std::sqrt(count);
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

double weight_revay_student(double count)
{
  return std::sqrt(count + 1.);
}

WeightedData::WeightedData(const std::vector<double>& channels,
                           const std::vector<double>& counts)
{
  if (channels.size() != counts.size())
    throw std::runtime_error("WeightedData::constructor x & y sizes don't match");
  if (channels.empty())
    throw std::runtime_error("WeightedData::constructor data empty");

  using namespace ranges;
  weights = counts | ranges::view::transform([](double i) { return weight_revay_student(i); });

  data.resize(channels.size());
  for (size_t i = 0; i < channels.size(); ++i)
  {
    auto& p = data[i];
    p.channel = channels[i];
    p.count = counts[i];
    p.weight_true = weight_true(p.count);
    p.weight_phillips_marlow = weight_phillips_marlow(counts, i);
    p.weight_revay = weight_revay_student(p.count);
  }
}

bool WeightedData::empty() const
{
  return data.empty();
}

double WeightedData::count_min() const
{
  double min = std::numeric_limits<double>::quiet_NaN();
  for (const auto& p : data)
  {
    if (std::isnan(min))
      min = p.count;
    min = std::min(min, p.count);
  }
  return min;
}

double WeightedData::count_max() const
{
  double max = std::numeric_limits<double>::quiet_NaN();
  for (const auto& p : data)
  {
    if (std::isnan(max))
      max = p.count;
    max = std::max(max, p.count);
  }
  return max;
}

WeightedData WeightedData::subset(double bound1, double bound2) const
{
  auto from = std::min(bound1, bound2);
  auto to = std::max(bound1, bound2);
  auto is_in_range = [from, to](WeightedDataPoint p)
      { return ((p.channel >= from) && (p.channel <= to)); };
  WeightedData ret;
  ret.data = ranges::view::filter(data, is_in_range);
  return ret;
}

WeightedData WeightedData::left(size_t size) const
{
  size = std::min(size, data.size());
  WeightedData ret;
  ret.data = ranges::view::take(data, size);
  return ret;
}

WeightedData WeightedData::right(size_t size) const
{
  size = std::min(size, data.size());
  WeightedData ret;
  ret.data = ranges::view::drop(data, data.size() - size);
  return ret;
}

void WeightedData::clear()
{
  data.clear();
}

}
