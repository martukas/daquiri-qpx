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

  if (k0 >= 25)
    return std::sqrt(k0);
  else
  {
    k0 = 1.0;
    if ((index > 0) && ((index + 1) < counts.size()))
      k0 = (counts[index - 1] + counts[index] + counts[index + 1]) / 3.0;
    // \todo what if on edges?
    return std::max(std::sqrt(k0), 1.0);
  }
}

double weight_revay_student(double count)
{
  return std::sqrt(count + 1.0);
}

WeightedData::WeightedData(const std::vector<double>& channels,
                           const std::vector<double>& counts)
{
  if (channels.size() != counts.size())
    throw std::runtime_error("WeightedData::constructor x & y sizes don't match");
  if (channels.empty())
    throw std::runtime_error("WeightedData::constructor data empty");
  data.resize(channels.size());
  count_min = count_max = counts[0];
  for (size_t i = 0; i < channels.size(); ++i)
  {
    auto& p = data[i];
    p.channel = channels[i];
    p.count = counts[i];
    p.weight_true = weight_true(p.count);
    p.weight_phillips_marlow = weight_phillips_marlow(counts, i);
    p.weight_revay = weight_revay_student(p.count);
    count_min = std::min(count_min, p.count);
    count_max = std::max(count_max, p.count);
  }
}

bool WeightedData::empty() const
{
  return data.empty();
}

WeightedData WeightedData::subset(double bound1, double bound2) const
{
  auto from = std::min(bound1, bound2);
  auto to = std::max(bound1, bound2);
  WeightedData ret;
  for (const auto& p : data)
    if ((p.channel >= from) && (p.channel <= to))
    {
      if (std::isnan(ret.count_min))
        ret.count_min = ret.count_max = p.count;
      ret.count_min = std::min(ret.count_min, p.count);
      ret.count_max = std::max(ret.count_max, p.count);
      ret.data.push_back(p);
    }
  return ret;
}

WeightedData WeightedData::left(size_t size) const
{
  size = std::min(size, data.size());
  WeightedData ret;
  ret.data = std::vector<WeightedDataPoint>(data.begin(), data.begin() + size);
  for (const auto& p : ret.data)
  {
    if (std::isnan(ret.count_min))
      ret.count_min = ret.count_max = p.count;
    ret.count_min = std::min(ret.count_min, p.count);
    ret.count_max = std::max(ret.count_max, p.count);
  }
  return ret;
}

WeightedData WeightedData::right(size_t size) const
{
  size = std::min(size, data.size());
  WeightedData ret;
  ret.data = std::vector<WeightedDataPoint>(data.begin() + (data.size() - size), data.end());
  for (const auto& p : ret.data)
  {
    if (std::isnan(ret.count_min))
      ret.count_min = ret.count_max = p.count;
    ret.count_min = std::min(ret.count_min, p.count);
    ret.count_max = std::max(ret.count_max, p.count);
  }
  return ret;
}

void WeightedData::clear()
{
  data.clear();
}

}
