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

WeightedData::WeightedData(const std::vector<double>& channels,
                           const std::vector<double>& counts)
{
  if (channels.size() != counts.size())
    throw std::runtime_error("WeightedData::constructor x & y sizes don't match");
  if (channels.empty())
    throw std::runtime_error("WeightedData::constructor data empty");

  chan = channels;
  count = counts;

//  data.resize(channels.size());
//  for (size_t i = 0; i < channels.size(); ++i)
//  {
//    auto& p = data[i];
//    p.channel = channels[i];
//    p.count = counts[i];
//  }
}

bool WeightedData::valid() const
{
  return (!chan.empty() &&
      (chan.size() == count.size()) &&
      (chan.size() == count_weight.size()));
}

//bool WeightedData::empty() const
//{
//  return data.empty();
//}

double WeightedData::count_min() const
{
  double min = std::numeric_limits<double>::quiet_NaN();
  for (const auto& p : count)
  {
    if (std::isnan(min))
      min = p;
    min = std::min(min, p);
  }
  return min;
}

double WeightedData::count_max() const
{
  double max = std::numeric_limits<double>::quiet_NaN();
  for (const auto& p : count)
  {
    if (std::isnan(max))
      max = p;
    max = std::max(max, p);
  }
  return max;
}

WeightedData WeightedData::subset(double bound1, double bound2) const
{
  auto from = std::min(bound1, bound2);
  auto to = std::max(bound1, bound2);
//  auto is_in_range = [from, to](WeightedDataPoint p)
//  {
//    return ((p.channel >= from) && (p.channel <= to));
//  };

  auto is_in_range2 = [from, to](ranges::v3::common_tuple<const double&, const double&, const double&> p)
  {
    return ((std::get<0>(p) >= from) && (std::get<0>(p) <= to));
  };

  auto zipped = ranges::view::filter(ranges::view::zip(chan, count, count_weight), is_in_range2);

  WeightedData ret;
//  ret.data = ranges::view::filter(data, is_in_range);
  for (const auto& z : zipped)
  {
    ret.chan.push_back(std::get<0>(z));
    ret.count.push_back(std::get<1>(z));
    ret.count_weight.push_back(std::get<2>(z));
  }
  return ret;
}

WeightedData WeightedData::left(size_t size) const
{
  size = std::min(size, chan.size());
  WeightedData ret;
//  ret.data = ranges::view::take(data, size);
  ret.chan = ranges::view::take(chan, size);
  ret.count = ranges::view::take(count, size);
  ret.count_weight = ranges::view::take(count_weight, size);
  return ret;
}

WeightedData WeightedData::right(size_t size) const
{
  size = std::min(size, chan.size());
  WeightedData ret;
//  ret.data = ranges::view::drop(data, data.size() - size);
  ret.chan = ranges::view::drop(chan, chan.size() - size);
  ret.count = ranges::view::drop(count, count.size() - size);
  ret.count_weight = ranges::view::drop(count_weight, count_weight.size() - size);
  return ret;
}

void WeightedData::clear()
{
  chan.clear();
  count.clear();
  count_weight.clear();
//  data.clear();
}

}
