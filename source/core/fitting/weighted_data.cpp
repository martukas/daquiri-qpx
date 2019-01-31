#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

WeightedData::WeightedData(const std::vector<double>& x,
                           const std::vector<double>& y)
{
  if (x.size() != y.size())
    throw std::runtime_error("SubSpectrum::set x & y sized don't match");
  data.resize(x.size());
  for (size_t i = 0; i < x.size(); ++i)
  {
    auto& p = data[i];
    p.x = x[i];
    p.y = y[i];
    p.weight_true = weight_true(y, i);
    p.weight_phillips_marlow = weight_phillips_marlow(y, i);
    p.weight_revay = weight_revay_student(y, i);
  }
}

bool WeightedData::empty() const
{
  return data.empty();
}

WeightedData WeightedData::subset(double b1, double b2) const
{
  auto from = std::min(b1, b2);
  auto to = std::max(b1, b2);
  WeightedData ret;
  for (const auto& p : data)
    if ((p.x >= from) && (p.x <= to))
      ret.data.push_back(p);
  return ret;
}

WeightedData WeightedData::left(size_t size) const
{
  size = std::min(size, data.size());
  WeightedData ret;
  ret.data = std::vector<WeightedDataPoint>(data.begin(), data.begin() + size);
  return ret;
}

WeightedData WeightedData::right(size_t size) const
{
  size = std::min(size, data.size());
  WeightedData ret;
  ret.data = std::vector<WeightedDataPoint>(data.begin() + (data.size() - size), data.end());
  return ret;
}

void WeightedData::clear()
{
  data.clear();
}

double WeightedData::weight_true(const std::vector<double>& y, size_t i) const
{
  return std::sqrt(y[i]);
}

double WeightedData::weight_phillips_marlow(const std::vector<double>& y, size_t i) const
{
  double k0 = y[i];

  if (k0 >= 25)
    return std::sqrt(k0);
  else
  {
    k0 = 1.0;
    if ((i > 0) && ((i + 1) < y.size()))
      k0 = y[i - 1] + y[i] + y[i + 1] / 3.0;
    return std::max(std::sqrt(k0), 1.0);
  }
}

double WeightedData::weight_revay_student(const std::vector<double>& y, size_t i) const
{
  double k0 = y[i] + 1;
  return std::sqrt(k0);
}

}
