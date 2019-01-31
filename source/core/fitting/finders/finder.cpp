#include <core/fitting/finders/finder.h>

namespace DAQuiri
{

double Finder::highest_residual(double l, double r) const
{
  auto li = find_index(l);
  auto ri = find_index(r);
  if ((li < 0) || (ri < 0))
    return 0.0;
  double ret{0.0};
  for (size_t j = li; j <= ri; ++j)
    ret = std::max(ret, y_[j]);
}

DetectedPeak Finder::tallest_detected() const
{
  DetectedPeak p;
  for (const auto& pp : filtered)
    if (pp.highest_y > p.highest_y)
      p = pp;
  return p;
}

int32_t Finder::find_index(double chan_val) const
{
  if (x_.empty())
    return -1;

  if (chan_val <= x_[0])
    return 0;

  if (chan_val >= x_[x_.size() - 1])
    return x_.size() - 1;

  size_t i = 0;
  while ((i < x_.size()) && (x_[i] < chan_val))
    i++;

  return i;
}

}
