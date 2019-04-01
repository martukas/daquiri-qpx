#include <core/gamma/finders/finder_kon_calibrated.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

CalibratedKON::CalibratedKON(const std::vector<double>& x,
                             const std::vector<double>& y,
                             const FCalibration& cal,
                             double sigma,
                             double edge_width_factor)
    : Finder(x, y)
      , sigma_(sigma)
      , edge_width_factor_(edge_width_factor)
{
  if (cal.cali_fwhm_.valid() && cal.cali_nrg_.valid())
  {
    for (size_t i=0; i < x_.size(); ++i)
    {
      double wchan = cal.bin_to_width(x_[i]);
      if (!std::isfinite(wchan))
        DBG("Kon {} - > {}", x_[i], wchan);
      fw_theoretical.push_back(wchan);
    }
  }

  calc_kon();
  find_peaks();
}

void CalibratedKON::calc_kon()
{
  uint16_t width;
  int32_t start {0};
  int32_t end {0};
  int32_t shift;

  for (size_t i = 0; i < fw_theoretical.size(); ++i)
    if (ceil(fw_theoretical[i]) < i)
    {
      start = i;
      break;
    }

  for (int32_t i = fw_theoretical.size() - 1; i >= 0; --i)
    if (2 * ceil(fw_theoretical[i]) + i + 1 < fw_theoretical.size())
    {
      end = i;
      break;
    }

  y_convolution.assign(y_.size(), 0.0);

  for (int32_t j = start; j < end; ++j)
  {
    width = static_cast<uint16_t>(std::floor(fw_theoretical[j]));
    shift = width / uint16_t(2);

    double kon{0.0};
    double avg{0.0};
    for (int32_t i = j; i <= (j + width + 1); ++i)
    {
      kon += 2 * y_[i] - y_[i - width] - y_[i + width];
      avg += y_[i];
    }
    avg = avg / width;
    y_convolution[j + shift] = kon / sqrt(6.0 * width * avg);
  }
}

void CalibratedKON::find_peaks()
{
  calc_kon();
  filtered.clear();

  std::vector<size_t> prelim;
  for (size_t j = 0; j < y_convolution.size(); ++j)
    if (y_convolution[j] > sigma_)
      prelim.push_back(j);

  if (prelim.empty())
    return;

  //find edges of contiguous peak areas
  std::vector<size_t> lefts;
  std::vector<size_t> rights;
  lefts.push_back(prelim[0]);
  size_t prev = prelim[0];
  for (const auto& current : prelim)
  {
    if ((current - prev) > 1)
    {
      rights.push_back(prev);
      lefts.push_back(current);
    }
    prev = current;
  }
  rights.push_back(prev);

  //assume center is bentween edges
  for (size_t i = 0; i < lefts.size(); ++i)
  {
    DetectedPeak p;
    size_t l = left_edge(lefts[i]);
    size_t r = right_edge(rights[i]);
    // \todo use function?
    for (size_t j = l; j <= r; ++j)
      p.highest_y = std::max(p.highest_y, y_[j]);
    p.left = x_[l];
    p.right = x_[r];
    p.center = 0.5 * (p.left + p.right);
    filtered.push_back(p);
  }
}

double CalibratedKON::find_left(double chan) const
{
  if (x_.empty())
    return 0;

  //assume x is monotone increasing

  if ((chan < x_[0]) || (chan >= x_[x_.size() - 1]))
    return x_.front();

  int i = x_.size() - 1;
  while ((i > 0) && (x_[i] > chan))
    i--;

  return x_[left_edge(i)];
}

double CalibratedKON::find_right(double chan) const
{
  if (x_.empty())
    return 0;

  //assume x is monotone increasing

  if ((chan < x_[0]) || (chan >= x_[x_.size() - 1]))
    return x_.back();

  size_t i = 0;
  while ((i < x_.size()) && (x_[i] < chan))
    i++;

  return x_[right_edge(i)];
}

size_t CalibratedKON::left_edge(size_t idx) const
{
  if (y_convolution.empty() || idx >= y_convolution.size())
    return 0;

  double width = floor(fw_theoretical[idx]);
  double goal = x_[idx] - 0.5 * width * edge_width_factor_;
  while ((idx > 0) && (x_[idx] > goal))
    idx--;
  return idx;
}

size_t CalibratedKON::right_edge(size_t idx) const
{
  if (y_convolution.empty() || idx >= y_convolution.size())
    return 0;

  double width = floor(fw_theoretical[idx]);
  double goal = x_[idx] + 0.5 * width * edge_width_factor_;
  while ((idx < x_.size()) && (x_[idx] < goal))
    idx++;
  return idx;
}

}
