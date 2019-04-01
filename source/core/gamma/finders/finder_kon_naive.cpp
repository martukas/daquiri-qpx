#include <core/gamma/finders/finder_kon_naive.h>

namespace DAQuiri
{

NaiveKON::NaiveKON(const std::vector<double>& x,
         const std::vector<double>& y,
         uint16_t kon_width,
         double sigma)
    : Finder(x, y)
    , kon_width_(std::max(kon_width, uint16_t(2)))
    , sigma_(sigma)
{
  if (x_.empty())
    return;
  calc_kon();
  find_peaks();
}

void NaiveKON::calc_kon()
{
  size_t start = kon_width_;
  size_t end = x_.size() - 1 - 2 * kon_width_;
  size_t shift = kon_width_ / uint16_t(2);

  y_convolution.assign(y_.size(), 0.0);

  for (size_t j = start; j < end; ++j)
  {
    double kon{0.0};
    double avg{0.0};
    for (size_t i = j; i <= (j + kon_width_ + 1); ++i)
    {
      kon += 2 * y_[i] - y_[i - kon_width_] - y_[i + kon_width_];
      avg += y_[i];
    }
    avg = avg / kon_width_;
    y_convolution[j + shift] = kon / sqrt(6.0 * kon_width_ * avg);
  }
}

void NaiveKON::find_peaks()
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

double NaiveKON::find_left(double chan) const
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

double NaiveKON::find_right(double chan) const
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

size_t NaiveKON::left_edge(size_t idx) const
{
  if (y_convolution.empty() || idx >= y_convolution.size())
    return 0;

  double edge_threshold = -0.5 * sigma_;

  while ((idx > 0) && (y_convolution[idx] >= 0))
    idx--;
  if (idx > 0)
    idx--;
  while ((idx > 0) && (y_convolution[idx] < edge_threshold))
    idx--;

  return idx;
}

size_t NaiveKON::right_edge(size_t idx) const
{
  if (y_convolution.empty() || idx >= y_convolution.size())
    return 0;

  double edge_threshold = -0.5 * sigma_;

  while ((idx < y_convolution.size()) && (y_convolution[idx] >= 0))
    idx++;
  if (idx < y_convolution.size())
    idx++;
  while ((idx < y_convolution.size()) && (y_convolution[idx] < edge_threshold))
    idx++;

  if (idx >= y_convolution.size())
    idx = y_convolution.size() - 1;

  return idx;
}

/*
//Laszlo's implementation

int32_t Region::L(int32_t i, int32_t j, int32_t m)
{
  // \todo what does this do?
  //m = FWHM
  int32_t ret{0};

  if (j - m <= i && i <= j - 1)
    ret = -1;
  if (j <= i && i <= j + m - 1)
    ret = 2;
  if (j + m <= i && i <= j + 2 * m - 1)
    ret = -1;

  return ret;
}

void Region::find_peaks(uint8_t threshold)
{
  auto m = static_cast<int32_t>(1.6551 * default_peak_.width_.val());
  if (m <= 0)
    m = 3;

  size_t i;
  for (size_t j = first_channel; j <= last_channel; ++j)
  {
    double val = 0;
    double var = 0;
    for (i = j - m; i <= (j + 2 * m - 1); ++j)
      if (i > 1)
        val = val + L(i, j, m) * spectrum.channels[i];
    var += square(L(i, j, m)) * spectrum.channels[i];

    //Conv(j - FirstChannel) = val / std::sqrt(Var)
    //if(((Conv(j - FirstChannel - 2) < Conv(j - FirstChannel - 1)) &&
    //(Conv(j - FirstChannel) < Conv(j - FirstChannel - 1)) &&
    //(Conv(j - FirstChannel - 1) > Threshold))) {
    //AddPeak(j - 1, j - 2, j, std::sqrt(spectrum.Channel[j]))
  }
}
*/


}
