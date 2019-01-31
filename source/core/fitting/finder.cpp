#include <core/fitting/finder.h>

namespace DAQuiri
{

Finder::Finder(const std::vector<double>& x, const std::vector<double>& y, const KONSettings& settings)
: settings_(settings)
{
  setNewData(SpectrumData(x, y));
}

bool Finder::cloneRange(const Finder& other, double l, double r)
{
  setNewData(other.weighted_data.subset(l, r));
}

void Finder::setNewData(const SpectrumData& d)
{
  clear();
  weighted_data = d;
  for (const auto& p : d.data)
  {
    x_.push_back(p.x);
    y_.push_back(p.y);
  }
  reset();
  calc_kon();
  find_peaks();
}

void Finder::clear()
{
  x_.clear();
  y_.clear();
  y_fit_.clear();
  y_background_.clear();
  y_resid_.clear();
  y_resid_on_background_.clear();
//  settings_.clear();

  weighted_data.clear();

  prelim.clear();
  filtered.clear();

  y_kon.clear();
  y_convolution.clear();
}

void Finder::reset()
{
  y_resid_on_background_ = y_resid_ = y_;
  y_fit_.assign(x_.size(), 0.0);
  y_background_.assign(x_.size(), 0.0);
}

bool Finder::empty() const
{
  return x_.empty();
}

void Finder::setFit(const std::vector<double>& y_fit,
                    const std::vector<double>& y_background)
{
  if ((y_fit.size() != y_background.size())
      || (y_fit.empty()))
    return;

  for (size_t i = 0; i < y_fit.size(); ++i)
  {
    y_fit_[i] = y_fit[i];
    y_background_[i] = y_background[i];
    double resid = y_[i] - y_fit[i];
    y_resid_[i] = resid;
    y_resid_on_background_[i] = y_background[i] + resid;
  }

  calc_kon();
  find_peaks();

//  if (y_fit.size() == y_.size()) {
//    y_fit_ = y_fit;

//    y_resid_ = y_;
//    y_resid_on_background_ = y_background;
//    for (int i=0; i < y_.size(); ++i) {
//      y_resid_[i] = y_[i] - y_fit_[i];
//      y_resid_on_background_[i] += y_resid_[i];
//    }

//    calc_kon();
//    find_peaks();
//  }
}

void Finder::calc_kon()
{
  fw_theoretical_bin.clear();

//  if (cal.cali_fwhm_.valid() && cal.cali_nrg_.valid())
//  {
//    for (size_t i=0; i < x_.size(); ++i)
//    {
//      double wchan = cal.bin_to_width(x_[i]);
//      if (!std::isfinite(wchan))
//        DBG("Kon {} - > {}", x_[i], wchan);
//      fw_theoretical_bin.push_back(wchan);
//    }
//  }

  uint16_t width = settings_.width;

  if (width < 2)
    width = 2;

  double sigma = settings_.sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.sigma_resid;
  }


//  DBG << "<Finder> width " << settings_.width;

  int start = width;
  int end = x_.size() - 1 - 2 * width;
  int shift = width / 2;

  if (!fw_theoretical_bin.empty())
  {
    for (size_t i = 0; i < fw_theoretical_bin.size(); ++i)
      if (ceil(fw_theoretical_bin[i]) < i)
      {
        start = i;
        break;
      }

    for (int i = fw_theoretical_bin.size() - 1; i >= 0; --i)
      if (2 * ceil(fw_theoretical_bin[i]) + i + 1 < fw_theoretical_bin.size())
      {
        end = i;
        break;
      }
  }

  y_kon.assign(y_resid_.size(), 0.0);
  y_convolution.assign(y_resid_.size(), 0.0);
  prelim.clear();

  for (int j = start; j < end; ++j)
  {
    if (!fw_theoretical_bin.empty())
    {
      width = floor(fw_theoretical_bin[j]);
      shift = width / 2;
    }

    double kon = 0;
    double avg = 0;
    for (int i = j; i <= (j + width + 1); ++i)
    {
      kon += 2 * y_resid_[i] - y_resid_[i - width] - y_resid_[i + width];
      avg += y_resid_[i];
    }
    avg = avg / width;
    y_kon[j + shift] = kon;
    y_convolution[j + shift] = kon / sqrt(6 * width * avg);

    if (y_convolution[j + shift] > sigma)
      prelim.push_back(j + shift);
  }
}

void Finder::find_peaks()
{
  calc_kon();
  filtered.clear();

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
    for (size_t j = l; j <= r; ++j)
      p.highest_y = std::max(p.highest_y, y_resid_[j]);
    p.left = x_[l];
    p.right = x_[r];
    p.center = 0.5 * (p.left + p.right);
    filtered.push_back(p);
  }
}

double Finder::highest_residual(double l, double r) const
{
  auto li = find_index(l);
  auto ri = find_index(r);
  if ((li < 0) || (ri < 0))
    return 0.0;
  double ret {0.0};
  for (size_t j = li; j <= ri; ++j)
    ret = std::max(ret, y_resid_[j]);
}

DetectedPeak Finder::tallest_detected() const
{
  DetectedPeak p;
  for (const auto& pp : filtered)
    if (pp.highest_y > p.highest_y)
      p = pp;
  return p;
}

double Finder::find_left(double chan) const
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

double Finder::find_right(double chan) const
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

size_t Finder::left_edge(size_t idx) const
{
  if (y_convolution.empty() || idx >= y_convolution.size())
    return 0;

  if (!fw_theoretical_bin.empty())
  {
    double width = floor(fw_theoretical_bin[idx]);
    double goal = x_[idx] - 0.5 * width * settings_.edge_width_factor;
    while ((idx > 0) && (x_[idx] > goal))
      idx--;
    return idx;
  }

  double sigma = settings_.sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.sigma_resid;
  }

  double edge_threshold = -0.5 * sigma;

  while ((idx > 0) && (y_convolution[idx] >= 0))
    idx--;
  if (idx > 0)
    idx--;
  while ((idx > 0) && (y_convolution[idx] < edge_threshold))
    idx--;

  return idx;
}

size_t Finder::right_edge(size_t idx) const
{
  if (y_convolution.empty() || idx >= y_convolution.size())
    return 0;

  if (!fw_theoretical_bin.empty())
  {
    double width = floor(fw_theoretical_bin[idx]);
    double goal = x_[idx] + 0.5 * width * settings_.edge_width_factor;
    while ((idx < x_.size()) && (x_[idx] < goal))
      idx++;
    return idx;
  }

  double sigma = settings_.sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.sigma_resid;
  }

  double edge_threshold = -0.5 * sigma;

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
