#include <core/fitting/finder.h>

namespace DAQuiri
{

Finder::Finder(const std::vector<double>& x, const std::vector<double>& y, const FitSettings& settings)
{
  settings_ = settings;
  setNewData(x, y);
}

void Finder::setNewData(const std::vector<double>& x, const std::vector<double>& y)
{
  clear();
  if (x.size() == y.size())
  {
    x_ = x;
    y_ = y;
    reset();
    calc_uncertainties();

    calc_kon();
    find_peaks();
  }
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

  y_weight_true.clear();
  y_weight_phillips_marlow.clear();
  y_weight_revay.clear();

  prelim.clear();
  filtered.clear();
  lefts.clear();
  rights.clear();

  y_kon.clear();
  y_convolution.clear();
}

void Finder::reset()
{
  y_resid_on_background_ = y_resid_ = y_;
  y_fit_.resize(x_.size(), 0);
  y_background_.resize(x_.size(), 0);
  y_weight_true.resize(x_.size(), 0);
  y_weight_phillips_marlow.resize(x_.size(), 0);
  y_weight_revay.resize(x_.size(), 0);
}

void Finder::calc_uncertainties()
{
  for (size_t i = 0; i < y_.size(); ++i)
  {
    y_weight_true[i] = weight_true(i);
    y_weight_phillips_marlow[i] = weight_phillips_marlow(i);
    y_weight_revay[i] = weight_revay_student(i);
  }
}

bool Finder::empty() const
{
  return x_.empty();
}

bool Finder::cloneRange(const Finder& other, double l, double r)
{
  if (other.x_.empty()
      || other.y_.empty()
      || other.x_.size() != other.y_.size())
    return false;

  size_t min = other.find_index(l);
  size_t max = other.find_index(r);

  if (min >= other.x_.size())
    min = other.x_.size() - 1;
  if (max >= other.x_.size())
    max = other.x_.size() - 1;

  std::vector<double> x_local, y_local;
  for (size_t i = min; i < max; ++i)
  {
    x_local.push_back(other.x_[i]);
    y_local.push_back(other.y_[i]);
  }
  setNewData(x_local, y_local);
  return true;
}

void Finder::setFit(const std::vector<double>& x_fit,
                    const std::vector<double>& y_fit,
                    const std::vector<double>& y_background)
{
  if ((x_fit.size() != y_fit.size())
      || (x_fit.size() != y_background.size())
      || (x_fit.empty()))
    return;

  size_t l = find_index(x_fit.front());
  size_t r = find_index(x_fit.back());

  if ((r - l + 1) != x_fit.size())
    return;

  for (size_t i = 0; i < x_fit.size(); ++i)
  {
    y_fit_[l + i] = y_fit[i];
    y_background_[l + i] = y_background[i];
    double resid = y_[l + i] - y_fit[i];
    y_resid_[l + i] = resid;
    y_resid_on_background_[l + i] = y_background[i] + resid;
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
  fw_theoretical_nrg.clear();
  fw_theoretical_bin.clear();

//  DBG << "Kon calc " << x_.front() << "-"  << x_.back();

  /*if (settings_.cali_fwhm_.valid() && settings_.cali_nrg_.valid())
  {
    for (size_t i=0; i < x_.size(); ++i)
    {
      double nrg = settings_.bin_to_nrg(x_[i]);
      double fw = settings_.cali_fwhm_.transform(nrg);
      double L = settings_.nrg_to_bin(nrg - fw/2);
      double R = settings_.nrg_to_bin(nrg + fw/2);
      double wchan = R-L;

      if (!std::isfinite(fw))
      {
        DBG << "Kon " << x_[i] << "->" << nrg << "  fw=" << fw
            << " L=" << L << " R=" << R << " wchan=" << wchan;
      }

      fw_theoretical_nrg.push_back(nrg);
      fw_theoretical_bin.push_back(wchan);
    }
  }*/

  uint16_t width = settings_.KON_width;

  if (width < 2)
    width = 2;

  double sigma = settings_.KON_sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.KON_sigma_resid;
  }


//  DBG << "<Finder> width " << settings_.KON_width;

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

  y_kon.resize(y_resid_.size(), 0);
  y_convolution.resize(y_resid_.size(), 0);
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
  lefts.clear();
  rights.clear();

  if (prelim.empty())
    return;

  //find edges of contiguous peak areas
  lefts.push_back(prelim[0]);
  size_t prev = prelim[0];
  for (size_t i = 0; i < prelim.size(); ++i)
  {
    size_t current = prelim[i];
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
    filtered.push_back((rights[i] + lefts[i]) / 2);

  for (size_t i = 0; i < filtered.size(); ++i)
  {
    lefts[i] = left_edge(lefts[i]);
    rights[i] = right_edge(rights[i]);
//    DBG << "<Finder> Peak " << lefts[i] << "-"  << filtered[i] << "-"  << rights[i];
  }
}

double Finder::find_left(double chan) const
{
  if (x_.empty())
    return 0;

  //assume x is monotone increasing

  double sigma = settings_.KON_sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.KON_sigma_resid;
  }

  double edge_threshold = -0.5 * sigma;

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

  double sigma = settings_.KON_sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.KON_sigma_resid;
  }

  //assume x is monotone increasing

  double edge_threshold = -0.5 * sigma;

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
    double goal = x_[idx] - width * settings_.ROI_extend_peaks / 2;
    while ((idx > 0) && (x_[idx] > goal))
      idx--;
    return idx;
  }

  double sigma = settings_.KON_sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.KON_sigma_resid;
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
    double goal = x_[idx] + width * settings_.ROI_extend_peaks / 2;
    while ((idx < x_.size()) && (x_[idx] < goal))
      idx++;
    return idx;
  }

  double sigma = settings_.KON_sigma_spectrum;
  if (y_resid_ != y_)
  {
//    DBG << "<Finder> Using sigma resid";
    sigma = settings_.KON_sigma_resid;
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

double Finder::weight_true(size_t i) const
{
  return std::sqrt(y_[i]);
}

double Finder::weight_phillips_marlow(size_t i) const
{
  double k0 = y_[i];

  if (k0 >= 25)
    return std::sqrt(k0);
  else
  {
    k0 = 1.0;
    if ((i > 0) && ((i + 1) < y_.size()))
      k0 = y_[i - 1] + y_[i] + y_[i + 1] / 3.0;
    return std::max(std::sqrt(k0), 1.0);
  }
}

double Finder::weight_revay_student(size_t i) const
{
  double k0 = y_[i] + 1;
  return std::sqrt(k0);
}

}
