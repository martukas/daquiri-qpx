#include <optimizerBFGS/Region.h>
#include <optimizerBFGS/more_math.h>

#include <core/util/custom_logger.h>

namespace Hypermet
{

Region::Region(CSpectrum& spe, size_t from_channel, size_t to_channel)
    : spectrum(spe)
      , first_channel(std::min(from_channel, to_channel))
      , last_channel(std::max(from_channel, to_channel))
{
  // \todo validate region bounds

  background_slope_.to_fit = true;
  background_curve_.to_fit = true;

  default_peak_.width_.bound(0.8, 4.0);

  default_peak_.short_tail.amplitude.bound(0.02, 1.5);
  default_peak_.short_tail.amplitude.to_fit = true;
  default_peak_.short_tail.slope.bound(0.2, 0.5);
  default_peak_.short_tail.slope.to_fit = true;

  default_peak_.right_tail.amplitude.bound(0.01, 0.9);
  default_peak_.right_tail.amplitude.to_fit = true;
  default_peak_.right_tail.slope.bound(0.3, 1.5);
  default_peak_.right_tail.slope.to_fit = true;

  default_peak_.long_tail.amplitude.bound(0.0001, 0.15);
  default_peak_.long_tail.amplitude.to_fit = true;
  default_peak_.long_tail.slope.bound(2.5, 50);
  default_peak_.long_tail.slope.to_fit = true;

  default_peak_.step.amplitude.bound(0.000001, 0.05);
  default_peak_.step.amplitude.to_fit = true;

  //  SearchPeaks(); avoid this for the sake of batch fitting
}

//bool Region::step() const
//{
//  return step_enabled_;
//}
//
//void Region::step(bool enable)
//{
//  step_enabled_ = enable;
//  step_amplitude_.to_fit = enable;
//}

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
  try
  {
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
  catch (...)
  {
    ERR("Search Peaks failed!");
  }
}

void Region::add_peak(double position, double min, double max, double amplitude)
{
  // \todo this will be wrong in case of first peak
  int32_t base_index = 0;
  if (peaks_.size())
    base_index = peaks_.back().amplitude.x_index;

  peaks_.emplace_back();

  peaks_.back().amplitude.x_index = base_index + 2;
  peaks_.back().amplitude.val(amplitude);
  peaks_.back().amplitude.uncert_value = 0;

  peaks_.back().position.x_index = base_index + 3;
  peaks_.back().position.val(position);
  peaks_.back().position.min(min);
  peaks_.back().position.max(max);
  peaks_.back().position.uncert_value = 0;
}

void Region::remove_peak(size_t index)
{
  if (index >= peaks_.size())
    throw std::runtime_error("Can't delete the peak! (invalid index)");
  if (peaks_.size() == 1)
    throw std::runtime_error("Can't delete the only peak!");

  for (size_t i = index; i < peaks_.size() - 1; ++i)
  {
    peaks_[i].position = peaks_[i + 1].position;
    peaks_[i].amplitude = peaks_[i + 1].amplitude;
  } //i
  peaks_.resize(peaks_.size() - 1);
}

double Region::peak_area(size_t index) const
{
  return peaks_[index].area();
}

double Region::peak_area_unc(size_t index) const
{
  // \todo make this more rigorous
  return peaks_[index].area_uncert(chi_sq_normalized());
}

double Region::peak_area_eff(size_t index, const Calibration& cal)
{
  double eff{1.0};
  if (cal.efficiency.initialized())
    eff = cal.efficiency.val(peaks_[index].peak_energy(cal));
  return peak_area(index) / eff;
}

double Region::peak_area_eff_unc(size_t index, const Calibration& cal)
{
  double area = peak_area(index);
  double eff{0.0};
  double sigrel_eff{0.0};
  if (cal.efficiency.initialized())
  {
    auto energy = peaks_[index].peak_energy(cal);
    eff = cal.efficiency.val(energy);
    sigrel_eff = cal.efficiency.sigma_rel(energy);
  }
  return (square(std::sqrt(std::sqrt(area) / area)) + square(sigrel_eff)) *
      (area / eff) * std::max(1.0, chi_sq_normalized());
}

size_t Region::fit_var_count() const
{
  return current_fit.size();
}

void Region::map_fit()
{
  current_fit.clear();

  size_t unique_widths {0};
  size_t unique_short_tails {0};
  size_t unique_right_tails {0};
  size_t unique_long_tails {0};
  size_t unique_steps {0};
  for (auto& p : peaks_)
  {
    if (p.width_override)
      unique_widths++;
    if (p.short_tail.override)
      unique_short_tails++;
    if (p.right_tail.override)
      unique_right_tails++;
    if (p.long_tail.override)
      unique_long_tails++;
    if (p.step.override)
      unique_steps++;
  }

  int32_t i = 0;
  background_base_.x_index = i++;

  if (slope_enabled_)
    background_slope_.x_index = i++;
  else
    background_slope_.x_index = -1;

  if (curve_enabled_)
    background_curve_.x_index = i++;
  else
    background_curve_.x_index = -1;

  if (unique_widths < peaks_.size())
    default_peak_.width_.x_index = i++;
  else
    default_peak_.width_.x_index = -1;

  if (default_peak_.short_tail.enabled &&
      (unique_short_tails < peaks_.size()))
    default_peak_.short_tail.update_indices(i);

  if (default_peak_.right_tail.enabled &&
      (unique_right_tails < peaks_.size()))
    default_peak_.right_tail.update_indices(i);

  if (default_peak_.long_tail.enabled &&
      (unique_long_tails < peaks_.size()))
    default_peak_.long_tail.update_indices(i);

  if (default_peak_.step.enabled &&
      (unique_steps < peaks_.size()))
    default_peak_.step.update_indices(i);

  for (auto& p : peaks_)
    p.update_indices(i);

  current_fit.resize(static_cast<size_t>(i), 0.0);
  //fit_gradients.resize(vars);
  //ReDim ChisqGradient(FitVars - 1)
}

void Region::load_fit()
{
  background_base_.put(current_fit);

  background_slope_.put(current_fit);
  background_curve_.put(current_fit);

  default_peak_.put(current_fit);

  // \todo copy defaults
  for (auto& p : peaks_)
    p.put(current_fit);
}

void Region::save_fit()
{
  background_base_.get(current_fit);

  background_slope_.get(current_fit);
  background_curve_.get(current_fit);

  default_peak_.get(current_fit);

  for (auto& p : peaks_)
    p.get(current_fit);
}

void Region::save_fit_uncerts()
{
  save_fit();

  std::vector<double> diagonals;
  diagonals.reserve(current_fit.size());

  double df = degrees_of_freedom();
  for (size_t i = 0; i < current_fit.size(); ++i)
    diagonals.push_back(inv_hessian.coeff(i,i) * df);

  double chisq_norm = std::max(chi_sq_normalized(), 1.0) * 0.5;

  background_base_.get_uncert(current_fit, chisq_norm);

  background_slope_.get_uncert(current_fit, chisq_norm);
  background_curve_.get_uncert(current_fit, chisq_norm);

  default_peak_.get_uncerts(current_fit, chisq_norm);

  for (auto& p : peaks_)
    p.get_uncerts(current_fit, chisq_norm);
}

double Region::chi_sq_normalized() const
{
  return calc_chi_sq() / static_cast<double>(degrees_of_freedom());
}

size_t Region::degrees_of_freedom() const
{
  // \todo should it not be (last - first +1)?
  // \todo what if channel range is < fit_var_count?
  return ((last_channel - first_channel) - fit_var_count());
}

double Region::grad_chi_sq(std::vector<double>& gradients) const
{
  //Calculates the Chi-square and its gradient

  /*if(DiffType = 2)
  {
      Call dfunc2(reg, XVector, XGradient, Chisq)
      Exit Sub
  }

  if(DiffType = 3)
  {
      Call dfunc3(reg, XVector, XGradient, Chisq)
      Exit Sub
  }*/

  //Dim XGradient2(XGradient.GetLength(0) - 1) As Double, Chisq2 As Double
  //dfunc2(reg, XVector, XGradient2, Chisq2)

  // zero-out arrays
  gradients.assign(gradients.size(), 0.0);
  auto chan_gradients = gradients;

  double Chisq = 0;

  for (size_t pos = first_channel; pos <= last_channel; ++pos)
  {
    chan_gradients.assign(chan_gradients.size(), 0.0);

    //--- Poly Background ---
    double FTotal = background_base_.val();
    chan_gradients[background_base_.x_index] = background_base_.grad();
    if (slope_enabled_)
    {
      FTotal += background_slope_.val() * (pos - first_channel);
      chan_gradients[background_slope_.x_index] = (pos - first_channel);
    }

    if (curve_enabled_)
    {
      FTotal += background_curve_.val() * square(pos - first_channel);
      chan_gradients[background_curve_.x_index] = square(pos - first_channel);
    }

    for (auto& p : peaks_)
    {
      auto ret = p.eval_grad(pos, chan_gradients);
      FTotal += ret.gaussian + ret.step + ret.short_tail + ret.right_tail + ret.long_tail;
    } //Peak

    double t3 = -2.0 * (spectrum.channels[pos] - FTotal) / square(spectrum.weight(pos));
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((spectrum.channels[pos] - FTotal) / spectrum.weight(pos));
  }
  //Chisq /= df

  return Chisq;
}

double Region::calc_chi_sq_at(const std::vector<double>& fit) const
{
  //Calculates the Chi-square over a region
  double ChiSq = 0;

  for (size_t pos = first_channel; pos <= last_channel; ++pos)
  {
    // Background
    double FTotal = background_base_.val_at(fit[background_base_.x_index]);
    if (slope_enabled_)
      FTotal += background_slope_.val_at(fit[background_slope_.x_index]) * (pos - first_channel);
    if (curve_enabled_)
      FTotal += background_curve_.val_at(fit[background_curve_.x_index]) * square(pos - first_channel);

    for (auto& p : peaks_)
    {
      auto ret = p.eval_at(pos, fit);
      FTotal += ret.gaussian + ret.step + ret.short_tail + ret.right_tail + ret.long_tail;
    }
    ChiSq += square((spectrum.channels[pos] - FTotal) /
        spectrum.weight(pos));
  } //Channel

  return ChiSq;
}

double Region::grad_chi_sq_at(const std::vector<double>& fit,
    std::vector<double>& gradients) const
{
  //Calculates the Chi-square and its gradient

  /*if(DiffType = 2)
  {
      Call dfunc2(reg, XVector, XGradient, Chisq)
      Exit Sub
  }

  if(DiffType = 3)
  {
      Call dfunc3(reg, XVector, XGradient, Chisq)
      Exit Sub
  }*/

  //Dim XGradient2(XGradient.GetLength(0) - 1) As Double, Chisq2 As Double
  //dfunc2(reg, XVector, XGradient2, Chisq2)

  // zero-out arrays
  gradients.assign(gradients.size(), 0.0);
  auto chan_gradients = gradients;

  double Chisq = 0;

  for (size_t pos = first_channel; pos <= last_channel; ++pos)
  {
    chan_gradients.assign(chan_gradients.size(), 0.0);

    //--- Poly Background ---
    double FTotal = background_base_.val_at(fit[background_base_.x_index]);
    chan_gradients[background_base_.x_index] = background_base_.grad_at(fit[background_base_.x_index]);
    if (slope_enabled_)
    {
      FTotal += background_slope_.val_at(fit[background_slope_.x_index]) * (pos - first_channel);
      chan_gradients[background_slope_.x_index] = (pos - first_channel);
    }

    if (curve_enabled_)
    {
      FTotal += background_curve_.val_at(fit[background_curve_.x_index]) * square(pos - first_channel);
      chan_gradients[background_curve_.x_index] = square(pos - first_channel);
    }

    for (auto& p : peaks_)
    {
      auto ret = p.eval_grad_at(pos, fit, chan_gradients);
      FTotal += ret.gaussian + ret.step + ret.short_tail + ret.right_tail + ret.long_tail;
    } //Peak

    double t3 = -2.0 * (spectrum.channels[pos] - FTotal) / square(spectrum.weight(pos));
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((spectrum.channels[pos] - FTotal) / spectrum.weight(pos));
  }
  //Chisq /= df

  return Chisq;
}

}
