#include <core/fitting/roi.h>
#include <core/util/custom_logger.h>
#include <core/util/timer.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri {

Fit::Fit(const SUM4Edge &lb, const SUM4Edge &rb,
         const std::map<double, Peak> &peaks,
         const Finder &finder,
         std::string descr)
         : finder_(finder)
{
  settings_ = finder.settings_;
  //background_ = backg;
  LB_ = lb;
  RB_ = rb;
  peaks_ = peaks;

  description.description = descr;
  description.peaknum = peaks_.size();
  if (peaks_.size()) {
    description.rsq = chi_sq_normalized();
    // \todo use uncertan for these 2
    double tot_gross {0.0};
    double tot_back {0.0};
    for (auto &p : peaks_)
    {
      tot_gross += p.second.sum4().gross_area();
      tot_back  += p.second.sum4().background_area();
      p.second.hr_peak_.clear();
      p.second.hr_fullfit_.clear();
    }
    auto tot_net = tot_gross - tot_back;
    // \todo reenable this
    //description.sum4aggregate = tot_net.error();
  }

  // \todo only upon default creation!
  background_slope_.to_fit = true;
  background_curve_.to_fit = true;

  default_peak_.hypermet_.width_.bound(0.8, 4.0);

  default_peak_.hypermet_.short_tail.amplitude.bound(0.02, 1.5);
  default_peak_.hypermet_.short_tail.amplitude.to_fit = true;
  default_peak_.hypermet_.short_tail.slope.bound(0.2, 0.5);
  default_peak_.hypermet_.short_tail.slope.to_fit = true;

  default_peak_.hypermet_.right_tail.amplitude.bound(0.01, 0.9);
  default_peak_.hypermet_.right_tail.amplitude.to_fit = true;
  default_peak_.hypermet_.right_tail.slope.bound(0.3, 1.5);
  default_peak_.hypermet_.right_tail.slope.to_fit = true;

  default_peak_.hypermet_.long_tail.amplitude.bound(0.0001, 0.15);
  default_peak_.hypermet_.long_tail.amplitude.to_fit = true;
  default_peak_.hypermet_.long_tail.slope.bound(2.5, 50);
  default_peak_.hypermet_.long_tail.slope.to_fit = true;

  default_peak_.hypermet_.step.amplitude.bound(0.000001, 0.05);
  default_peak_.hypermet_.step.amplitude.to_fit = true;
}


//bool Fit::step() const
//{
//  return step_enabled_;
//}
//
//void Fit::step(bool enable)
//{
//  step_enabled_ = enable;
//  step_amplitude_.to_fit = enable;
//}


//void Fit::add_peak(double position, double min, double max, double amplitude)
//{
//  peaks_.emplace_back();
//
//  peaks_.back().amplitude.val(amplitude);
//  peaks_.back().amplitude.uncert_value = 0;
//
//  peaks_.back().position.val(position);
//  peaks_.back().position.min(min);
//  peaks_.back().position.max(max);
//  peaks_.back().position.uncert_value = 0;
//}

double Fit::peak_area(size_t index) const
{
  return peaks_[index].area();
}

double Fit::peak_area_unc(size_t index) const
{
  // \todo make this more rigorous
  return peaks_[index].area_uncert(chi_sq_normalized());
}

double Fit::peak_area_eff(size_t index, const Calibration& cal)
{
  double eff{1.0};
  if (cal.efficiency.initialized())
    eff = cal.efficiency.val(peaks_[index].hypermet().peak_energy(cal));
  return peak_area(index) / eff;
}

double Fit::peak_area_eff_unc(size_t index, const Calibration& cal)
{
  double area = peak_area(index);
  double eff{0.0};
  double sigrel_eff{0.0};
  if (cal.efficiency.initialized())
  {
    auto energy = peaks_[index].hypermet().peak_energy(cal);
    eff = cal.efficiency.val(energy);
    sigrel_eff = cal.efficiency.sigma_rel(energy);
  }
  return (square(std::sqrt(std::sqrt(area) / area)) + square(sigrel_eff)) *
      (area / eff) * std::max(1.0, chi_sq_normalized());
}

void Fit::map_fit()
{
  size_t unique_widths{0};
  size_t unique_short_tails{0};
  size_t unique_right_tails{0};
  size_t unique_long_tails{0};
  size_t unique_steps{0};
  for (auto& p : peaks_)
  {
    if (p.second.hypermet_.width_override)
      unique_widths++;
    if (p.second.hypermet_.short_tail.override)
      unique_short_tails++;
    if (p.second.hypermet_.right_tail.override)
      unique_right_tails++;
    if (p.second.hypermet_.long_tail.override)
      unique_long_tails++;
    if (p.second.hypermet_.step.override)
      unique_steps++;
  }

  var_count_ = 0;
  background_base_.x_index = var_count_++;

  if (slope_enabled_)
    background_slope_.x_index = var_count_++;
  else
    background_slope_.x_index = -1;

  if (curve_enabled_)
    background_curve_.x_index = var_count_++;
  else
    background_curve_.x_index = -1;

  if (unique_widths < peaks_.size())
    default_peak_.hypermet_.width_.x_index = var_count_++;
  else
    default_peak_.hypermet_.width_.x_index = -1;

  if (default_peak_.hypermet_.short_tail.enabled &&
      (unique_short_tails < peaks_.size()))
    default_peak_.hypermet_.short_tail.update_indices(var_count_);

  if (default_peak_.hypermet_.right_tail.enabled &&
      (unique_right_tails < peaks_.size()))
    default_peak_.hypermet_.right_tail.update_indices(var_count_);

  if (default_peak_.hypermet_.long_tail.enabled &&
      (unique_long_tails < peaks_.size()))
    default_peak_.hypermet_.long_tail.update_indices(var_count_);

  if (default_peak_.hypermet_.step.enabled &&
      (unique_steps < peaks_.size()))
    default_peak_.hypermet_.step.update_indices(var_count_);

  for (auto& p : peaks_)
    p.second.hypermet_.update_indices(var_count_);
}

size_t Fit::fit_var_count() const
{
  return static_cast<size_t>(var_count_);
}

std::vector<double> Fit::variables() const
{
  std::vector<double> ret;
  ret.resize(fit_var_count(), 0.0);

  background_base_.put(ret);

  background_slope_.put(ret);
  background_curve_.put(ret);

  default_peak_.hypermet_.put(ret);

  // \todo copy defaults
  for (auto& p : peaks_)
    p.second.hypermet_.put(ret);

  return ret;
}

void Fit::save_fit(const std::vector<double>& variables)
{
  background_base_.get(variables);

  background_slope_.get(variables);
  background_curve_.get(variables);

  default_peak_.hypermet_.get(variables);

  for (auto& p : peaks_)
    p.second.hypermet_.get(variables);
}

void Fit::save_fit_uncerts(const FitResult& result)
{
  save_fit(result.variables);

  std::vector<double> diagonals;
  diagonals.reserve(result.variables.size());

  double df = degrees_of_freedom();
  for (size_t i = 0; i < result.variables.size(); ++i)
    diagonals.push_back(result.inv_hessian.coeff(i, i) * df);

  double chisq_norm = std::max(chi_sq_normalized(), 1.0) * 0.5;

  background_base_.get_uncert(result.variables, chisq_norm);

  background_slope_.get_uncert(result.variables, chisq_norm);
  background_curve_.get_uncert(result.variables, chisq_norm);

  default_peak_.hypermet_.get_uncerts(result.variables, chisq_norm);

  for (auto& p : peaks_)
    p.second.hypermet_.get_uncerts(result.variables, chisq_norm);
}

double Fit::chi_sq_normalized() const
{
  return chi_sq() / degrees_of_freedom();
}

double Fit::degrees_of_freedom() const
{
  // \todo what if channel range is < fit_var_count?
  return (finder_.x_.size() - fit_var_count());
}

double Fit::chi_sq() const
{
  //Calculates the Chi-square over a region
  double ChiSq = 0;

  for (size_t pos = 0; pos < finder_.x_.size(); ++pos)
  {
    // Background
    double FTotal = background_base_.val();
    if (slope_enabled_)
      FTotal += background_slope_.val() * (pos - finder_.x_.front());
    if (curve_enabled_)
      FTotal += background_curve_.val() * square(pos - finder_.x_.front());

    for (auto& p : peaks_)
    {
      auto ret = p.second.hypermet_.eval(pos);
      FTotal += ret.gaussian + ret.step + ret.short_tail + ret.right_tail + ret.long_tail;
    }
    ChiSq += square((finder_.y_[pos] - FTotal) / finder_.y_weight_true[pos]);
  } //Channel

  return ChiSq;
}


double Fit::grad_chi_sq(std::vector<double>& gradients) const
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

  for (size_t pos = 0; pos < finder_.x_.size(); ++pos)
  {
    chan_gradients.assign(chan_gradients.size(), 0.0);

    //--- Poly Background ---
    double FTotal = background_base_.val();
    chan_gradients[background_base_.x_index] = background_base_.grad();
    if (slope_enabled_)
    {
      FTotal += background_slope_.val() * (pos - finder_.x_.front());
      chan_gradients[background_slope_.x_index] = (pos - finder_.x_.front());
    }

    if (curve_enabled_)
    {
      FTotal += background_curve_.val() * square(pos - finder_.x_.front());
      chan_gradients[background_curve_.x_index] = square(pos - finder_.x_.front());
    }

    for (auto& p : peaks_)
    {
      auto ret = p.second.hypermet_.eval_grad(pos, chan_gradients);
      FTotal += ret.gaussian + ret.step + ret.short_tail + ret.right_tail + ret.long_tail;
    } //Peak

    double t3 = -2.0 * (finder_.y_[pos] - FTotal) / square(finder_.y_weight_true[pos]);
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((finder_.y_[pos] - FTotal) / finder_.y_weight_true[pos]);
  }
  //Chisq /= df

  return Chisq;
}

double Fit::chi_sq(const std::vector<double>& fit) const
{
  //Calculates the Chi-square over a region
  double ChiSq = 0;

  for (size_t pos = 0; pos < finder_.x_.size(); ++pos)
  {
    // Background
    double FTotal = background_base_.val_at(fit[background_base_.x_index]);
    if (slope_enabled_)
      FTotal += background_slope_.val_at(fit[background_slope_.x_index]) * (pos - finder_.x_.front());
    if (curve_enabled_)
      FTotal += background_curve_.val_at(fit[background_curve_.x_index]) * square(pos - finder_.x_.front());

    for (auto& p : peaks_)
    {
      auto ret = p.second.hypermet_.eval_at(pos, fit);
      FTotal += ret.gaussian + ret.step + ret.short_tail + ret.right_tail + ret.long_tail;
    }
    ChiSq += square((finder_.y_[pos] - FTotal) / finder_.y_weight_phillips_marlow[pos]);
  } //Channel

  return ChiSq;
}

double Fit::grad_chi_sq(const std::vector<double>& fit,
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

  for (size_t pos = 0; pos < finder_.x_.size(); ++pos)
  {
    chan_gradients.assign(chan_gradients.size(), 0.0);

    //--- Poly Background ---
    double FTotal = background_base_.val_at(fit[background_base_.x_index]);
    chan_gradients[background_base_.x_index] = background_base_.grad_at(fit[background_base_.x_index]);
    if (slope_enabled_)
    {
      FTotal += background_slope_.val_at(fit[background_slope_.x_index]) * (pos - finder_.x_.front());
      chan_gradients[background_slope_.x_index] = (pos - finder_.x_.front());
    }

    if (curve_enabled_)
    {
      FTotal += background_curve_.val_at(fit[background_curve_.x_index]) * square(pos - finder_.x_.front());
      chan_gradients[background_curve_.x_index] = square(pos - finder_.x_.front());
    }

    for (auto& p : peaks_)
    {
      auto ret = p.second.hypermet_.eval_grad_at(pos, fit, chan_gradients);
      FTotal += ret.gaussian + ret.step + ret.short_tail + ret.right_tail + ret.long_tail;
    } //Peak

    double t3 = -2.0 * (finder_.y_[pos] - FTotal) / square(finder_.y_weight_phillips_marlow[pos]);
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((finder_.y_[pos] - FTotal) / finder_.y_weight_phillips_marlow[pos]);
  }
  //Chisq /= df

  return Chisq;
}



ROI::ROI(const Finder &parentfinder, double min, double max)
{
  finder_.settings_ = parentfinder.settings_;
  set_data(parentfinder, min, max);
}

double ROI::ID() const
{
  return left_bin();
}

double ROI::left_bin() const
{
  if (finder_.x_.empty())
    return -1;
  else
    return finder_.x_.front();
}

double ROI::right_bin() const
{
  if (finder_.x_.empty())
    return -1;
  else
    return finder_.x_.back();
}

double ROI::left_nrg() const
{
  if (hr_x_nrg.empty())
    return std::numeric_limits<double>::quiet_NaN();
  else
    return hr_x_nrg.front();
}

double ROI::right_nrg() const
{
  if (hr_x_nrg.empty())
    return std::numeric_limits<double>::quiet_NaN();
  else
    return hr_x_nrg.back();
}

double ROI::width() const
{
  if (finder_.x_.empty())
    return 0;
  else
    return right_bin() - left_bin() + 1;
}

void ROI::set_data(const Finder &parentfinder, double l, double r)
{
  if (!finder_.cloneRange(parentfinder, l, r))
  {
    finder_.clear();
    return;
  }

  hr_x.clear();
  hr_background.clear();
  hr_back_steps.clear();
  hr_fullfit.clear();
  init_edges();
  init_background();
  render();
}

bool ROI::refit(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  if (peaks_.empty())
    return auto_fit(optimizer, interruptor);

  if (!rebuild(optimizer, interruptor))
    return false;

  save_current_fit("Refit");
  return true;
}


bool ROI::auto_fit(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  peaks_.clear();
  finder_.y_resid_ = finder_.y_;
  finder_.find_peaks();  //assumes default params!!!

  if (finder_.filtered.empty())
    return false;

  if ((LB_.width() == 0) || (RB_.width() == 0))
  {
    init_edges();
    init_background();
  }

  if (!finder_.settings_.sum4_only)
  {
    std::vector<double> y_nobkg = remove_background();

    for (size_t i=0; i < finder_.filtered.size(); ++i)
    {
      std::vector<double> x_pk = std::vector<double>(finder_.x_.begin() + finder_.lefts[i], finder_.x_.begin() + finder_.rights[i] + 1);
      std::vector<double> y_pk = std::vector<double>(y_nobkg.begin() + finder_.lefts[i], y_nobkg.begin() + finder_.rights[i] + 1);

      Gaussian gaussian;
      optimizer->fit(gaussian, x_pk, y_pk);

      if (
          std::isfinite(gaussian.height().value()) && (gaussian.height().value() > 0) &&
              std::isfinite(gaussian.hwhm().value()) && (gaussian.hwhm().value() > 0) &&
          (finder_.x_[finder().lefts[i]] < gaussian.center().value()) &&
          (gaussian.center().value() < finder_.x_[finder_.rights[i]])
          )
      {
        Peak fitted(Hypermet(gaussian, finder_.settings_), SUM4(), finder_.settings_);
        peaks_[fitted.center()] = fitted;
      }
    }
    if (peaks_.empty())
      finder_.settings_.sum4_only = true;
  }

  if (!rebuild(optimizer, interruptor))
    return false;

  save_current_fit("Autofit");

  if (finder_.settings_.resid_auto)
    iterative_fit(optimizer, interruptor);

  return true;
}

void ROI::iterative_fit(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  if (!finder_.settings_.cali_fwhm_.valid() || peaks_.empty())
    return;

  double prev_rsq = peaks_.begin()->second.hypermet().chi2();

  for (int i=0; i < finder_.settings_.resid_max_iterations; ++i)
  {
    ROI new_fit = *this;

    DBG("Attempting add from resid with {} peaks", peaks_.size());

    if (!new_fit.add_from_resid(optimizer, interruptor, -1)) {
      //      DBG << "    failed add from resid";
      break;
    }
    double new_rsq = new_fit.peaks_.begin()->second.hypermet().chi2();
    double improvement = (prev_rsq - new_rsq) / prev_rsq * 100;
    DBG("X2 new={} previous={} improved by {}", new_rsq, prev_rsq, improvement);

    if ((new_rsq > prev_rsq) || std::isnan(new_rsq)) {
      DBG(" not improved. reject new refit");
      break;
    }

    new_fit.save_current_fit("Iterative " + std::to_string(new_fit.peaks().size()));
    prev_rsq = new_rsq;
    *this = new_fit;

    if (interruptor.load())
      break;
  }
}

bool ROI::add_from_resid(BFGS& optimizer, std::atomic<bool>& interruptor, int32_t centroid_hint)
{
  if (finder_.filtered.empty())
    return false;

  int64_t target_peak = -1;
  if (centroid_hint == -1) {
    double biggest = 0;
    for (size_t j=0; j < finder_.filtered.size(); ++j)
    {
      std::vector<double> x_pk = std::vector<double>(finder_.x_.begin() + finder_.lefts[j],
                                                     finder_.x_.begin() + finder_.rights[j] + 1);
      std::vector<double> y_pk = std::vector<double>(finder_.y_resid_.begin() + finder_.lefts[j],
                                                     finder_.y_resid_.begin() + finder_.rights[j] + 1);
      Gaussian gaussian;
      optimizer->fit(gaussian, x_pk, y_pk);

      bool too_close = false;

      double lateral_slack = finder_.settings_.resid_too_close * gaussian.hwhm().value() * 2;
      for (auto &p : peaks_)
      {
        if ((p.second.center() > (gaussian.center().value() - lateral_slack))
            && (p.second.center() < (gaussian.center().value() + lateral_slack)))
          too_close = true;
      }

      //      if (too_close)
      //        DBG << "Too close at " << settings_.cali_nrg_.transform(gaussian.center_, settings_.bits_);

      if ( !too_close &&
          std::isfinite(gaussian.height().value()) && (gaussian.height().value() > 0) &&
          std::isfinite(gaussian.hwhm().value()) && (gaussian.hwhm().value() > 0) &&
           (finder_.x_[finder_.lefts[j]] < gaussian.center().value()) &&
           (gaussian.center().value() < finder_.x_[finder_.rights[j]]) &&
           (gaussian.height().value() > finder_.settings_.resid_min_amplitude) &&
           (gaussian.area() > biggest)
           )
      {
        target_peak = j;
        biggest = gaussian.area();
      }
    }
    //    DBG << "    biggest potential add at " << finder_.x_[finder_.filtered[target_peak]] << " with area=" << biggest;
  } else {

    //THIS NEVER HAPPENS
    double diff = std::abs(finder_.x_[finder_.filtered[target_peak]] - centroid_hint);
    for (size_t j=0; j < finder_.filtered.size(); ++j)
      if (std::abs(finder_.x_[finder_.filtered[j]] - centroid_hint) < diff) {
        target_peak = j;
        diff = std::abs(finder_.x_[finder_.filtered[j]] - centroid_hint);
      }
  }

  if (target_peak == -1)
    return false;

  std::vector<double> x_pk = std::vector<double>(finder_.x_.begin() + finder_.lefts[target_peak],
                                                 finder_.x_.begin() + finder_.rights[target_peak] + 1);
  std::vector<double> y_pk = std::vector<double>(finder_.y_resid_.begin() + finder_.lefts[target_peak],
                                                 finder_.y_resid_.begin() + finder_.rights[target_peak] + 1);
  Gaussian gaussian;
  optimizer->fit(gaussian, x_pk, y_pk);

  if (
      std::isfinite(gaussian.height().value()) && (gaussian.height().value() > 0) &&
          std::isfinite(gaussian.hwhm().value()) && (gaussian.hwhm().value() > 0) &&
      (finder_.x_[finder_.lefts[target_peak]] < gaussian.center().value()) &&
      (gaussian.center().value() < finder_.x_[finder_.rights[target_peak]])
      )
  {
    Peak fitted(Hypermet(gaussian, finder_.settings_), SUM4(), finder_.settings_);
    peaks_[fitted.center()] = fitted;

    rebuild(optimizer, interruptor);
    return true;
  }
  else
    return false;
}


bool ROI::contains(double peakID) const
{
  return (peaks_.count(peakID) > 0);
}

Peak ROI::peak(double peakID) const
{
  if (contains(peakID))
    return peaks_.at(peakID);
  else
    return Peak();
}


bool ROI::overlaps(double bin) const
{
  if (!width())
    return false;
  return ((bin >= left_bin()) && (bin <= right_bin()));
}

bool ROI::overlaps(double Lbin, double Rbin) const
{
  if (finder_.x_.empty())
    return false;
  if (overlaps(Lbin) || overlaps(Rbin))
    return true;
  if ((Lbin <= left_bin()) && (Rbin >= right_bin()))
    return true;
  return false;
}

bool ROI::overlaps(const ROI& other) const
{
  if (!other.width())
    return false;
  return overlaps(other.left_bin(), other.right_bin());
}

size_t ROI::peak_count() const
{
  return peaks_.size();
}

const std::map<double, Peak> &ROI::peaks() const
{
  return peaks_;
}


bool ROI::adjust_sum4(double &peakID, double left, double right)
{
  if (!contains(peakID))
    return false;

  Peak pk = peaks_.at(peakID);
  SUM4 new_sum4(left, right, finder_, LB_, RB_);
  pk = Peak(pk.hypermet(), new_sum4, finder_.settings_);
  remove_peak(peakID);
  peakID = pk.center();
  peaks_[peakID] = pk;
  render();
  save_current_fit("SUM4 adjusted on " + std::to_string(pk.energy()));
  return true;
}


bool ROI::replace_hypermet(double &peakID, Peak hyp)
{
  if (!contains(peakID))
    return false;

  Peak pk = peaks_.at(peakID);
  pk = Peak(hyp, pk.sum4(), finder_.settings_);
  remove_peak(peakID);
  peakID = pk.center();
  peaks_[peakID] = pk;
  //set rsq to 0 for all peaks?

  render();
  save_current_fit("Hypermet adjusted on " + std::to_string(pk.energy()));
  return true;
}

bool ROI::override_energy(double peakID, double energy)
{
  if (!contains(peakID))
    return false;

   peaks_[peakID].override_energy(energy);

   render();
   save_current_fit("Peak energy override " + std::to_string(peaks_.at(peakID).center())
                    + "->" + std::to_string(peaks_.at(peakID).energy()));
   return true;
}


bool ROI::add_peak(const Finder &parentfinder,
                   double left, double right,
                   BFGS& optimizer,
                   std::atomic<bool>& interruptor)
{
  uint16_t center_prelim = (left+right) * 0.5; //assume down the middle

  if (overlaps(left) && overlaps(right)) {
    ROI new_fit = *this;

    if (!finder_.settings_.sum4_only && new_fit.add_from_resid(optimizer, interruptor, center_prelim))
    {
      *this = new_fit;
      save_current_fit("Added from residuals");
      return true;
    }
    else
    {
      Peak fitted(Hypermet(), SUM4(left, right, finder_, LB_, RB_), finder_.settings_);
      peaks_[fitted.center()] = fitted;
      render();
      save_current_fit("Manually added " + std::to_string(fitted.energy()));
      return true;
    }
  }
  else if (width()) //changing region bounds
  {
    if (!finder_.x_.empty())
    {
      left  = std::min(left, left_bin());
      right = std::max(right, right_bin());
    }
    if (!finder_.cloneRange(parentfinder, left, right))
      return false;

    init_edges();
    init_background();
    finder_.y_resid_ = remove_background();
    render();
    finder_.find_peaks();  //assumes default params!!!

    ROI new_fit = *this;
    if (finder_.settings_.sum4_only)
    {
      Peak fitted(Hypermet(), SUM4(left, right, finder_, LB_, RB_), finder_.settings_);
      peaks_[fitted.center()] = fitted;
      render();
      save_current_fit("Manually added " + std::to_string(fitted.energy()));
      return true;
    }
    else if (new_fit.add_from_resid(optimizer, interruptor, finder_.find_index(center_prelim)))
    {
      *this = new_fit;
      save_current_fit("Added from residuals");
      return true;
    }
    else
      return auto_fit(optimizer, interruptor);
  }

  DBG("<ROI> cannot add to empty ROI");
  return false;
}

bool ROI::remove_peaks(const std::set<double> &pks, BFGS& optimizer, std::atomic<bool>& interruptor)
{
  bool found = false;
  for (auto &q : pks)
    if (remove_peak(q))
      found = true;

  if (!found)
    return false;

  if (peaks_.size() && !rebuild(optimizer, interruptor))
    return false;

  render();
  save_current_fit("Peaks removed");
  return true;
}

bool ROI::remove_peak(double bin)
{
  if (peaks_.count(bin)) {
    peaks_.erase(bin);
    return true;
  }
  return false;
}

bool ROI::override_settings(const FitSettings &fs, std::atomic<bool>& interruptor)
{
  finder_.settings_ = fs;
  finder_.settings_.overriden = true; //do this in fitter if different?
  save_current_fit("Fit settings overriden");

  //propagate to peaks

  //render if calibs changed?
  return true;
}

void ROI::save_current_fit(std::string description)
{
  Fit thisfit(LB_, RB_, background_,
              peaks_, finder_, description);
  fits_.push_back(thisfit);
  current_fit_ = fits_.size() - 1;
}

bool ROI::rebuild(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  hr_x.clear();
  hr_x_nrg.clear();
  hr_background.clear();
  hr_back_steps.clear();
  hr_fullfit.clear();

  bool hypermet_fit = false;
  for (auto &q : peaks_)
    if (!q.second.hypermet().gaussian_only())
    {
      hypermet_fit = true;
      break;
    }

  bool success = false;
  if (hypermet_fit)
    success = rebuild_as_hypermet(optimizer, interruptor);
  else
    success = rebuild_as_gaussian(optimizer, interruptor);

  if (!success)
    return false;

  render();
  return true;
}

bool ROI::rebuild_as_hypermet(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  Timer timer(true);

  std::map<double, Peak> new_peaks;

  std::vector<Hypermet> old_hype;
  for (auto &q : peaks_) {
    if (q.second.hypermet().height().value())
      old_hype.push_back(q.second.hypermet());
    else if (q.second.sum4().peak_width())
    {
      Peak s4only(Hypermet(),
                  SUM4(q.second.sum4().left(), q.second.sum4().right(),
                       finder_, LB_, RB_),
                  finder_.settings_);
      new_peaks[s4only.center()] = s4only;
    }
  }

  if (old_hype.empty())
    return false;

//  std::vector<Hypermet> hype = fit_multi(finder_.x_, finder_.y_,
//                                                   old_hype, background_,
//                                                   finder_.settings_);

  std::vector<Hypermet> hype = optimizer->fit_multiplet(finder_.x_, finder_.y_,
                                                        old_hype, background_,
                                                        finder_.settings_);

  for (size_t i=0; i < hype.size(); ++i) {
    double edge =  hype[i].width().value() * sqrt(log(2)) * 3; //use const from settings
    double left = hype[i].center().value() - edge;
    double right = hype[i].center().value() + edge;
    Peak one(hype[i], SUM4(left, right, finder_, LB_, RB_), finder_.settings_);
    new_peaks[one.center()] = one;
  }

  peaks_ = new_peaks;
  return true;
}

bool ROI::rebuild_as_gaussian(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  Timer timer(true);

  std::map<double, Peak> new_peaks;

  std::vector<Gaussian> old_gauss;
  for (auto &q : peaks_)
  {
    if (q.second.hypermet().height().value())
      old_gauss.push_back(q.second.hypermet().gaussian());
    else if (q.second.sum4().peak_width())
    {
      Peak s4only(Hypermet(), SUM4(q.second.sum4().left(),
                                   q.second.sum4().right(),
                                   finder_, LB_, RB_),
                  finder_.settings_);
      //      q.second.construct(settings_);
      new_peaks[s4only.center()] = s4only;
    }
  }

  if (old_gauss.empty())
    return false;

  std::vector<Gaussian> gauss = optimizer->fit_multiplet(finder_.x_, finder_.y_,
                                                         old_gauss, background_,
                                                         finder_.settings_);

  for (size_t i=0; i < gauss.size(); ++i)
  {
    double edge =  gauss[i].hwhm().value() * 3; //use const from settings
    double left = gauss[i].center().value() - edge;
    double right = gauss[i].center().value() + edge;
    Peak one(Hypermet(gauss[i], finder_.settings_),
             SUM4(left, right, finder_, LB_, RB_),
             finder_.settings_);
    new_peaks[one.center()] = one;
  }

  peaks_ = new_peaks;
  return true;
}


void ROI::render()
{
  hr_x.clear();
  hr_background.clear();
  hr_back_steps.clear();
  hr_fullfit.clear();
  Polynomial sum4back = SUM4::sum4_background(LB_, RB_, finder_);

  for (double i = 0; i < finder_.x_.size(); i += 0.1)
  {
    hr_x.push_back(finder_.x_[0] + i);
    hr_fullfit.push_back(finder_.y_[std::floor(i)]);
  }
  hr_background = background_.eval(hr_x);
  hr_sum4_background_ = sum4back.eval(hr_x);
  hr_x_nrg = finder_.settings_.cali_nrg_.transform(hr_x);

  std::vector<double> lowres_backsteps = sum4back.eval(finder_.x_);
  std::vector<double> lowres_fullfit   = sum4back.eval(finder_.x_);

  for (auto &p : peaks_)
    p.second.hr_fullfit_ = p.second.hr_peak_ = hr_fullfit;

  if (!finder_.settings_.sum4_only) {
    hr_fullfit    = hr_background;
    hr_back_steps = hr_background;
    lowres_backsteps = background_.eval(finder_.x_);
    lowres_fullfit = background_.eval(finder_.x_);
    for (auto &p : peaks_) {
      for (int32_t j = 0; j < static_cast<int32_t>(hr_x.size()); ++j) {
        double step = p.second.hypermet().eval_step_tail(hr_x[j]);
        hr_back_steps[j] += step;
        hr_fullfit[j]    += step + p.second.hypermet().eval_peak(hr_x[j]);
      }

      for (int32_t j = 0; j < static_cast<int32_t>(finder_.x_.size()); ++j) {
        double step = p.second.hypermet().eval_step_tail(finder_.x_[j]);
        lowres_backsteps[j] += step;
        lowres_fullfit[j]   += step + p.second.hypermet().eval_peak(finder_.x_[j]);
      }
    }

    for (auto &p : peaks_) {
      p.second.hr_fullfit_ = hr_back_steps;
      p.second.hr_peak_.resize(hr_back_steps.size());
      for (int32_t j = 0; j < static_cast<int32_t>(hr_x.size()); ++j) {
        p.second.hr_peak_[j]     = p.second.hypermet().eval_peak(hr_x[j]);
        p.second.hr_fullfit_[j] += p.second.hr_peak_[j];
      }
    }
  }

  finder_.reset();
  finder_.setFit(finder_.x_, lowres_fullfit, lowres_backsteps);
}

std::vector<double> ROI::remove_background()
{
  std::vector<double> y_background = background_.eval(finder_.x_);

  std::vector<double> y_nobkg(finder_.x_.size());
  for (int32_t i = 0; i < static_cast<int32_t>(finder_.y_.size()); ++i)
    y_nobkg[i] = finder_.y_[i] - y_background[i];

  return y_nobkg;
}

bool ROI::adjust_LB(const Finder &parentfinder, double left, double right,
                     BFGS& optimizer, std::atomic<bool>& interruptor)
{
  size_t Lidx = parentfinder.find_index(left);
  size_t Ridx = parentfinder.find_index(right);
  if (Lidx >= Ridx)
    return false;

  SUM4Edge edge(parentfinder.x_, parentfinder.y_, Lidx, Ridx);
  if (!edge.width() || (edge.right() >= RB_.left()))
    return false;

  if ((edge.left() != left_bin()) && !finder_.cloneRange(parentfinder, left, right_bin()))
    return false;

  LB_ = edge;
  init_background();
  cull_peaks();
  render();
  rebuild(optimizer, interruptor);
  save_current_fit("Left baseline adjusted");
  return true;
}

bool ROI::adjust_RB(const Finder &parentfinder, double left, double right,
                    BFGS& optimizer, std::atomic<bool>& interruptor) {
  size_t Lidx = parentfinder.find_index(left);
  size_t Ridx = parentfinder.find_index(right);
  if (Lidx >= Ridx)
    return false;

  SUM4Edge edge(parentfinder.x_, parentfinder.y_, Lidx, Ridx);
  if (!edge.width() || (edge.left() <= LB_.right()))
    return false;

  if ((edge.right() != right_bin()) && !finder_.cloneRange(parentfinder, left_bin(), right))
    return false;

  RB_ = edge;
  init_background();
  cull_peaks();
  render();
  rebuild(optimizer, interruptor);
  save_current_fit("Right baseline adjusted");
  return true;
}

void ROI::init_edges()
{
  init_LB();
  init_RB();
}

void ROI::init_LB()
{
  int32_t LBend = 0;
  uint16_t samples = finder_.settings_.background_edge_samples;

  if ((samples > 0) && (finder_.y_.size() > samples*3))
    LBend += samples;

  LB_ = SUM4Edge(finder_.x_, finder_.y_, 0, LBend);
}

void ROI::init_RB()
{
  int32_t RBstart = finder_.y_.size() - 1;
  uint16_t samples = finder_.settings_.background_edge_samples;
  if ((samples > 0) && (finder_.y_.size() > samples*3))
    RBstart -= samples;

  RB_ = SUM4Edge(finder_.x_, finder_.y_, RBstart, finder_.y_.size() - 1);
}

void ROI::init_background()
{
  if (finder_.x_.empty())
    return;

  background_ = Polynomial();
  auto xoffset = background_.x_offset();

  //by default, linear
  double run = RB_.left() - LB_.right();
  xoffset.constrain(LB_.left(), LB_.left());

  double minslope = 0, maxslope = 0;
  double ymin, ymax, yav;
  if (LB_.average() < RB_.average())
  {
    run = RB_.right() - LB_.right();
    xoffset.constrain(LB_.right(), LB_.right());
    minslope = (RB_.min() - LB_.max()) / (RB_.right() - LB_.left());
    maxslope = (RB_.max() - LB_.min()) / (RB_.left() - LB_.right());
    ymin = LB_.min();
    ymax = RB_.max();
    yav = LB_.average();
  }

  if (RB_.average() < LB_.average())
  {
    run = RB_.left() - LB_.left();
    xoffset.constrain(LB_.left(), LB_.left());
    minslope = (RB_.min() - LB_.max()) / (RB_.left() - LB_.right());
    maxslope = (RB_.max() - LB_.min()) / (RB_.right() - LB_.left());
    ymin = RB_.min();
    ymax = LB_.max();
    yav = RB_.average();
  }

  double maxcurve = (run*run - std::min(LB_.min(), RB_.min()))
                    / std::max(LB_.max(), RB_.max());

  double slope = (RB_.average() - LB_.average()) / run ;

  background_.x_offset(xoffset);

  background_.set_coeff(0, {ymin, ymax, yav});

  background_.set_coeff(1, {0.5 * minslope, 2 * maxslope, slope});

  background_.set_coeff(2, {-maxcurve, maxcurve, 0});
}

size_t ROI::current_fit() const
{
  return current_fit_;
}

size_t ROI::history_size() const
{
  return fits_.size();
}

std::vector<FitDescription> ROI::history() const
{
  std::vector<FitDescription> ret;
  for (auto &f : fits_)
    ret.push_back(f.description);
  return ret;
}


bool ROI::rollback(const Finder &parent_finder, size_t i)
{
  if (i >= fits_.size())
    return false;

  finder_.settings_ = fits_[i].settings_;
  set_data(parent_finder, fits_[i].LB_.left(), fits_[i].RB_.right());
  background_ = fits_[i].background_;
  LB_ = fits_[i].LB_;
  RB_ = fits_[i].RB_;
  peaks_ = fits_[i].peaks_;
  render();

  current_fit_ = i;

  return true;
}

void ROI::cull_peaks()
{
  std::map<double, Peak> peaks;
  for (auto &p : peaks_) {
    if ((p.first > LB_.right()) &&
        (p.first < RB_.left()))
      peaks[p.first] = p.second;
  }
  peaks_ = peaks;
}

nlohmann::json ROI::to_json(const Finder &parent_finder) const
{
  nlohmann::json j;

  if (fits_.empty())
    return j;

  j["current_fit"] = current_fit_;

  ROI temp(*this);

  for (size_t i=0; i < temp.fits_.size(); ++i)
  {
    nlohmann::json jj;

    jj["description"] = temp.fits_[i].description.description;
    temp.rollback(parent_finder, i);

    if (finder_.settings_.overriden)
      jj["settings"] = finder_.settings_;

    jj["background_left"] = temp.LB();
    jj["background_right"] = temp.RB();
    jj["background_poly"] = temp.background_;

    for (auto &p : temp.peaks_)
      jj["peaks"].push_back(p.second);

    j["fits"].push_back(jj);
  }

  return j;
}

ROI::ROI(const nlohmann::json& j, const Finder &finder)
{
  if (finder.x_.empty() || (finder.x_.size() != finder.y_.size()))
    return;

  if (j.count("fits"))
  {
    nlohmann::json o = j["fits"];
    for (nlohmann::json::iterator it = o.begin(); it != o.end(); ++it)
    {
      SUM4Edge LB(it.value()["background_left"], finder);
      SUM4Edge RB(it.value()["background_right"], finder);

      if (!LB.width() || !RB.width())
        return;

      if (it.value().count("settings"))
        finder_.settings_ = it.value()["settings"];
      else
        finder_.settings_ = finder.settings_;

      //validate background and edges?
      set_data(finder, LB.left(), RB.left());

      LB_ = LB;
      RB_ = RB;

      if (it.value().count("peaks"))
      {
        nlohmann::json p = it.value()["peaks"];
        for (nlohmann::json::iterator it2 = p.begin(); it2 != p.end(); ++it2)
        {
          Peak newpeak(it2.value(), finder_, LB_, RB_);
          peaks_[newpeak.center()] = newpeak;
        }
      }

      background_ = it.value()["background_poly"];
      render();
      save_current_fit(it.value()["description"]);
      peaks_.clear();
    }
  }

  rollback(finder, j["current_fit"]);
}

}
