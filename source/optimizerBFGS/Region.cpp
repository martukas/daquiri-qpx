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

  width_.max(4);
  width_.min(0.8);

  short_tail_amplitude_.max(1.5);
  short_tail_amplitude_.min(0.02);
  short_tail_amplitude_.to_fit = true;

  short_tail_slope_.max(0.5);
  short_tail_slope_.min(0.2);
  short_tail_slope_.to_fit = true;

  long_tail_amplitude_.max(0.15);
  long_tail_amplitude_.min(0.0001);
  long_tail_amplitude_.to_fit = true;

  long_tail_slope_.max(50);
  long_tail_slope_.min(2.5);
  long_tail_slope_.to_fit = true;

  right_tail_amplitude_.max(0.9);
  right_tail_amplitude_.min(0.01);
  right_tail_amplitude_.to_fit = true;

  right_tail_slope_.max(1.5);
  right_tail_slope_.min(0.3);
  right_tail_slope_.to_fit = true;

  step_amplitude_.max(0.05);
  step_amplitude_.min(0.000001);
  step_amplitude_.to_fit = true;

  background_slope_.to_fit = true;
  background_curve_.to_fit = true;

  //  SearchPeaks(); avoid this for the sake of batch fitting
}

bool Region::left_tail() const
{
  return left_tail_enabled_;
}

void Region::left_tail(bool enable)
{
  left_tail_enabled_ = enable;
  long_tail_amplitude_.to_fit = enable;
  long_tail_slope_.to_fit = enable;
}

bool Region::right_tail() const
{
  return right_tail_enabled_;
}
void Region::right_tail(bool enable)
{
  right_tail_enabled_ = enable;
  right_tail_amplitude_.to_fit = enable;
  right_tail_slope_.to_fit = enable;
}

bool Region::slope() const
{
  return slope_enabled_;
}

void Region::slope(bool enable)
{
  slope_enabled_ = enable;
  background_slope_.to_fit = enable;
}

bool Region::curve() const
{
  return curve_enabled_;
}

void Region::curve(bool enable)
{
  curve_enabled_ = enable;
  background_curve_.to_fit = enable;
}

bool Region::step() const
{
  return step_enabled_;
}

void Region::step(bool enable)
{
  step_enabled_ = enable;
  step_amplitude_.to_fit = enable;
}

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
  auto m = static_cast<int32_t>(1.6551 * width_.val());
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
  return peaks_[index].amplitude.val() * width_.val() * (std::sqrt(M_PI) +
      short_tail_amplitude_.val() * short_tail_slope_.val() *
          std::exp(-0.25 / square(short_tail_slope_.val())) +
      right_tail_amplitude_.val() * right_tail_slope_.val() *
          std::exp(-0.25 / square(right_tail_slope_.val())));
}

double Region::peak_area_unc(size_t index) const
{
  // \todo make this more rigorous
  double t = peak_area(index);
  //, i, j As Integer
  //Dim cs As Double = ChisqNorm() * 0.5
  //for( i = 0 To FitVars - 1
  //for( j = 0 To i - 1
  //t += fit_gradients(i) * fit_gradients(j) * Hessinv.coeff(i, j) * cs
  //} //j
  //} //i
  //
  //Dim Bracket As Double = (std::sqrt(M_PI) + AST.val() * BST.val() * std::exp(-0.25 / BST.val() ^ 2) + ART.val() * BRT.val() * std::exp(-0.25 / BRT.val() ^ 2))
  //(dGAM/dX*dArea/dGAM)^2*Var(X)*Chisq
  //t = (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) ^ 2 * Hessinv.coeff(DEL.x_index, DEL.x_index) * cs
  //(dGAM/dX*dArea/dGAM)*(dDEL/dY*dArea/dDEL)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) * (DEL.GradAt(DEL.X) * Peak(PeakIndex).GAM.val() * Bracket) * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, DEL.x_index) * cs
  //if(AST.to_fit = true) {
  ////(dGAM/dX*dArea/dGAM)*(dAST/dY*dArea/dAST)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) * (AST.GradAt(AST.X) * Peak(PeakIndex).GAM.val() * DEL.val() * BST.val() * std::exp(-0.25 / BST.val() ^ 2)) * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, AST.x_index) * cs
  //}
  //if(BST.to_fit = true) {
  //(dGAM/dX*dArea/dGAM)*(dBST/dY*dArea/dBST)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) * (BST.GradAt(BST.X) * Peak(PeakIndex).GAM.val() * DEL.val() * AST.val() * (1 + 0.5 / BST.val() ^ 2) * std::exp(-0.25 / BST.val() ^ 2)) * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, AST.x_index) * cs
  //}
  return std::sqrt(t) * std::max(1.0, chi_sq_normalized());
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
  size_t n = 2;    //BLN,DEL: always on!
  if (short_tail_amplitude_.to_fit)
    n += 1; //AST,BST
  if (short_tail_slope_.to_fit)
    n += 1;
  if (left_tail_enabled_)
    n += 2; //ALT,BLT
  if (right_tail_enabled_)
    n += 2; //ART,BRT
  if (step())
    n += 1;
  if (slope_enabled_)
    n += 1;
  if (curve_enabled_)
    n += 1;
  n += 2 * peaks_.size(); //GAM, POS
  return n;
}

void Region::setup_fit()
{
  auto vars = fit_var_count(); // - 1 ?
  current_fit.resize(vars);
  //fit_gradients.resize(vars);
  //ReDim ChisqGradient(FitVars - 1)

  int32_t shift = 0;

  current_fit[0] = background_base_.x();
  background_base_.x_index = 0;
  current_fit[1] = width_.x();
  width_.x_index = 1;

  if (short_tail_amplitude_.to_fit)
  {
    current_fit[2] = short_tail_amplitude_.x();
    short_tail_amplitude_.x_index = 2;
    shift += 1;
  }

  if (short_tail_slope_.to_fit)
  {
    current_fit[2 + shift] = short_tail_slope_.x();
    short_tail_slope_.x_index = 2 + shift;
    shift += 1;
  }

  if (left_tail_enabled_)
  {
    current_fit[2 + shift] = long_tail_amplitude_.x();
    long_tail_amplitude_.x_index = 2 + shift;
    current_fit[3 + shift] = long_tail_slope_.x();
    long_tail_slope_.x_index = 3 + shift;
    shift += 2;
  }

  if (right_tail_enabled_)
  {
    current_fit[2 + shift] = right_tail_amplitude_.x();
    right_tail_amplitude_.x_index = 2 + shift;
    current_fit[3 + shift] = right_tail_slope_.x();
    right_tail_slope_.x_index = 3 + shift;
    shift += 2;
  }

  if (step_enabled_)
  {
    current_fit[2 + shift] = step_amplitude_.x();
    step_amplitude_.x_index = 2 + shift;
    shift += 1;
  }

  if (slope_enabled_)
  {
    current_fit[2 + shift] = background_slope_.x();
    background_slope_.x_index = 2 + shift;
    shift += 1;
  }

  if (curve_enabled_)
  {
    current_fit[2 + shift] = background_curve_.x();
    background_curve_.x_index = 2 + shift;
    shift += 1;
  }

  for (auto& p : peaks_)
  {
    current_fit[2 + shift] = p.amplitude.x();
    p.amplitude.x_index = 2 + shift;
    current_fit[3 + shift] = p.position.x();
    p.position.x_index = 3 + shift;
    shift += 2;
  }
}

void Region::store_fit()
{
  double chisq_norm = std::max(chi_sq_normalized(), 1.0) * 0.5;
  int32_t shift{0};

  double df = degrees_of_freedom();
  for (size_t i = 0; i < fit_var_count(); ++i)
    for (size_t j = 0; j < fit_var_count(); ++j)
      inv_hessian.coeffRef(i, j) *= df;

  background_base_.x(current_fit[0]);
  background_base_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(0, 0) *
      background_base_.grad_at(square(current_fit[0])) * chisq_norm));
  width_.x(current_fit[1]);
  width_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(1, 1) *
      width_.grad_at(square(current_fit[1])) * chisq_norm));

  if (short_tail_amplitude_.to_fit)
  {
    short_tail_amplitude_.x(current_fit[2]);
    short_tail_amplitude_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2, 2) *
        short_tail_amplitude_.grad_at(square(current_fit[2])) * chisq_norm));
    shift += 1;
  }

  if (short_tail_slope_.to_fit)
  {
    short_tail_slope_.x(current_fit[2 + shift]);
    short_tail_slope_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                                            * short_tail_slope_.grad_at(square(current_fit[2 + shift]))
                                                            * chisq_norm));
    shift += 1;
  }

  if (left_tail_enabled_)
  {
    long_tail_amplitude_.x(current_fit[2 + shift]);
    long_tail_amplitude_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                                               * long_tail_amplitude_.grad_at(square(current_fit[2
                                                                   + shift]))
                                                               * chisq_norm));
    long_tail_slope_.x(current_fit[3 + shift]);
    long_tail_slope_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(3 + shift, 3 + shift)
                                                           * long_tail_slope_.grad_at(square(current_fit[3 + shift]))
                                                           * chisq_norm));
    shift += 2;
  }

  if (right_tail_enabled_)
  {
    right_tail_amplitude_.x(current_fit[2 + shift]);
    right_tail_amplitude_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                                                * right_tail_amplitude_.grad_at(square(current_fit[2
                                                                    + shift]))
                                                                * chisq_norm));
    right_tail_slope_.x(current_fit[3 + shift]);
    right_tail_slope_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(3 + shift, 3 + shift)
                                                            * right_tail_slope_.grad_at(square(current_fit[3 + shift]))
                                                            * chisq_norm));
    shift += 2;
  }

  if (step_enabled_)
  {
    step_amplitude_.x(current_fit[2 + shift]);
    step_amplitude_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                                          * step_amplitude_.grad_at(square(current_fit[2 + shift]))
                                                          * chisq_norm));
    shift += 1;
  }

  if (slope_enabled_)
  {
    background_slope_.x(current_fit[2 + shift]);
    background_slope_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                                            * background_slope_.grad_at(square(current_fit[2 + shift]))
                                                            * chisq_norm));
    shift += 1;
  }

  if (curve_enabled_)
  {
    background_curve_.x(current_fit[2 + shift]);
    background_curve_.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                                            * background_curve_.grad_at(square(current_fit[2 + shift]))
                                                            * chisq_norm));
    shift += 1;
  }

  for (auto& p : peaks_)
  {
    p.amplitude.x(current_fit[2 + shift]);
    p.amplitude.uncert_value =
        std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift) *
            p.amplitude.grad_at(square(current_fit[2 + shift])) * chisq_norm));
    p.position.x(current_fit[3 + shift]);
    p.position.uncert_value =
        std::sqrt(std::abs(inv_hessian.coeff(3 + shift, 3 + shift) *
            p.position.grad_at(square(current_fit[3 + shift])) * chisq_norm));
    shift += 2;
  }
}

void Region::eval_fit(double pos, std::vector<double>& ret) const
{
  //returns the value of the fitted curve and the background at position pos
  //Dim FTotal, FBkg0, FBkg, FPeak As Double
  ret.resize(peaks_.size() + 2);
  try
  {
    //val(1):bkg
    ret[1] = background_base_.val();
    if (slope())
      ret[1] += background_slope_.val() * (pos - first_channel);
    if (curve())
      ret[1] += background_curve_.val() * square(pos - first_channel);

    double right_slope = right_tail_slope_.val();
    double short_slope = short_tail_slope_.val();
    double width = width_.val();
    double left_slope = long_tail_slope_.val();

    for (size_t i = 0; i < peaks_.size(); ++i)
    {
      double x_rel = pos - peaks_[i].position.val();
      double ampl = peaks_[i].amplitude.val();
      double half_ampl = 0.5 * ampl;
      double spread = x_rel / width;

      // background
      if (left_tail())
        ret[1] += half_ampl * long_tail_amplitude_.val() *
            std::exp(spread / left_slope) *
            std::erfc(0.5 / left_slope + spread);
      if (step())
        ret[1] += step_amplitude_.val() * 0.5 * ampl *
            std::erfc(peaks_[i].step_type() * spread);

      // peak
      ret[i + 2] = ampl * std::exp(-1.0 * square(spread));
      ret[i + 2] += half_ampl * short_tail_amplitude_.val() *
          std::exp(spread / short_slope) *
          std::erfc(0.5 / short_slope + spread);
      if (right_tail())
        ret[i + 2] += half_ampl * right_tail_amplitude_.val() *
            std::exp(-1.0 * spread / right_slope) *
            std::erfc(0.5 / right_slope - spread);
    }
    //val(0):FTotal
    ret[0] = ret[1];
    for (size_t i = 0; i < peaks_.size(); ++i)
      ret[0] += ret[i + 2];
  }
  catch (...)
  {
    ERR("Failed to eval fit");
  }
}

double Region::calc_chi_sq(const std::vector<double>& fit) const
{
  //Calculates the normalized Chi-square over a region
  try
  {
    chi_squared = 0;

    double short_ampl = short_tail_amplitude_.to_fit ?
                        short_tail_amplitude_.val_at(fit[short_tail_amplitude_.x_index]) :
                        short_tail_amplitude_.val();
    double short_slope = short_tail_slope_.to_fit ?
                         short_tail_slope_.val_at(fit[short_tail_slope_.x_index]) :
                         short_tail_slope_.val();
    double width = width_.val_at(fit[width_.x_index]);
    double left_ampl = left_tail_enabled_ ?
                       long_tail_amplitude_.val_at(fit[long_tail_amplitude_.x_index]) :
                       0.0;
    double left_slope = left_tail_enabled_ ?
                        long_tail_slope_.val_at(fit[long_tail_slope_.x_index]) :
                        0.0;
    double step_ampl = step_enabled_ ?
                       step_amplitude_.val_at(fit[step_amplitude_.x_index]) :
                       0.0;
    double right_ampl = right_tail_enabled_ ?
                        right_tail_amplitude_.val_at(fit[right_tail_amplitude_.x_index]) :
                        0.0;
    double right_slope = right_tail_enabled_ ?
                         right_tail_slope_.val_at(fit[right_tail_slope_.x_index]) :
                         0.0;

    for (size_t pos = first_channel; pos <= last_channel; ++pos)
    {
      // Background
      double FTotal = background_base_.val_at(fit[background_base_.x_index]);
      if (slope_enabled_)
        FTotal += background_slope_.val_at(fit[background_slope_.x_index])
            * (pos - first_channel);
      if (curve_enabled_)
        FTotal += background_curve_.val_at(fit[background_curve_.x_index])
            * square(pos - first_channel);

      for (auto& p : peaks_)
      {
        double x_rel = pos - p.position.val_at(fit[p.position.x_index]);
        double ampl = p.amplitude.val_at(fit[p.amplitude.x_index]);
        double half_ampl = 0.5 * ampl;
        double spread = x_rel / width;

        if (left_tail_enabled_)
        {
          FTotal += half_ampl * left_ampl *
              std::exp(spread / left_slope) *
              std::erfc(0.5 / left_slope + spread);
        }
        if (step_enabled_)
        {
          FTotal += half_ampl * step_ampl *
              std::erfc(p.step_type() * spread);
        }
        FTotal = std::max(FTotal, 0.0);
        //--- Peak components ---
        // Gaussian
        FTotal += ampl * std::exp(-1.0 * square(spread));
        // Short tail
        FTotal += half_ampl * short_ampl *
            std::exp(spread / short_slope) *
            std::erfc(0.5 / short_slope + spread);
        if (right_tail_enabled_)
        {
          FTotal += half_ampl * right_ampl *
              std::exp(-1.0 * spread / right_slope) *
              std::erfc(0.5 / right_slope - spread);
        }
      }
      chi_squared += square((spectrum.channels[pos] - FTotal) /
          spectrum.weight(pos));
    } //Channel

    return chi_squared;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("CalcChiSq failed"));
  }
}

double Region::chi_sq_normalized() const
{
  return chi_squared / static_cast<double>(degrees_of_freedom());
}

size_t Region::degrees_of_freedom() const
{
  return ((last_channel - first_channel) - fit_var_count());
}

void Region::grad_chi_sq(const std::vector<double>& fit,
                         std::vector<double>& gradients, double& Chisq) const
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
  try
  {
    // zero-out arrays
    gradients.assign(gradients.size(), 0.0);
    auto chan_gradients = gradients;

    double t2;

    Chisq = 0;

    double short_ampl = short_tail_amplitude_.to_fit ?
                        short_tail_amplitude_.val_at(fit[short_tail_amplitude_.x_index]) :
                        short_tail_amplitude_.val();
    double short_slope = short_tail_slope_.to_fit ?
                         short_tail_slope_.val_at(fit[short_tail_slope_.x_index]) :
                         short_tail_slope_.val();
    double width = width_.val_at(fit[width_.x_index]);
    double left_ampl = left_tail_enabled_ ?
                       long_tail_amplitude_.val_at(fit[long_tail_amplitude_.x_index]) :
                       0.0;
    double left_slope = left_tail_enabled_ ?
                        long_tail_slope_.val_at(fit[long_tail_slope_.x_index]) :
                        0.0;
    double step_ampl = step_enabled_ ?
                       step_amplitude_.val_at(fit[step_amplitude_.x_index]) :
                       0.0;
    double right_ampl = right_tail_enabled_ ?
                        right_tail_amplitude_.val_at(fit[right_tail_amplitude_.x_index]) :
                        0.0;
    double right_slope = right_tail_enabled_ ?
                         right_tail_slope_.val_at(fit[right_tail_slope_.x_index]) :
                         0.0;

    for (size_t pos = first_channel; pos <= last_channel; ++pos)
    {
      chan_gradients.assign(chan_gradients.size(), 0.0);

      //--- Poly Background ---
      double FTotal = background_base_.val_at(fit[background_base_.x_index]);
      chan_gradients[background_base_.x_index] = background_base_.grad_at(fit[background_base_.x_index]);
      if (slope_enabled_)
      {
        FTotal += background_slope_.val_at(fit[background_slope_.x_index])
            * (pos - first_channel);
        chan_gradients[background_slope_.x_index] = (pos - first_channel);
      }

      if (curve_enabled_)
      {
        FTotal += background_curve_.val_at(fit[background_curve_.x_index])
            * square(pos - first_channel);
        chan_gradients[background_curve_.x_index] = square(pos - first_channel);
      }

      for (auto& p : peaks_)
      {

        double x_rel = pos - p.position.val_at(fit[p.position.x_index]);
        double ampl = p.amplitude.val_at(fit[p.amplitude.x_index]);
        double half_ampl = 0.5 * ampl;
        double spread = x_rel / width;
        //---Left Tail---
        if (left_tail_enabled_)
        {
          double long_tail = half_ampl * left_ampl *
              std::exp(spread / left_slope) *
              std::erfc(spread + 0.5 / left_slope);

          FTotal += long_tail;

          //t2 = (ampl * left_ampl * std::exp(spread / left_slope) / ::sqrt(M_PI) * std::exp(-(1.0 / (2.0 * left_slope) + spread) ^ 2) * spread / width)
          t2 = (ampl * left_ampl * std::exp(spread / left_slope) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * left_slope) + spread)) / width);
          chan_gradients[width_.x_index] += width_.grad_at(fit[width_.x_index])
              * (-1.0 * spread / (width * left_slope) * long_tail + t2 * spread);
          chan_gradients[p.position.x_index] += -1.0 / (left_slope * width) *
              long_tail + t2;
          chan_gradients[p.amplitude.x_index] += long_tail / ampl;

          chan_gradients[long_tail_amplitude_.x_index] += long_tail / left_ampl *
              long_tail_amplitude_.grad_at(fit[long_tail_amplitude_.x_index]);
          chan_gradients[long_tail_slope_.x_index] += long_tail_slope_.grad_at(fit[long_tail_slope_.x_index])
              * ((-1.0 * spread / square(left_slope)) *
                  long_tail + (width / (2.0 * square(left_slope)) * t2));

        }
        //---Step---
        if (step_enabled_)
        {
          double step =  half_ampl * step_ampl * std::erfc(p.step_type() * spread);
          FTotal += step;

          chan_gradients[width_.x_index] += width_.grad_at(fit[width_.x_index]) *
              (ampl * step_ampl * p.step_type() / std::sqrt(M_PI) *
                  std::exp(-x_rel / width * spread) * spread / width);
          chan_gradients[p.amplitude.x_index] += step / ampl;
          chan_gradients[step_amplitude_.x_index] += step / step_ampl *
              step_amplitude_.grad_at(fit[step_amplitude_.x_index]);
        }
        FTotal = std::max(FTotal, 0.0);

        //---Gaussian---
        double gauss = ampl * std::exp(-1.0 * square(spread));
        FTotal += gauss;
        chan_gradients[width_.x_index] += width_.grad_at(fit[width_.x_index]) *
            (2.0 * square(spread) / width * gauss);
        chan_gradients[p.position.x_index] += 2.0 * spread / width * gauss;
        chan_gradients[p.amplitude.x_index] += gauss / ampl;

        //---Short Tail---
        double short_tail = half_ampl * short_ampl * std::exp(spread / short_slope) *
            std::erfc(spread + 0.5 / short_slope);
        FTotal += short_tail;
        //t2 = (ampl * short_ampl * std::exp(spread / short_slope) / std::sqrt(M_PI) *
        //    std::exp(-1.0 * square(1.0 / (2.0 * short_slope) + spread)) * spread / width)
        t2 = (ampl * short_ampl * std::exp(spread / short_slope) / std::sqrt(M_PI) *
            std::exp(-1.0 * square(1.0 / (2.0 * short_slope) + spread)) / width);
        chan_gradients[width_.x_index] += width_.grad_at(fit[width_.x_index]) *
            (-1.0 * spread / (width * short_slope) * short_tail + t2 * spread);
        chan_gradients[p.position.x_index] += -1.0 / (short_slope * width) *
            short_tail + t2;
        chan_gradients[p.amplitude.x_index] += short_tail / ampl;

        if (short_tail_amplitude_.to_fit)
          chan_gradients[short_tail_amplitude_.x_index] += short_tail / short_ampl *
              short_tail_amplitude_.grad_at(fit[short_tail_amplitude_.x_index]);
        if (short_tail_slope_.to_fit)
          chan_gradients[short_tail_slope_.x_index] += short_tail_slope_.grad_at(fit[short_tail_slope_.x_index]) *
              ((-1.0 * spread / square(short_slope)) *
                  short_tail + (width / (2.0 * square(short_slope)) * t2));

        //---Right Tail---
        if (right_tail_enabled_)
        {
          double right_tail = half_ampl * right_ampl *
              std::exp(-1.0 * spread / right_slope) *
              std::erfc(0.5 / right_slope - spread);
          FTotal += right_tail;

          //t2 = (ampl * right_ampl * std::exp(-1.0 * spread / right_slope) / std::sqrt(M_PI) *
          //  std::exp(-square(1.0 / (2.0 * right_slope) - spread)) * spread / width)
          t2 = (ampl * right_ampl * std::exp(-1.0 * spread / right_slope) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * right_slope) - spread)) / width);
          chan_gradients[width_.x_index] += width_.grad_at(fit[width_.x_index]) *
              ((spread / (width * right_slope) * right_tail - t2 * spread));

          chan_gradients[p.position.x_index] += 1.0 / (right_slope * width) *
              right_tail - t2;
          chan_gradients[p.amplitude.x_index] +=
              right_tail / ampl;

          chan_gradients[right_tail_amplitude_.x_index] += right_tail / right_ampl *
              right_tail_amplitude_.grad_at(fit[right_tail_amplitude_.x_index]);
          chan_gradients[right_tail_slope_.x_index] += right_tail_slope_.grad_at(fit[right_tail_slope_.x_index])
              * ((spread / square(right_slope)) * right_tail + (width / (2.0 *
                  square(right_slope)) * t2));
        }

        //XXGradient(DEL.x_index) *= DEL.GradAt(XVector(DEL.x_index))
        chan_gradients[p.amplitude.x_index] *= p.amplitude.grad_at(fit[p.amplitude.x_index]);
        chan_gradients[p.position.x_index] *= p.position.grad_at(fit[p.position.x_index]);
      } //Peak

      double t3 = -2.0 * (spectrum.channels[pos] - FTotal) /
          square(spectrum.weight(pos));

      for (size_t var = 0; var < fit_var_count(); ++var)
        gradients[var] += chan_gradients[var] * t3;
      Chisq += square((spectrum.channels[pos] - FTotal) / spectrum.weight(pos));
    }
    //Chisq /= df
  }
  catch (...)
  {
    ERR("Error in GradChiSq");
  }
}

}
