#include <optimizerBFGS/Region.h>
#include <optimizerBFGS/more_math.h>

#include <core/util/custom_logger.h>

namespace Hypermet
{

Region::Region(CSpectrum& spe, double from_channel, double to_channel)
    : spectrum(spe)
      , first_channel(std::min(from_channel, to_channel))
      , last_channel(std::max(from_channel, to_channel))
{
  try
  {
    width.max(4);
    width.min(0.8);

    short_tail_amplitude.max(1.5);
    short_tail_amplitude.min(0.02);
    short_tail_amplitude.to_fit = true;

    short_tail_slope.max(0.5);
    short_tail_slope.min(0.2);
    short_tail_slope.to_fit = true;

    long_tail_amplitude.max(0.15);
    long_tail_amplitude.min(0.0001);
    long_tail_amplitude.to_fit = true;

    long_tail_slope.max(50);
    long_tail_slope.min(2.5);
    long_tail_slope.to_fit = true;

    right_tail_amplitude.max(0.9);
    right_tail_amplitude.min(0.01);
    right_tail_amplitude.to_fit = true;

    right_tail_slope.max(1.5);
    right_tail_slope.min(0.3);
    right_tail_slope.to_fit = true;

    step_amplitude.max(0.05);
    step_amplitude.min(0.000001);
    step_amplitude.to_fit = true;

    background_slope.to_fit = true;
    background_curve.to_fit = true;

    //  SearchPeaks();

  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("Object Region failed!"));
  }
}

bool Region::left_tail() const
{
  //Get/Set Left Long Tail Status
  return left_tail_enabled_;
}

void Region::left_tail(bool enable)
{
  left_tail_enabled_ = enable;
  long_tail_amplitude.to_fit = enable;
  long_tail_slope.to_fit = enable;
}

bool Region::right_tail() const
{
  return right_tail_enabled_;
}
void Region::right_tail(bool enable)
{
  right_tail_enabled_ = enable;
  right_tail_amplitude.to_fit = enable;
  right_tail_slope.to_fit = enable;
}

bool Region::slope() const
{
  return slope_enabled_;
}

void Region::slope(bool enable)
{
  slope_enabled_ = enable;
  background_slope.to_fit = enable;
}

bool Region::curve() const
{
  return curve_enabled_;
}

void Region::curve(bool enable)
{
  curve_enabled_ = enable;
  background_curve.to_fit = enable;
}

bool Region::step() const
{
  //Get/Set Step Status
  return step_enabled_;
}

void Region::step(bool enable)
{
  step_enabled_ = enable;
  step_amplitude.to_fit = enable;
}

int32_t Region::L(int32_t i, int32_t j, int32_t m)
{
  //M = FWHM
  int32_t ret{0};

  if (j - m <= i && i <= j - 1)
    ret = -1;
  if (j <= i && i <= j + m - 1)
    ret = 2;
  if (j + m <= i && i <= j + 2 * m - 1)
    ret = -1;

  return ret;
}

void Region::find_peaks(uint8_t Threshold)
{
  int32_t m;
  try
  {
    m = static_cast<int32_t>(1.6551 * width.val());
  }
  catch (...)
  {
    m = 3;
  }

  size_t i;
  try
  {
    for (size_t j = first_channel; j <= last_channel; ++j)
    {
      double val = 0;
      double Var = 0;
      for (i = j - m; i <= (j + 2 * m - 1); ++j)
        if (i > 1)
          val = val + L(i, j, m) * spectrum.channels[i];
      Var += square(L(i, j, m)) * spectrum.channels[i];

      //Conv(j - FirstChannel) = val / std::sqrt(Var)
      //if(((Conv(j - FirstChannel - 2) < Conv(j - FirstChannel - 1)) && _
      //(Conv(j - FirstChannel) < Conv(j - FirstChannel - 1)) && _
      //(Conv(j - FirstChannel - 1) > Threshold))) {
      //AddPeak(j - 1, j - 2, j, std::sqrt(spectrum.Channel[j]))
    }
  }
  catch (...)
  {
    ERR("Search Peaks failed!");
  }
}

void Region::add_peak(double Position, double Min, double Max, double Gamma)
{
  try
  {
    peaks.emplace_back();

    peaks.back().GAM.x_index = peaks[std::max(peaks.size() - 2, size_t(0))].GAM.x_index + 2;
    peaks.back().position.x_index = peaks[std::max(peaks.size() - 2, size_t(0))].position.x_index + 2;

    peaks.back().position.min(Min);
    peaks.back().position.max(Max);
    peaks.back().position.val(Position);
    peaks.back().position.uncert_value = 0;

    peaks.back().GAM.val(Gamma);
    peaks.back().GAM.uncert_value = 0;

  }
  catch (...)
  {
    ERR("Add Peak failed!");
  }
}

void Region::remove_peak(size_t index)
{
  try
  {
    if (index >= peaks.size())
    {
      ERR("Can't delete the peak! (invalid index)");
      return;
    }
    if (peaks.size() == 1)
    {
      ERR("Can't delete the only peak!");
      return;
    }
    for (size_t i = index; i < peaks.size() - 1; ++i)
    {
      peaks[i].position = peaks[i + 1].position;
      peaks[i].GAM = peaks[i + 1].GAM;
    } //i
    peaks.resize(peaks.size() - 1);
  }
  catch (...)
  {
    ERR("Delete Peak failed!");
  }
}

double Region::peak_area(size_t PeakIndex) const
{
  return peaks[PeakIndex].GAM.val() * width.val() * (std::sqrt(M_PI) +
      short_tail_amplitude.val() * short_tail_slope.val() * std::exp(-0.25 / square(short_tail_slope.val())) +
      right_tail_amplitude.val() * right_tail_slope.val() * std::exp(-0.25 / square(right_tail_slope.val())));
}

double Region::peak_area_unc(size_t PeakIndex) const
{
  double t = peak_area(PeakIndex);
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

double Region::peak_area_eff(size_t PeakIndex, const Calibration& cal)
{
  double
      eff = cal.efficiency.val(cal.channel_to_energy(peaks[PeakIndex].position.val()));
  //eff = 1 if uninitialized
  return (peaks[PeakIndex].GAM.val() * width.val() * (std::sqrt(M_PI) +
      short_tail_amplitude.val() * short_tail_slope.val() * std::exp(-0.25 / square(short_tail_slope.val())) +
      right_tail_amplitude.val() * right_tail_slope.val() * std::exp(-0.25 / square(right_tail_slope.val())))) / eff;
}

double Region::peak_area_eff_unc(size_t PeakIndex, const Calibration& cal)
{
  double t = peak_area(PeakIndex);
  double
      eff = cal.efficiency.val(cal.channel_to_energy(peaks[PeakIndex].position.val()));
  double sigrel_eff =
      cal.efficiency.sigma_rel(cal.channel_to_energy(peaks[PeakIndex].position.val()));
  //sigrel_eff = 0 if uninitialized
  return (square(std::sqrt(std::sqrt(t) / t)) + square(sigrel_eff)) *
      (t / eff) * std::max(1.0, chi_sq_normalized());
}

size_t Region::fit_var_count() const
{
  size_t n = 2;    //BLN,DEL: always on!
  if (short_tail_amplitude.to_fit)
    n += 1; //AST,BST
  if (short_tail_slope.to_fit)
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
  n += 2 * peaks.size(); //GAM, POS
  return n;
}

void Region::setup_fit()
{
  auto vars = fit_var_count(); // - 1 ?
  fit.resize(vars);
  //fit_gradients.resize(vars);
  //ReDim ChisqGradient(FitVars - 1)
  inv_hessian.resize(vars, vars);

  int32_t shift = 0;

  fit[0] = background_base.x();
  background_base.x_index = 0;
  fit[1] = width.x();
  width.x_index = 1;

  if (short_tail_amplitude.to_fit)
  {
    fit[2] = short_tail_amplitude.x();
    short_tail_amplitude.x_index = 2;
    shift += 1;
  }

  if (short_tail_slope.to_fit)
  {
    fit[2 + shift] = short_tail_slope.x();
    short_tail_slope.x_index = 2 + shift;
    shift += 1;
  }

  if (left_tail_enabled_)
  {
    fit[2 + shift] = long_tail_amplitude.x();
    long_tail_amplitude.x_index = 2 + shift;
    fit[3 + shift] = long_tail_slope.x();
    long_tail_slope.x_index = 3 + shift;
    shift += 2;
  }

  if (right_tail_enabled_)
  {
    fit[2 + shift] = right_tail_amplitude.x();
    right_tail_amplitude.x_index = 2 + shift;
    fit[3 + shift] = right_tail_slope.x();
    right_tail_slope.x_index = 3 + shift;
    shift += 2;
  }

  if (step_enabled_)
  {
    fit[2 + shift] = step_amplitude.x();
    step_amplitude.x_index = 2 + shift;
    shift += 1;
  }

  if (slope_enabled_)
  {
    fit[2 + shift] = background_slope.x();
    background_slope.x_index = 2 + shift;
    shift += 1;
  }

  if (curve_enabled_)
  {
    fit[2 + shift] = background_curve.x();
    background_curve.x_index = 2 + shift;
    shift += 1;
  }

  for (auto& p : peaks)
  {
    fit[2 + shift] = p.GAM.x();
    p.GAM.x_index = 2 + shift;
    fit[3 + shift] = p.position.x();
    p.position.x_index = 3 + shift;
    shift += 2;
  }
}

void Region::store_fit()
{
  double Chisq_norm = std::max(chi_sq_normalized(), 1.0) * 0.5;
  int32_t shift{0};

  double df = degrees_of_freedom();
  for (size_t i = 0; i < fit_var_count(); ++i)
    for (size_t j = 0; j < fit_var_count(); ++j)
      inv_hessian.coeffRef(i, j) *= df;

  background_base.x(fit[0]);
  background_base.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(0, 0) *
      background_base.grad_at(square(fit[0])) * Chisq_norm));
  width.x(fit[1]);
  width.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(1, 1) *
      width.grad_at(square(fit[1])) * Chisq_norm));

  if (short_tail_amplitude.to_fit)
  {
    short_tail_amplitude.x(fit[2]);
    short_tail_amplitude.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2, 2) *
        short_tail_amplitude.grad_at(square(fit[2])) * Chisq_norm));
    shift += 1;
  }

  if (short_tail_slope.to_fit)
  {
    short_tail_slope.x(fit[2 + shift]);
    short_tail_slope.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                              * short_tail_slope.grad_at(square(fit[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (left_tail_enabled_)
  {
    long_tail_amplitude.x(fit[2 + shift]);
    long_tail_amplitude.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                              * long_tail_amplitude.grad_at(square(fit[2 + shift])) * Chisq_norm));
    long_tail_slope.x(fit[3 + shift]);
    long_tail_slope.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(3 + shift, 3 + shift)
                                              * long_tail_slope.grad_at(square(fit[3 + shift])) * Chisq_norm));
    shift += 2;
  }

  if (right_tail_enabled_)
  {
    right_tail_amplitude.x(fit[2 + shift]);
    right_tail_amplitude.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                              * right_tail_amplitude.grad_at(square(fit[2 + shift])) * Chisq_norm));
    right_tail_slope.x(fit[3 + shift]);
    right_tail_slope.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(3 + shift, 3 + shift)
                                              * right_tail_slope.grad_at(square(fit[3 + shift])) * Chisq_norm));
    shift += 2;
  }

  if (step_enabled_)
  {
    step_amplitude.x(fit[2 + shift]);
    step_amplitude.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                              * step_amplitude.grad_at(square(fit[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (slope_enabled_)
  {
    background_slope.x(fit[2 + shift]);
    background_slope.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                              * background_slope.grad_at(square(fit[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (curve_enabled_)
  {
    background_curve.x(fit[2 + shift]);
    background_curve.uncert_value = std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift)
                                              * background_curve.grad_at(square(fit[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  for (auto& p : peaks)
  {
    p.GAM.x(fit[2 + shift]);
    p.GAM.uncert_value =
        std::sqrt(std::abs(inv_hessian.coeff(2 + shift, 2 + shift) *
            p.GAM.grad_at(square(fit[2 + shift])) * Chisq_norm));
    p.position.x(fit[3 + shift]);
    p.position.uncert_value =
        std::sqrt(std::abs(inv_hessian.coeff(3 + shift, 3 + shift) *
            p.position.grad_at(square(fit[3 + shift])) * Chisq_norm));
    shift += 2;
  }
}

void Region::eval_fit(double E, std::vector<double>& ret) const
{
  //returns the value of the fitted curve and the background at Energy E
  //Dim FTotal, FBkg0, FBkg, FPeak As Double
  ret.resize(peaks.size() + 2);
  double _DE, _GAM, _DEL, _BST, _BRT, _BLT;
  try
  {
    //val(1):bkg
    ret[1] = background_base.val();
    if (slope())
      ret[1] += background_slope.val() * (E - first_channel);
    if (curve())
      ret[1] += background_curve.val() * square(E - first_channel);

    _BRT = right_tail_slope.val();
    _BST = short_tail_slope.val();
    _DEL = width.val();
    _BLT = long_tail_slope.val();

    for (size_t i = 0; i < peaks.size(); ++i)
    {
      _DE = E - peaks[i].position.val();
      _GAM = peaks[i].GAM.val();

      if (left_tail())
        ret[1] += _GAM * 0.5 * long_tail_amplitude.val() *
            std::exp(_DE / (_BLT * _DEL)) *
            std::erfc(_DE / _DEL + 0.5 / _BLT);
      if (step())
        ret[1] += step_amplitude.val() * 0.5 * _GAM *
            std::erfc(peaks[i].step_type() * _DE / _DEL);

      ret[i + 2] = _GAM * std::exp(-1.0 * (square(_DE) / square(_DEL)));
      ret[i + 2] += _GAM * 0.5 * short_tail_amplitude.val() *
          std::exp(_DE / (_BST * _DEL)) *
          std::erfc(_DE / _DEL + 0.5 / _BST);
      if (right_tail())
        ret[i + 2] += _GAM * 0.5 * right_tail_amplitude.val() *
            std::exp(-1.0 * _DE / (_BRT * _DEL)) *
            std::erfc(-1.0 * _DE / _DEL + 0.5 / _BRT);
    }
    //val(0):FTotal
    ret[0] = ret[1];
    for (size_t i = 0; i < peaks.size(); ++i)
      ret[0] += ret[i + 2];
  }
  catch (...)
  {
    ERR("Expection");
  }
}

double Region::calc_chi_sq(const std::vector<double>& XVector) const
{
  //Calculates the normalized Chi-square over a region
  try
  {
    chi_squared = 0;

    double _AST = short_tail_amplitude.to_fit ? short_tail_amplitude.val_at(XVector[short_tail_amplitude.x_index]) : short_tail_amplitude.val();
    double _BST = short_tail_slope.to_fit ? short_tail_slope.val_at(XVector[short_tail_slope.x_index]) : short_tail_slope.val();
    double _DEL = width.val_at(XVector[width.x_index]);
    double _ALT = left_tail_enabled_ ? long_tail_amplitude.val_at(XVector[long_tail_amplitude.x_index]) : 0.0;
    double _BLT = left_tail_enabled_ ? long_tail_slope.val_at(XVector[long_tail_slope.x_index]) : 0.0;
    double _SIG = step_enabled_ ? step_amplitude.val_at(XVector[step_amplitude.x_index]) : 0.0;
    double _ART = right_tail_enabled_ ? right_tail_amplitude.val_at(XVector[right_tail_amplitude.x_index]) : 0.0;
    double _BRT = right_tail_enabled_ ? right_tail_slope.val_at(XVector[right_tail_slope.x_index]) : 0.0;

    for (size_t j = first_channel; j <= last_channel; ++j)
    {
      // Background
      double FTotal = background_base.val_at(XVector[background_base.x_index]);
      if (slope_enabled_)
        FTotal += background_slope.val_at(XVector[background_slope.x_index]) * (j - first_channel);
      if (curve_enabled_)
        FTotal += background_curve.val_at(XVector[background_curve.x_index]) * square(j - first_channel);

      for (auto& p : peaks)
      {
        double DE = j - p.position.val_at(XVector[p.position.x_index]);
        double _GAM = p.GAM.val_at(XVector[p.GAM.x_index]);
        if (left_tail_enabled_)
        {
          FTotal += _GAM * 0.5 * _ALT *
              std::exp(DE / (_BLT * _DEL)) *
              std::erfc(DE / _DEL + 0.5 / _BLT);
        }
        if (step_enabled_)
        {
          FTotal += _GAM * 0.5 * _SIG *
              std::erfc(p.step_type() * DE / _DEL);
        }
        FTotal = std::max(FTotal, 0.0);
        //--- Peak components ---
        // Gaussian
        FTotal += _GAM * std::exp(-1.0 * square(DE / _DEL));
        // Short tail
        FTotal += _GAM * 0.5 * _AST *
            std::exp(DE / (_BST * _DEL)) *
            std::erfc(DE / _DEL + 0.5 / _BST);
        if (right_tail_enabled_)
        {
          FTotal += _GAM * 0.5 * _ART *
              std::exp(-1.0 * DE / (_BRT * _DEL)) *
              std::erfc(0.5 / _BRT - DE / _DEL);
        }
      }
      chi_squared += square((spectrum.channels[j] - FTotal) /
          spectrum.weight(j));
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
  return chi_squared / ((last_channel - first_channel) - fit_var_count());
}

size_t Region::degrees_of_freedom() const
{
  return ((last_channel - first_channel) - fit_var_count());
}

void Region::grad_chi_sq(const std::vector<double>& XVector,
                         std::vector<double>& XGradient, double& Chisq) const
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
    std::vector<double> XXGradient(XGradient.size(), 0.0);
    XGradient.assign(XGradient.size(), 0.0);

    double t2;

    Chisq = 0;

    double _AST = short_tail_amplitude.to_fit ? short_tail_amplitude.val_at(XVector[short_tail_amplitude.x_index]) : short_tail_amplitude.val();
    double _BST = short_tail_slope.to_fit ? short_tail_slope.val_at(XVector[short_tail_slope.x_index]) : short_tail_slope.val();
    double _DEL = width.val_at(XVector[width.x_index]);
    double _ALT = left_tail_enabled_ ? long_tail_amplitude.val_at(XVector[long_tail_amplitude.x_index]) : 0.0;
    double _BLT = left_tail_enabled_ ? long_tail_slope.val_at(XVector[long_tail_slope.x_index]) : 0.0;
    double _SIG = step_enabled_ ? step_amplitude.val_at(XVector[step_amplitude.x_index]) : 0.0;
    double _ART = right_tail_enabled_ ? right_tail_amplitude.val_at(XVector[right_tail_amplitude.x_index]) : 0.0;
    double _BRT = right_tail_enabled_ ? right_tail_slope.val_at(XVector[right_tail_slope.x_index]) : 0.0;

    for (size_t j = first_channel; j <= last_channel; ++j)
    {
      //--- Poly Background ---
      double FTotal = background_base.val_at(XVector[background_base.x_index]);
      XXGradient[background_base.x_index] = background_base.grad_at(XVector[background_base.x_index]);
      if (slope_enabled_)
      {
        FTotal += background_slope.val_at(XVector[background_slope.x_index]) *
            (j - first_channel);
        XXGradient[background_slope.x_index] = (j - first_channel);
      }

      if (curve_enabled_)
      {
        FTotal += background_curve.val_at(XVector[background_curve.x_index]) *
            square(j - first_channel);
        XXGradient[background_curve.x_index] = square(j - first_channel);
      }

      for (auto& p : peaks)
      {

        double DE = j - p.position.val_at(XVector[p.position.x_index]);
        double _GAM = p.GAM.val_at(XVector[p.GAM.x_index]);
        double t1 = DE / _DEL;
        //---Left Tail---
        if (left_tail_enabled_)
        {
          double _LongTail = _GAM * 0.5 * _ALT * std::exp(t1 / _BLT) *
              std::erfc(t1 + 0.5 / _BLT);

          FTotal += _LongTail;

          //t2 = (_GAM * _ALT * std::exp(t1 / _BLT) / M_PI ^ (0.5) * std::exp(-(1.0 / (2.0 * _BLT) + t1) ^ 2) * t1 / _DEL)
          t2 = (_GAM * _ALT * std::exp(t1 / _BLT) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * _BLT) + t1)) / _DEL);
          XXGradient[width.x_index] += width.grad_at(XVector[width.x_index])
              * (-1.0 * t1 / (_DEL * _BLT) * _LongTail + t2 * t1);
          XXGradient[p.position.x_index] += -1.0 / (_BLT * _DEL) *
              _LongTail + t2;
          XXGradient[p.GAM.x_index] += _LongTail / _GAM;

          XXGradient[long_tail_amplitude.x_index] += _LongTail / _ALT *
              long_tail_amplitude.grad_at(XVector[long_tail_amplitude.x_index]);
          XXGradient[long_tail_slope.x_index] += long_tail_slope.grad_at(XVector[long_tail_slope.x_index])
              * ((-1.0 * t1 / square(_BLT)) *
                  _LongTail + (_DEL / (2.0 * square(_BLT)) * t2));

        }
        //---Step---
        if (step_enabled_)
        {
          double _StepBkg = _SIG * 0.5 * _GAM *
              std::erfc(p.step_type() * t1);
          FTotal += _StepBkg;

          XXGradient[width.x_index] += width.grad_at(XVector[width.x_index]) *
              (_GAM * _SIG * p.step_type() / std::sqrt(M_PI) *
                  std::exp(-DE / _DEL * t1) * t1 / _DEL);
          XXGradient[p.GAM.x_index] += _StepBkg / _GAM;
          XXGradient[step_amplitude.x_index] += _StepBkg / _SIG *
              step_amplitude.grad_at(XVector[step_amplitude.x_index]);
        }
        FTotal = std::max(FTotal, 0.0);

        //---Gaussian---
        double _Gauss = _GAM * std::exp(-1.0 * square(t1));
        FTotal += _Gauss;

        XXGradient[width.x_index] += width.grad_at(XVector[width.x_index]) *
            (2.0 * square(t1) / _DEL * _Gauss);

        XXGradient[p.position.x_index] += 2.0 * t1 / _DEL * _Gauss;
        XXGradient[p.GAM.x_index] += _Gauss / _GAM;

        //---Short Tail---

        double _ShortTail = _GAM * 0.5 * _AST * std::exp(t1 / _BST) *
            std::erfc(t1 + 0.5 / _BST);
        FTotal += _ShortTail;

        //t2 = (_GAM * _AST * std::exp(t1 / _BST) / M_PI ^ (0.5) * std::exp(-1.0 * (1.0 / (2.0 * _BST) + t1) ^ 2) * t1 / _DEL)
        t2 = (_GAM * _AST * std::exp(t1 / _BST) / std::sqrt(M_PI) *
            std::exp(-1.0 * square(1.0 / (2.0 * _BST) + t1)) / _DEL);
        XXGradient[width.x_index] += width.grad_at(XVector[width.x_index]) *
            (-1.0 * t1 / (_DEL * _BST) * _ShortTail + t2 * t1);

        XXGradient[p.position.x_index] += -1.0 / (_BST * _DEL) *
            _ShortTail + t2;
        XXGradient[p.GAM.x_index] += _ShortTail / _GAM;

        if (short_tail_amplitude.to_fit)
          XXGradient[short_tail_amplitude.x_index] += _ShortTail / _AST *
              short_tail_amplitude.grad_at(XVector[short_tail_amplitude.x_index]);
        if (short_tail_slope.to_fit)
          XXGradient[short_tail_slope.x_index] += short_tail_slope.grad_at(XVector[short_tail_slope.x_index]) *
              ((-1.0 * t1 / square(_BST)) *
                  _ShortTail + (_DEL / (2.0 * square(_BST)) * t2));

        //---Right Tail---
        if (right_tail_enabled_)
        {
          double _RightTail = _GAM * 0.5 * _ART *
              std::exp(-1.0 * t1 / _BRT) *
              std::erfc(0.5 / _BRT - t1);
          FTotal += _RightTail;

          //t2 = (_GAM * _ART * std::exp(-1.0 * t1 / _BRT) / M_PI ^ (0.5) * std::exp(-(1.0 / (2.0 * _BRT) - t1) ^ 2) * t1 / _DEL)
          t2 = (_GAM * _ART * std::exp(-1.0 * t1 / _BRT) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * _BRT) - t1)) / _DEL);
          XXGradient[width.x_index] += width.grad_at(XVector[width.x_index]) *
              ((t1 / (_DEL * _BRT) * _RightTail - t2 * t1));

          XXGradient[p.position.x_index] += 1.0 / (_BRT * _DEL) *
              _RightTail - t2;
          XXGradient[p.GAM.x_index] +=
              _RightTail / _GAM;

          XXGradient[right_tail_amplitude.x_index] += _RightTail / _ART *
              right_tail_amplitude.grad_at(XVector[right_tail_amplitude.x_index]);
          XXGradient[right_tail_slope.x_index] += right_tail_slope.grad_at(XVector[right_tail_slope.x_index])
              * ((t1 / square(_BRT)) * _RightTail + (_DEL / (2.0 *
                  square(_BRT)) * t2));
        }

        //XXGradient(DEL.x_index) *= DEL.GradAt(XVector(DEL.x_index))
        XXGradient[p.GAM.x_index] *= p.GAM.grad_at(XVector[p.GAM.x_index]);
        XXGradient[p.position.x_index] *= p.position.grad_at(XVector[p.position.x_index]);
      } //Peak

      double t3 = -2 * (spectrum.channels[j] - FTotal) /
          square(spectrum.weight(j));

      for (size_t k = 0; k < fit_var_count(); ++k)
      {
        XGradient[k] += XXGradient[k] * t3;
        XXGradient[k] = 0.0;
      }
      Chisq += square((spectrum.channels[j] - FTotal) / spectrum.weight(j));
    } //j //Channel
    //Chisq /= df
  }
  catch (...)
  {
    ERR("Error in GradChiSq");
  }
}

}
