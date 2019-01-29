#include <core/fitting/hypermet/Hypermet.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

double Hypermet::Components::peak_skews() const
{
  return gaussian + short_tail + right_tail;
}

double Hypermet::Components::step_tail() const
{
  return long_tail + step;
}

double Hypermet::Components::all() const
{
  return peak_skews() + step_tail();
}

Hypermet::Hypermet()
{
  width_.bound(0.8, 4.0);

  short_tail.amplitude.bound(0.02, 1.5);
  short_tail.slope.bound(0.2, 0.5);

  right_tail.amplitude.bound(0.01, 0.9);
  right_tail.slope.bound(0.3, 1.5);

  long_tail.amplitude.bound(0.0001, 0.15);
  long_tail.slope.bound(2.5, 50);

  step.amplitude.bound(0.000001, 0.05);
}

void Hypermet::apply_defaults(const Hypermet& other)
{
  if (!width_override)
    width_ = other.width_;
  if (!short_tail.override)
    short_tail = other.short_tail;
  if (!right_tail.override)
    right_tail = other.right_tail;
  if (!long_tail.override)
    long_tail = other.long_tail;
}

void Hypermet::force_defaults(const Hypermet& other)
{
  width_ = other.width_;
  short_tail = other.short_tail;
  right_tail = other.right_tail;
  long_tail = other.long_tail;
}

Hypermet Hypermet::gaussian_only() const
{
  Hypermet ret = *this;
  ret.short_tail.override = true;
  ret.short_tail.enabled = false;
  ret.right_tail.override = true;
  ret.right_tail.enabled = false;
  ret.long_tail.override = true;
  ret.long_tail.enabled = false;
  ret.step.override = true;
  ret.step.enabled = false;
  return ret;
}

bool Hypermet::is_gaussian_only() const
{
  return (!short_tail.enabled && !right_tail.enabled && !long_tail.enabled && !step.enabled);
}

bool Hypermet::sanity_check(double min_x, double max_x) const
{
  auto amp = amplitude.val();
  auto wid = width_.val();
  auto pos = position.val();
  return std::isfinite(amp) && (amp > 0.0) &&
      std::isfinite(wid) && (wid > 0.0) &&
      std::isfinite(pos) && (min_x < pos) && (pos < max_x);
}

bool Hypermet::full_energy_peak() const
{
  return (step.flip(1.0) > 0);
}

void Hypermet::full_energy_peak(bool flag)
{
  if (flag)
    step.side = Side::left;
  else
    step.side = Side::right;
}

bool Hypermet::operator<(const Hypermet& other) const
{
  return position.val() < other.position.val();
}

UncertainDouble Hypermet::peak_position() const
{
  return {position.val(), position.uncert_value};
}

UncertainDouble Hypermet::peak_energy(const HCalibration& cal) const
{
  // \todo if there is a curve, does the slope not differ with x?
  return {cal.channel_to_energy(position.val()),
          cal.energy_slope() * position.uncert_value};
}

UncertainDouble Hypermet::peak_energy(const Calibration& cal) const
{
  return {cal.transform(position.val()),
          cal.function()->derivative(position.val()) * position.uncert_value};
}

UncertainDouble Hypermet::area() const
{
  auto a = amplitude.val() * width_.val() * (std::sqrt(M_PI) +
      short_tail.amplitude.val() * short_tail.slope.val() *
          std::exp(-0.25 / square(short_tail.slope.val())) +
      right_tail.amplitude.val() * right_tail.slope.val() *
          std::exp(-0.25 / square(right_tail.slope.val())));

  // \todo make this more rigorous
//  double cs = chi_sq_norm * 0.5;
//  int i, j;
//  for( i = 0 To FitVars - 1)
//  for( j = 0 To i - 1)
//  t += fit_gradients(i) * fit_gradients(j) * Hessinv.coeff(i, j) * cs;

//  //(dGAM/dX*dArea/dGAM)^2*Var(X)*Chisq
//  t = square(amplitude.grad() * width_.val() * a) * cs;
//  // * Hessinv.coeff(DEL.x_index, DEL.x_index);
//  //(dGAM/dX*dArea/dGAM)*(dDEL/dY*dArea/dDEL)*Covar(X,Y)*Chisq
//  t += square(amplitude.grad() * width_.val() * a) * cs;
//  // * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, DEL.x_index)
//  if (short_tail.amplitude.to_fit)
//  {
//    //(dGAM/dX*dArea/dGAM)*(dAST/dY*dArea/dAST)*Covar(X,Y)*Chisq
//    t += (amplitude.grad() * width_.val() * a) *
//        (short_tail.amplitude.grad() * amplitude.val() * width_.val() *
//            short_tail.slope.val() * std::exp(-0.25 / square(short_tail.slope.val()))) * cs;
//    // * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, AST.x_index)
//  }
//  if (short_tail.slope.to_fit)
//  {
//    // (dGAM/dX*dArea/dGAM)*(dBST/dY*dArea/dBST)*Covar(X,Y)*Chisq
//    t += (amplitude.grad() * width_.val() * a) *
//        (short_tail.slope.grad() * amplitude.val() * width_.val()
//            * short_tail.amplitude.val() * (1 + 0.5 / square(short_tail.slope.val()))
//            * std::exp(-0.25 / square(short_tail.slope.val()))) * cs;
//    // * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, AST.x_index)
//  }
  return {a, std::sqrt(a) * std::max(1.0, chi_sq_norm)};
}

UncertainDouble Hypermet::peak_area_eff(const HCalibration& cal) const
{
  // \todo should this also be more rigorous, like area() ?
  auto e = peak_energy(cal).value();
  auto a = area().value();

  double eff{0.0};
  double sigrel_eff{0.0};
  if (cal.efficiency.initialized())
  {
    eff = cal.efficiency.val(e);
    sigrel_eff = cal.efficiency.sigma_rel(e);
  }
  return {a / eff,
          square(std::sqrt(std::sqrt(a) / a)) + square(sigrel_eff) *
              (a / eff) * std::max(1.0, chi_sq_norm)};
}

UncertainDouble Hypermet::fwhm() const
{
  // \todo multiply with sqrt(log(2.0)) ?
  return {width_.val(), width_.uncert_value};
}

UncertainDouble Hypermet::fwhm_energy(const HCalibration& cal) const
{
  double width = width_.val() * sqrt(log(2.0));
  double max_width = (width_.val() + width_.uncert_value) * sqrt(log(2.0));
  double min_width = (width_.val() - width_.uncert_value) * sqrt(log(2.0));

  double val = cal.channel_to_energy(position.val() + width)
      - cal.channel_to_energy(position.val() - width);
  double max = cal.channel_to_energy(position.val() + max_width)
      - cal.channel_to_energy(position.val() - max_width);
  double min = cal.channel_to_energy(position.val() + min_width)
      - cal.channel_to_energy(position.val() - min_width);

  return {val, 0.5 * (max - min)};
}

UncertainDouble Hypermet::fwhm_energy(const Calibration& cal) const
{
  double width = width_.val() * sqrt(log(2.0));
  double max_width = (width_.val() + width_.uncert_value) * sqrt(log(2.0));
  double min_width = (width_.val() - width_.uncert_value) * sqrt(log(2.0));

  double val = cal.transform(position.val() + width) - cal.transform(position.val() - width);
  double max = cal.transform(position.val() + max_width) - cal.transform(position.val() - max_width);
  double min = cal.transform(position.val() + min_width) - cal.transform(position.val() - min_width);

  return {val, 0.5 * (max - min)};
}


void Hypermet::update_indices(int32_t& i)
{
  amplitude.x_index = i++;
  position.x_index = i++;

  if (width_override)
    width_.x_index = i++;
  else
    width_.x_index = -1;

  if (short_tail.override && short_tail.enabled)
    short_tail.update_indices(i);

  if (right_tail.override && right_tail.enabled)
    right_tail.update_indices(i);

  if (long_tail.override && long_tail.enabled)
    long_tail.update_indices(i);

  if (step.override && step.enabled)
    step.update_indices(i);
}

void Hypermet::put(std::vector<double>& fit) const
{
  position.put(fit);
  amplitude.put(fit);
  width_.put(fit);
  short_tail.put(fit);
  right_tail.put(fit);
  long_tail.put(fit);
  step.put(fit);
}

void Hypermet::get(const std::vector<double>& fit)
{
  position.get(fit);
  amplitude.get(fit);
  width_.get(fit);
  short_tail.get(fit);
  right_tail.get(fit);
  long_tail.get(fit);
  step.get(fit);
}

void Hypermet::get_uncerts(const std::vector<double>& diagonals, double chisq_norm)
{
  chi_sq_norm = chisq_norm;
  position.get_uncert(diagonals, chisq_norm);
  amplitude.get_uncert(diagonals, chisq_norm);
  width_.get_uncert(diagonals, chisq_norm);
  short_tail.get_uncerts(diagonals, chisq_norm);
  right_tail.get_uncerts(diagonals, chisq_norm);
  long_tail.get_uncerts(diagonals, chisq_norm);
  step.get_uncerts(diagonals, chisq_norm);
}

PrecalcVals Hypermet::precalc_vals(double chan) const
{
  PrecalcVals ret;
  ret.ampl = amplitude.val();
  ret.half_ampl = 0.5 * ret.ampl;
  ret.width = width_.val();
  ret.spread = (chan - position.val()) / ret.width;
  return ret;
}

PrecalcVals Hypermet::precalc_vals_at(double chan, const std::vector<double>& fit) const
{
  PrecalcVals ret;
  ret.ampl = amplitude.val_at(fit[amplitude.x_index]);
  ret.half_ampl = 0.5 * ret.ampl;
  ret.width = width_.val_at(fit[width_.x_index]);
  ret.spread = (chan - position.val_at(fit[position.x_index])) / ret.width;
  return ret;
}

Hypermet::Components Hypermet::eval(double chan) const
{
  Hypermet::Components ret;

  auto pre = precalc_vals(chan);

  ret.gaussian = amplitude.val() * std::exp(-square(pre.spread));

  if (short_tail.enabled)
    ret.short_tail = short_tail.eval(pre);
  if (right_tail.enabled)
    ret.right_tail = right_tail.eval(pre);
  if (long_tail.enabled)
    ret.long_tail = long_tail.eval(pre);
  if (step.enabled)
    ret.step = step.eval(pre);
  return ret;
}

Hypermet::Components Hypermet::eval_at(double chan, const std::vector<double>& fit) const
{
  Hypermet::Components ret;

  auto pre = precalc_vals_at(chan, fit);

  ret.gaussian = amplitude.val_at(fit[amplitude.x_index]) * std::exp(-square(pre.spread));

  if (short_tail.enabled)
    ret.short_tail = short_tail.eval_at(pre, fit);
  if (right_tail.enabled)
    ret.right_tail = right_tail.eval_at(pre, fit);
  if (long_tail.enabled)
    ret.long_tail = long_tail.eval_at(pre, fit);
  if (step.enabled)
    ret.step = step.eval_at(pre, fit);
  return ret;
}

Hypermet::Components Hypermet::eval_grad(double chan, std::vector<double>& grads) const
{
  Hypermet::Components ret;

  auto pre = precalc_vals(chan);

  ret.gaussian = amplitude.val() * std::exp(-square(pre.spread));

  grads[width_.x_index] += ret.gaussian * 2.0 * square(pre.spread) / pre.width;
  grads[position.x_index] += ret.gaussian * 2.0 * pre.spread / pre.width;
  grads[amplitude.x_index] += ret.gaussian / pre.ampl;

  if (short_tail.enabled)
    ret.short_tail = short_tail.eval_grad(pre, grads,
                                          width_.x_index, position.x_index, amplitude.x_index);

  if (right_tail.enabled)
    ret.right_tail = short_tail.eval_grad(pre, grads,
                                          width_.x_index, position.x_index, amplitude.x_index);

  if (long_tail.enabled)
    ret.long_tail = long_tail.eval_grad(pre, grads,
                                        width_.x_index, position.x_index, amplitude.x_index);

  if (step.enabled)
    ret.step = step.eval_grad(pre, grads, width_.x_index, position.x_index, amplitude.x_index);

  grads[width_.x_index] *= width_.grad();
  grads[amplitude.x_index] *= amplitude.grad();
  grads[position.x_index] *= position.grad();

  return ret;
}

Hypermet::Components Hypermet::eval_grad_at(double chan, const std::vector<double>& fit,
                                            std::vector<double>& grads) const
{
  Hypermet::Components ret;

  auto pre = precalc_vals_at(chan, fit);

  ret.gaussian = amplitude.val_at(fit[amplitude.x_index]) * std::exp(-square(pre.spread));

  grads[width_.x_index] += ret.gaussian * 2.0 * square(pre.spread) / pre.width;
  grads[position.x_index] += ret.gaussian * 2.0 * pre.spread / pre.width;
  grads[amplitude.x_index] += ret.gaussian / pre.ampl;

  if (short_tail.enabled)
    ret.short_tail = short_tail.eval_grad_at(pre, fit, grads,
                                             width_.x_index, position.x_index, amplitude.x_index);

  if (right_tail.enabled)
    ret.right_tail = short_tail.eval_grad_at(pre, fit, grads,
                                             width_.x_index, position.x_index, amplitude.x_index);

  if (long_tail.enabled)
    ret.long_tail = long_tail.eval_grad_at(pre, fit, grads,
                                           width_.x_index, position.x_index, amplitude.x_index);

  if (step.enabled)
    ret.step = step.eval_grad_at(pre, fit, grads,
                                 width_.x_index, position.x_index, amplitude.x_index);

  grads[width_.x_index] *= width_.grad_at(fit[width_.x_index]);
  grads[amplitude.x_index] *= amplitude.grad_at(fit[amplitude.x_index]);
  grads[position.x_index] *= position.grad_at(fit[position.x_index]);

  return ret;
}

std::string Hypermet::to_string() const
{
  std::stringstream ss;
  ss << "pos = " << position.to_string() << "\n";
  ss << "amp = " << amplitude.to_string() << "\n";
  ss << "width" << (width_override ? "(OVERRIDEN) = " : " = ")
     << position.to_string() << "\n";
  ss << "left_skew = " << short_tail.to_string() << "\n";
  ss << "right_skew = " << right_tail.to_string() << "\n";
  ss << "long_tail = " << long_tail.to_string() << "\n";
  ss << "step = " << step.to_string() << "\n";
  return ss.str();
}

void to_json(nlohmann::json& j, const Hypermet& s)
{
  j["position"] = s.position;
  j["amplitude"] = s.amplitude;
  j["width_override"] = s.width_override;
  j["width"] = s.width_;
  j["short_tail"] = s.short_tail;
  j["right_tail"] = s.right_tail;
  j["long_tail"] = s.long_tail;
  j["step"] = s.step;
}

void from_json(const nlohmann::json& j, Hypermet& s)
{
  s.position = j["position"];
  s.amplitude = j["amplitude"];
  s.width_override = j["width_override"];
  s.width_ = j["width"];
  s.short_tail = j["short_tail"];
  s.right_tail = j["right_tail"];
  s.long_tail = j["long_tail"];
  s.step = j["step"];
}

}

