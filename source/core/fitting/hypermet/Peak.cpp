#include <core/fitting/hypermet/Peak.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

double Peak::Components::peak_skews() const
{
  return gaussian + short_tail + right_tail;
}

double Peak::Components::step_tail() const
{
  return long_tail + step;
}

double Peak::Components::all() const
{
  return peak_skews() + step_tail();
}

Peak::Peak()
{
  //amplitude.bound(0, 1000);

  width.bound(0.8, 5.0);

  short_tail.amplitude.bound(0.02, 1.5);
  short_tail.slope.bound(0.2, 0.5);

  right_tail.amplitude.bound(0.01, 0.9);
  right_tail.slope.bound(0.3, 1.5);

  long_tail.amplitude.bound(0.0001, 0.15);
  long_tail.slope.bound(2.5, 50);

  step.amplitude.bound(0.000001, 0.05);
}

void Peak::apply_defaults(const Peak& other)
{
  if (!width_override)
    width = other.width;
  if (!short_tail.override)
    short_tail = other.short_tail;
  if (!right_tail.override)
    right_tail = other.right_tail;
  if (!long_tail.override)
    long_tail = other.long_tail;
  if (!step.override)
    step = other.step;
}

void Peak::force_defaults(const Peak& other)
{
  width = other.width;
  short_tail = other.short_tail;
  right_tail = other.right_tail;
  long_tail = other.long_tail;
  step = other.step;
}

Peak Peak::gaussian_only() const
{
  Peak ret = *this;
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

bool Peak::is_gaussian_only() const
{
  return (!short_tail.enabled && !right_tail.enabled && !long_tail.enabled && !step.enabled);
}

bool Peak::sanity_check(double min_x, double max_x) const
{
  auto amp = amplitude.val();
  auto wid = width.val();
  auto pos = position.val();
  return std::isfinite(amp) && (amp > 0.0) &&
      std::isfinite(wid) && (wid > 0.0) &&
      std::isfinite(pos) && (min_x < pos) && (pos < max_x);
}

bool Peak::full_energy_peak() const
{
  return (step.side == Side::left);
}

void Peak::full_energy_peak(bool flag)
{
  if (flag)
    step.side = Side::left;
  else
    step.side = Side::right;
}

bool Peak::operator<(const Peak& other) const
{
  return position.val() < other.position.val();
}

double Peak::id() const
{
  return position.val();
}

UncertainDouble Peak::peak_position() const
{
  return {position.val(), position.uncert()};
}

UncertainDouble Peak::peak_energy(const HCalibration& cal) const
{
  // \todo if there is a curve, does the slope not differ with x?
  return {cal.channel_to_energy(position.val()),
          cal.energy_slope() * position.uncert()};
}

UncertainDouble Peak::peak_energy(const Calibration& cal) const
{
  return {cal.transform(position.val()),
          cal.function()->derivative(position.val()) * position.uncert()};
}

UncertainDouble Peak::area() const
{
  auto a = amplitude.val() * width.val() * (std::sqrt(M_PI) +
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
//  // * Hessinv.coeff(DEL.index(), DEL.index());
//  //(dGAM/dX*dArea/dGAM)*(dDEL/dY*dArea/dDEL)*Covar(X,Y)*Chisq
//  t += square(amplitude.grad() * width_.val() * a) * cs;
//  // * Hessinv.coeff(Peak(PeakIndex).GAM.index(), DEL.index())
//  if (short_tail.amplitude.to_fit)
//  {
//    //(dGAM/dX*dArea/dGAM)*(dAST/dY*dArea/dAST)*Covar(X,Y)*Chisq
//    t += (amplitude.grad() * width_.val() * a) *
//        (short_tail.amplitude.grad() * amplitude.val() * width_.val() *
//            short_tail.slope.val() * std::exp(-0.25 / square(short_tail.slope.val()))) * cs;
//    // * Hessinv.coeff(Peak(PeakIndex).GAM.index(), AST.index())
//  }
//  if (short_tail.slope.to_fit)
//  {
//    // (dGAM/dX*dArea/dGAM)*(dBST/dY*dArea/dBST)*Covar(X,Y)*Chisq
//    t += (amplitude.grad() * width_.val() * a) *
//        (short_tail.slope.grad() * amplitude.val() * width_.val()
//            * short_tail.amplitude.val() * (1 + 0.5 / square(short_tail.slope.val()))
//            * std::exp(-0.25 / square(short_tail.slope.val()))) * cs;
//    // * Hessinv.coeff(Peak(PeakIndex).GAM.index(), AST.index())
//  }
  return {a, std::sqrt(a) * std::max(1.0, chi_sq_norm)};
}

UncertainDouble Peak::peak_area_eff(const HCalibration& cal) const
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

UncertainDouble Peak::fwhm() const
{
  // \todo multiply with sqrt(log(2.0)) ?
  return {width.val(), width.uncert()};
}

UncertainDouble Peak::fwhm_energy(const HCalibration& cal) const
{
  double bin_width = width.val() * sqrt(log(2.0));
  double max_width = (width.val() + width.uncert()) * sqrt(log(2.0));
  double min_width = (width.val() - width.uncert()) * sqrt(log(2.0));

  double val = cal.channel_to_energy(position.val() + bin_width)
      - cal.channel_to_energy(position.val() - bin_width);
  double max = cal.channel_to_energy(position.val() + max_width)
      - cal.channel_to_energy(position.val() - max_width);
  double min = cal.channel_to_energy(position.val() + min_width)
      - cal.channel_to_energy(position.val() - min_width);

  return {val, 0.5 * (max - min)};
}

UncertainDouble Peak::fwhm_energy(const Calibration& cal) const
{
  double bin_width = width.val() * sqrt(log(2.0));
  double max_width = (width.val() + width.uncert()) * sqrt(log(2.0));
  double min_width = (width.val() - width.uncert()) * sqrt(log(2.0));

  double val = cal.transform(position.val() + bin_width) - cal.transform(position.val() - bin_width);
  double max = cal.transform(position.val() + max_width) - cal.transform(position.val() - max_width);
  double min = cal.transform(position.val() + min_width) - cal.transform(position.val() - min_width);

  return {val, 0.5 * (max - min)};
}


void Peak::update_indices(int32_t& i)
{
  position.update_index(i);
  amplitude.update_index(i);

  if (width_override)
    width.update_index(i);

  if (short_tail.override)
    short_tail.update_indices(i);

  if (right_tail.override)
    right_tail.update_indices(i);

  if (long_tail.override)
    long_tail.update_indices(i);

  if (step.override)
    step.update_indices(i);
}

void Peak::put(Eigen::VectorXd& fit) const
{
  position.put(fit);
  amplitude.put(fit);
  width.put(fit);
  short_tail.put(fit);
  right_tail.put(fit);
  long_tail.put(fit);
  step.put(fit);
}

void Peak::get(const Eigen::VectorXd& fit)
{
  position.get(fit);
  amplitude.get(fit);
  width.get(fit);
  short_tail.get(fit);
  right_tail.get(fit);
  long_tail.get(fit);
  step.get(fit);
}

void Peak::get_uncerts(const Eigen::VectorXd& diagonals, double chisq_norm)
{
  chi_sq_norm = chisq_norm;
  position.get_uncert(diagonals, chisq_norm);
  amplitude.get_uncert(diagonals, chisq_norm);
  width.get_uncert(diagonals, chisq_norm);
  short_tail.get_uncerts(diagonals, chisq_norm);
  right_tail.get_uncerts(diagonals, chisq_norm);
  long_tail.get_uncerts(diagonals, chisq_norm);
  step.get_uncerts(diagonals, chisq_norm);
}

PrecalcVals Peak::precalc_vals(double chan) const
{
  PrecalcVals ret;
  ret.ampl = amplitude.val();
  ret.half_ampl = 0.5 * ret.ampl;
  ret.width = width.val();
  ret.spread = (chan - position.val()) / ret.width;

  ret.amp_grad = amplitude.grad();
  ret.width_grad = width.grad();
  ret.pos_grad = position.grad();

  ret.i_amp = amplitude.index();
  ret.i_width = width.index();
  ret.i_pos = position.index();
  return ret;
}

PrecalcVals Peak::precalc_vals_at(double chan, const Eigen::VectorXd& fit) const
{
  PrecalcVals ret;
  ret.ampl = amplitude.val_from(fit);
  ret.half_ampl = 0.5 * ret.ampl;
  ret.width = width.val_from(fit);
  ret.spread = (chan - position.val_from(fit)) / ret.width;

  ret.amp_grad = amplitude.grad_from(fit);
  ret.width_grad = width.grad_from(fit);
  ret.pos_grad = position.grad_from(fit);

  ret.i_amp = amplitude.index();
  ret.i_width = width.index();
  ret.i_pos = position.index();
  return ret;
}

Peak::Components Peak::eval(double chan) const
{
  Peak::Components ret;

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

Peak::Components Peak::eval_at(double chan, const Eigen::VectorXd& fit) const
{
  Peak::Components ret;

  auto pre = precalc_vals_at(chan, fit);

  ret.gaussian = amplitude.val_from(fit) * std::exp(-square(pre.spread));

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

Peak::Components Peak::eval_grad(double chan, Eigen::VectorXd& grads) const
{
  Peak::Components ret;

  auto pre = precalc_vals(chan);

  ret.gaussian = amplitude.val() * std::exp(-square(pre.spread));

  if (width.to_fit)
    grads[width.index()] += pre.width_grad * ret.gaussian * 2.0 * square(pre.spread) / pre.width;
  if (position.to_fit)
    grads[position.index()] += pre.pos_grad * (ret.gaussian * 2.0 * pre.spread / pre.width);
  if (amplitude.to_fit)
    grads[amplitude.index()] += pre.amp_grad * ret.gaussian / pre.ampl;

  if (short_tail.enabled)
    ret.short_tail = short_tail.eval_grad(pre, grads);
  if (right_tail.enabled)
    ret.right_tail = short_tail.eval_grad(pre, grads);
  if (long_tail.enabled)
    ret.long_tail = long_tail.eval_grad(pre, grads);
  if (step.enabled)
    ret.step = step.eval_grad(pre, grads);

  return ret;
}

Peak::Components Peak::eval_grad_at(double chan, const Eigen::VectorXd& fit,
                                            Eigen::VectorXd& grads) const
{
  Peak::Components ret;

  auto pre = precalc_vals_at(chan, fit);

  ret.gaussian = amplitude.val_from(fit) * std::exp(-square(pre.spread));

  if (width.to_fit)
    grads[width.index()] += pre.width_grad * ret.gaussian * 2.0 * square(pre.spread) / pre.width;
  if (position.to_fit)
    grads[position.index()] += pre.pos_grad * (ret.gaussian * 2.0 * pre.spread / pre.width);
  if (amplitude.to_fit)
    grads[amplitude.index()] += pre.amp_grad * ret.gaussian / pre.ampl;

  if (short_tail.enabled)
    ret.short_tail = short_tail.eval_grad_at(pre, fit, grads);
  if (right_tail.enabled)
    ret.right_tail = short_tail.eval_grad_at(pre, fit, grads);
  if (long_tail.enabled)
    ret.long_tail = long_tail.eval_grad_at(pre, fit, grads);
  if (step.enabled)
    ret.step = step.eval_grad_at(pre, fit, grads);

  return ret;
}

std::string Peak::to_string(std::string prepend) const
{
  std::stringstream ss;
  ss << prepend << "pos = " << position.to_string() << "\n";
  ss << prepend << "amp = " << amplitude.to_string() << "\n";
  ss << prepend << "width"
     << (width_override ? "(OVERRIDEN) = " : " = ")
     << width.to_string() << "\n";
  ss << prepend << "left_skew = " << short_tail.to_string() << "\n";
  ss << prepend << "right_skew = " << right_tail.to_string() << "\n";
  ss << prepend << "long_tail = " << long_tail.to_string() << "\n";
  ss << prepend << "step = " << step.to_string() << "\n";
  return ss.str();
}

void to_json(nlohmann::json& j, const Peak& s)
{
  j["position"] = s.position;
  j["amplitude"] = s.amplitude;
  j["width_override"] = s.width_override;
  j["width"] = s.width;
  j["short_tail"] = s.short_tail;
  j["right_tail"] = s.right_tail;
  j["long_tail"] = s.long_tail;
  j["step"] = s.step;
  j["sum4"] = s.sum4;

  // \todo should this be here?
  j["chi_sq_norm"] = s.chi_sq_norm;
}

void from_json(const nlohmann::json& j, Peak& s)
{
  s.position = j["position"];
  s.amplitude = j["amplitude"];
  s.width_override = j["width_override"];
  s.width = j["width"];
  s.short_tail = j["short_tail"];
  s.right_tail = j["right_tail"];
  s.long_tail = j["long_tail"];
  s.step = j["step"];
  s.sum4 = j["sum4"];

  // \todo should this be here?
  s.chi_sq_norm = j["chi_sq_norm"];
}

}

