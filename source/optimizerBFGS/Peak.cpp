#include <optimizerBFGS/Peak.h>

#include <core/util/custom_logger.h>
#include <optimizerBFGS/more_math.h>

namespace Hypermet
{

double Tail::eval_with(const PrecalcVals& pre, double ampl, double slp) const
{
  return pre.half_ampl * ampl *
      std::exp(flip(pre.spread) / slp) * std::erfc(0.5 / slp + flip(pre.spread));
}

double Tail::eval(const PrecalcVals& pre) const
{
  return eval_with(pre, amplitude.val(), slope.val());
}

double Tail::eval_grad(const PrecalcVals& pre, std::vector<double>& grads,
                       size_t i_width, size_t i_pos, size_t i_amp) const
{
  double ampl = amplitude.val();
  double slp = slope.val();
  double ret = eval_with(pre, ampl, slp);
  double spread = flip(pre.spread);
  double t2 = (pre.ampl * ampl * std::exp(spread / slp) / std::sqrt(M_PI) *
      std::exp(-1.0 * square(1.0 / (2.0 * slp) + spread)) / pre.width);
  grads[i_width] += -spread / (pre.width * slp) * ret + t2 * spread;
  grads[i_pos] += -1.0 / (slp * pre.width) * ret + t2;
  grads[i_amp] += ret / ampl;

  if (amplitude.to_fit)
    grads[amplitude.x_index] += ret / ampl * amplitude.grad();
  if (slope.to_fit)
    grads[slope.x_index] += slope.grad() * ((-spread / square(slp)) *
        ret + (pre.width / (2.0 * square(slp)) * t2));
}

void Tail::update_indices(int32_t& i)
{
  if (amplitude.to_fit)
    amplitude.x_index = i++;
  else
    amplitude.x_index = -1;

  if (slope.to_fit)
    slope.x_index = i++;
  else
    slope.x_index = -1;
}

void Tail::put(std::vector<double>& fit) const
{
  amplitude.put(fit);
  slope.put(fit);
}

void Tail::get(const std::vector<double>& fit)
{
  amplitude.get(fit);
  slope.get(fit);
}

void Tail::get_uncerts(const std::vector<double>& diagonals, double chisq_norm)
{
  amplitude.get_uncert(diagonals, chisq_norm);
  slope.get_uncert(diagonals, chisq_norm);
}





double Step::eval_with(const PrecalcVals& pre, double ampl) const
{
  return pre.half_ampl * amplitude.val() * std::erfc(flip(pre.spread));

}

double Step::eval(const PrecalcVals& pre) const
{
  return eval_with(pre, amplitude.val());
}

double Step::eval_grad(const PrecalcVals& pre, std::vector<double>& grads,
                       size_t i_width, size_t i_pos, size_t i_amp) const
{
  double ampl = amplitude.val();
  double ret = eval_with(pre, ampl);

  grads[i_width] += (pre.ampl * flip(ampl) / std::sqrt(M_PI) *
      std::exp(-square(pre.spread)) * pre.spread / pre.width);
  grads[i_amp] += ret / pre.ampl;
  if (amplitude.to_fit)
    grads[amplitude.x_index] += ret / ampl * amplitude.grad();
}

void Step::update_indices(int32_t& i)
{
  if (amplitude.to_fit)
    amplitude.x_index = i++;
  else
    amplitude.x_index = -1;
}

void Step::put(std::vector<double>& fit) const
{
  amplitude.put(fit);
}

void Step::get(const std::vector<double>& fit)
{
  amplitude.get(fit);
}

void Step::get_uncerts(const std::vector<double>& diagonals, double chisq_norm)
{
  amplitude.get_uncert(diagonals, chisq_norm);
}



int32_t Peak::step_type() const
{
  return static_cast<int32_t>(step.flip(1.0));
}

double Peak::peak_position() const
{
  return position.val();
}

double Peak::peak_position_unc() const
{
  return position.uncert_value;
}

double Peak::peak_energy(const Calibration& cal) const
{
  return cal.channel_to_energy(position.val());
}

double Peak::peak_energy_unc(const Calibration& cal) const
{
  return cal.energy_slope() * position.uncert_value;
}

bool Peak::full_energy_peak() const
{
  return (step.flip(1.0) > 0);
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

PrecalcVals Peak::precalc_vals(double chan) const
{
  PrecalcVals ret;
  ret.ampl = amplitude.val();
  ret.half_ampl = 0.5 * ret.ampl;
  ret.width = width_.val();
  ret.spread = (chan - position.val()) / ret.width;
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

double Peak::area() const
{
  return amplitude.val() * width_.val() * (std::sqrt(M_PI) +
      short_tail.amplitude.val() * short_tail.slope.val() *
          std::exp(-0.25 / square(short_tail.slope.val())) +
      right_tail.amplitude.val() * right_tail.slope.val() *
          std::exp(-0.25 / square(right_tail.slope.val())));
}

double Peak::area_uncert(double chisq_norm) const
{
  // \todo make this more rigorous
  double t = area();
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
  return std::sqrt(t) * std::max(1.0, chisq_norm);
}

Peak::Components Peak::eval_grad(double chan, std::vector<double>& grads) const
{
  Peak::Components ret;

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

void Peak::update_indices(int32_t& i)
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

void Peak::put(std::vector<double>& fit) const
{
  position.put(fit);
  amplitude.put(fit);
  width_.put(fit);
  short_tail.put(fit);
  right_tail.put(fit);
  long_tail.put(fit);
  step.put(fit);
}

void Peak::get(const std::vector<double>& fit)
{
  position.get(fit);
  amplitude.get(fit);
  width_.get(fit);
  short_tail.get(fit);
  right_tail.get(fit);
  long_tail.get(fit);
  step.get(fit);
}

void Peak::get_uncerts(const std::vector<double>& diagonals, double chisq_norm)
{
  position.get_uncert(diagonals, chisq_norm);
  amplitude.get_uncert(diagonals, chisq_norm);
  width_.get_uncert(diagonals, chisq_norm);
  short_tail.get_uncerts(diagonals, chisq_norm);
  right_tail.get_uncerts(diagonals, chisq_norm);
  long_tail.get_uncerts(diagonals, chisq_norm);
  step.get_uncerts(diagonals, chisq_norm);
}



}
