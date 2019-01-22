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

}
