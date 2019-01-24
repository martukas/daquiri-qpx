#include <optimizerBFGS/Tail.h>

#include <core/util/custom_logger.h>
#include <optimizerBFGS/more_math.h>

namespace Hypermet
{

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

double Tail::eval_with(const PrecalcVals& pre, double ampl, double slp) const
{
  return pre.half_ampl * ampl *
      std::exp(flip(pre.spread) / slp) * std::erfc(0.5 / slp + flip(pre.spread));
}

double Tail::eval(const PrecalcVals& pre) const
{
  return eval_with(pre, amplitude.val(), slope.val());
}

double Tail::eval_at(const PrecalcVals& pre, const std::vector<double> fit) const
{
  return eval_with(pre,
                   amplitude.val_at(fit[amplitude.x_index]),
                   slope.val_at(fit[slope.x_index])
  );
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
  return ret;
}

double Tail::eval_grad_at(const PrecalcVals& pre, const std::vector<double> fit,
                          std::vector<double>& grads,
                          size_t i_width, size_t i_pos, size_t i_amp) const
{
  double ampl = amplitude.val_at(fit[amplitude.x_index]);
  double slp = slope.val_at(fit[slope.x_index]);
  double ret = eval_with(pre, ampl, slp);
  double spread = flip(pre.spread);
  double t2 = (pre.ampl * ampl * std::exp(spread / slp) / std::sqrt(M_PI) *
      std::exp(-1.0 * square(1.0 / (2.0 * slp) + spread)) / pre.width);
  grads[i_width] += -spread / (pre.width * slp) * ret + t2 * spread;
  grads[i_pos] += -1.0 / (slp * pre.width) * ret + t2;
  grads[i_amp] += ret / ampl;

  if (amplitude.to_fit)
    grads[amplitude.x_index] += ret / ampl * amplitude.grad_at(fit[amplitude.x_index]);
  if (slope.to_fit)
    grads[slope.x_index] += slope.grad_at(fit[slope.x_index]) * ((-spread / square(slp)) *
        ret + (pre.width / (2.0 * square(slp)) * t2));
  return ret;
}



}
