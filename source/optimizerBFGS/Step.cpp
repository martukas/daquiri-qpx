#include <optimizerBFGS/Step.h>

#include <core/util/custom_logger.h>
#include <optimizerBFGS/more_math.h>

namespace Hypermet
{

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


}
