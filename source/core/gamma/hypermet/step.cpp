#include <core/gamma/hypermet/step.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void Step::reset_indices()
{
  amplitude.reset_index();
}

void Step::update_indices(int32_t& i)
{
  if (enabled)
    amplitude.update_index(i);
  else
    reset_indices();
}

void Step::put(Eigen::VectorXd& fit) const
{
  amplitude.put(fit);
}

void Step::get(const Eigen::VectorXd& fit)
{
  amplitude.get(fit);
}

void Step::get_uncerts(const Eigen::VectorXd& diagonals, double chisq_norm)
{
  amplitude.get_uncert(diagonals, chisq_norm);
}

double Step::eval_with(const PrecalcVals& pre, double ampl) const
{
  return pre.half_ampl * ampl * std::erfc(flip(side, pre.spread));
}

double Step::eval(const PrecalcVals& pre) const
{
  return eval_with(pre, amplitude.val());
}

double Step::eval_at(const PrecalcVals& pre, const Eigen::VectorXd& fit) const
{
  return eval_with(pre, amplitude.val_from(fit));
}

double Step::eval_grad(const PrecalcVals& pre, Eigen::VectorXd& grads) const
{
  double ampl = amplitude.val();
  double ret = eval_with(pre, ampl);

  if (pre.i_width > AbstractParam::InvalidIndex)
    grads[pre.i_width] += pre.width_grad * (pre.ampl * flip(side, ampl) / std::sqrt(M_PI) *
        std::exp(-square(pre.spread)) * pre.spread / pre.width);
  if (pre.i_amp > AbstractParam::InvalidIndex)
    grads[pre.i_amp] += pre.amp_grad * ret / pre.ampl;
  // \todo pos unused?

  if (amplitude.valid_index())
    grads[amplitude.index()] += ret / ampl * amplitude.grad();

  return ret;
}

double Step::eval_grad_at(const PrecalcVals& pre, const Eigen::VectorXd& fit,
                          Eigen::VectorXd& grads) const
{
  double ampl = amplitude.val_from(fit);
  double ret = eval_with(pre, ampl);

  if (pre.i_width > AbstractParam::InvalidIndex)
    grads[pre.i_width] += pre.width_grad * (pre.ampl * flip(side, ampl) / std::sqrt(M_PI) *
        std::exp(-square(pre.spread)) * pre.spread / pre.width);
  if (pre.i_amp > AbstractParam::InvalidIndex)
    grads[pre.i_amp] += pre.amp_grad * ret / pre.ampl;
  // \todo pos unused?

  if (amplitude.valid_index())
    grads[amplitude.index()] += ret / ampl * amplitude.grad_from(fit);

  return ret;
}

bool Step::sane(double amp_min_epsilon, double amp_max_epsilon) const
{
  if (amplitude.to_fit && amplitude.at_extremum(amp_min_epsilon, amp_max_epsilon))
    return false;
  return true;
}

std::string Step::to_string() const
{
  return fmt::format("{}{:<9} {:<5}  amp={}",
                     enabled ? "ON " : "OFF",
                     override ? " OVERRIDE" : "",
                     DAQuiri::to_string(side),
                     amplitude.to_string());

}

void to_json(nlohmann::json& j, const Step& s)
{
  j["enabled"] = s.enabled;
  j["override"] = s.override;
  j["side"] = to_string(s.side);
  j["amplitude"] = s.amplitude;
}

void from_json(const nlohmann::json& j, Step& s)
{
  s.enabled = j["enabled"];
  s.override = j["override"];
  s.side = to_side(j["side"]);
  s.amplitude = j["amplitude"];
}

}
