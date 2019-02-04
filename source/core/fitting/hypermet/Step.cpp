#include <core/fitting/hypermet/Step.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void Step::update_indices(int32_t& i)
{
  if (amplitude.to_fit)
    amplitude.x_index = i++;
  else
    amplitude.x_index = -1;
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
  return pre.half_ampl * ampl * std::erfc(flip(pre.spread));
}

double Step::eval(const PrecalcVals& pre) const
{
  return eval_with(pre, amplitude.val());
}

double Step::eval_at(const PrecalcVals& pre, const Eigen::VectorXd& fit) const
{
  return eval_with(pre, amplitude.val_from(fit));
}

double Step::eval_grad(const PrecalcVals& pre, Eigen::VectorXd& grads,
                       size_t i_width, size_t i_pos, size_t i_amp) const
{
  double ampl = amplitude.val();
  double ret = eval_with(pre, ampl);

  grads[i_width] += pre.width_grad * (pre.ampl * flip(ampl) / std::sqrt(M_PI) *
      std::exp(-square(pre.spread)) * pre.spread / pre.width);
  grads[i_amp] += pre.pos_grad * ret / pre.ampl;
  if (amplitude.to_fit)
    grads[amplitude.x_index] += ret / ampl * amplitude.grad();

  // \todo pos unused?
  return ret;
}

double Step::eval_grad_at(const PrecalcVals& pre, const Eigen::VectorXd& fit,
                          Eigen::VectorXd& grads,
                          size_t i_width, size_t i_pos, size_t i_amp) const
{
  double ampl = amplitude.val_from(fit);
  double ret = eval_with(pre, ampl);

  grads[i_width] += pre.width_grad * (pre.ampl * flip(ampl) / std::sqrt(M_PI) *
      std::exp(-square(pre.spread)) * pre.spread / pre.width);
  grads[i_amp] += pre.pos_grad * ret / pre.ampl;
  if (amplitude.to_fit)
    grads[amplitude.x_index] += ret / ampl * amplitude.grad_from(fit);

  // \todo pos unused?
  return ret;
}

std::string Step::to_string() const
{
  return fmt::format("{}{} {}  amp={}",
                     enabled ? "ON" : "OFF",
                     override ? " OVERRIDE" : "",
                     side_to_string(side),
                     amplitude.to_string());

}

void to_json(nlohmann::json& j, const Step& s)
{
  j["enabled"] = s.enabled;
  j["override"] = s.override;
  j["side"] = side_to_string(s.side);
  j["amplitude"] = s.amplitude;
}

void from_json(const nlohmann::json& j, Step& s)
{
  s.enabled = j["enabled"];
  s.override = j["override"];
  s.side = side_from_string(j["side"]);
  s.amplitude = j["amplitude"];
}

}
