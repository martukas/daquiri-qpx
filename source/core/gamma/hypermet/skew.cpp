#include <core/gamma/hypermet/skew.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void Tail::reset_indices()
{
  amplitude.reset_index();
  slope.reset_index();
}

void Tail::update_indices(int32_t& i)
{
  if (enabled)
  {
    amplitude.update_index(i);
    slope.update_index(i);
  }
  else
    reset_indices();
}

void Tail::put(Eigen::VectorXd& fit) const
{
  amplitude.put(fit);
  slope.put(fit);
}

void Tail::get(const Eigen::VectorXd& fit)
{
  amplitude.get(fit);
  slope.get(fit);
}

void Tail::get_uncerts(const Eigen::VectorXd& diagonals, double chisq_norm)
{
  amplitude.get_uncert(diagonals, chisq_norm);
  slope.get_uncert(diagonals, chisq_norm);
}

double Tail::eval_with(const PrecalcVals& pre, double ampl, double slp) const
{
  // \todo make this param:
  double spread = flip(side, pre.spread);
  return pre.half_ampl * ampl * std::exp(spread / slp) * std::erfc(0.5 / slp + spread);
}

double Tail::eval(const PrecalcVals& pre) const
{
  return eval_with(pre, amplitude.val(), slope.val());
}

double Tail::eval_at(const PrecalcVals& pre, const Eigen::VectorXd& fit) const
{
  return eval_with(pre, amplitude.val_from(fit), slope.val_from(fit));
}

double Tail::eval_grad(const PrecalcVals& pre, Eigen::VectorXd& grads) const
{
  double ampl = amplitude.val();
  double slp = slope.val();
  double ret = eval_with(pre, ampl, slp);
  double spread = flip(side, pre.spread);
  double t2 = (pre.ampl * ampl * std::exp(spread / slp) / std::sqrt(M_PI) *
      std::exp(-1.0 * square(1.0 / (2.0 * slp) + spread)) / pre.width);
  if (pre.i_width > AbstractValue::InvalidIndex)
    grads[pre.i_width] += pre.width_grad * spread * (t2  - ret / (pre.width * slp));
  if (pre.i_pos > AbstractValue::InvalidIndex)
    grads[pre.i_pos] += pre.pos_grad * (-ret / (slp * pre.width) + t2);
  if (pre.i_amp > AbstractValue::InvalidIndex)
    grads[pre.i_amp] += pre.amp_grad * ret / pre.ampl;

  if (amplitude.valid_index())
    grads[amplitude.index()] += amplitude.grad() * ret / ampl;
  if (slope.valid_index())
    grads[slope.index()] += slope.grad() * ((-spread / square(slp)) *
        ret + (pre.width / (2.0 * square(slp)) * t2));
  return ret;
}

double Tail::eval_grad_at(const PrecalcVals& pre, const Eigen::VectorXd& fit,
                          Eigen::VectorXd& grads) const
{
  double ampl = amplitude.val_from(fit);
  double slp = slope.val_from(fit);
  double ret = eval_with(pre, ampl, slp);
  double spread = flip(side, pre.spread);
  double t2 = (pre.ampl * ampl * std::exp(spread / slp) / std::sqrt(M_PI) *
      std::exp(-1.0 * square(1.0 / (2.0 * slp) + spread)) / pre.width);
  if (pre.i_width > AbstractValue::InvalidIndex)
    grads[pre.i_width] += pre.width_grad * spread * (t2 - ret / (pre.width * slp));
  if (pre.i_pos > AbstractValue::InvalidIndex)
    grads[pre.i_pos] += pre.pos_grad * (-ret / (slp * pre.width) + t2);
  if (pre.i_amp > AbstractValue::InvalidIndex)
    grads[pre.i_amp] += pre.amp_grad * ret / pre.ampl;

  if (amplitude.valid_index())
    grads[amplitude.index()] += amplitude.grad_from(fit) * ret / ampl;
  if (slope.valid_index())
    grads[slope.index()] += slope.grad_from(fit) * ((-spread / square(slp)) *
        ret + (pre.width / (2.0 * square(slp)) * t2));
  return ret;
}

bool Tail::sane(double amp_min_epsilon, double amp_max_epsilon, double slope_epsilon) const
{
  if (amplitude.to_fit && amplitude.at_extremum(amp_min_epsilon, amp_max_epsilon))
    return false;
  if (slope.to_fit && slope.at_extremum(slope_epsilon, slope_epsilon))
    return false;
  return true;
}

std::string Tail::to_string() const
{
  return fmt::format("{}{:<9} {:<5}  amp={}  slope={}",
                     enabled ? "ON " : "OFF",
                     override ? " OVERRIDE" : "",
                     DAQuiri::to_string(side),
                     amplitude.to_string(),
                     slope.to_string());

}

void to_json(nlohmann::json& j, const Tail& s)
{
  j["enabled"] = s.enabled;
  j["override"] = s.override;
  j["side"] = to_string(s.side);
  j["amplitude"] = s.amplitude;
  j["slope"] = s.slope;
}

void from_json(const nlohmann::json& j, Tail& s)
{
  s.enabled = j["enabled"];
  s.override = j["override"];
  s.side = to_side(j["side"]);
  s.amplitude = j["amplitude"];
  s.slope = j["slope"];
}

}