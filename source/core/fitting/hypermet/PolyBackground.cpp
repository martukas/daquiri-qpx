#include <core/fitting/hypermet/PolyBackground.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void PolyBackground::update_indices(int32_t& i)
{
  background_base_.x_index = i++;

  if (slope_enabled_)
    background_slope_.x_index = i++;
  else
    background_slope_.x_index = -1;

  if (curve_enabled_)
    background_curve_.x_index = i++;
  else
    background_curve_.x_index = -1;
}

void PolyBackground::put(std::vector<double>& fit) const
{
  background_base_.put(fit);
  background_slope_.put(fit);
  background_curve_.put(fit);
}

void PolyBackground::get(const std::vector<double>& fit)
{
  background_base_.get(fit);
  background_slope_.get(fit);
  background_curve_.get(fit);
}

void PolyBackground::get_uncerts(const std::vector<double>& diagonals, double chisq_norm)
{
  background_base_.get_uncert(diagonals, chisq_norm);
  background_slope_.get_uncert(diagonals, chisq_norm);
  background_curve_.get_uncert(diagonals, chisq_norm);
}

double PolyBackground::eval(double bin) const
{
  double ret = background_base_.val();
  if (slope_enabled_)
    ret += background_slope_.val() * (bin - bin_offset);
  if (curve_enabled_)
    ret += background_curve_.val() * square(bin - bin_offset);
  return ret;
}

double PolyBackground::eval_at(double bin, const std::vector<double>& fit) const
{
  double ret = background_base_.val_at(fit[background_base_.x_index]);
  if (slope_enabled_)
    ret += background_slope_.val_at(fit[background_slope_.x_index]) * (bin - bin_offset);
  if (curve_enabled_)
    ret += background_curve_.val_at(fit[background_curve_.x_index]) * square(bin - bin_offset);
  return ret;
}

double PolyBackground::eval_grad(double bin, std::vector<double>& gradients) const
{
  double ret = background_base_.val();
  gradients[background_base_.x_index] = background_base_.grad();
  if (slope_enabled_)
  {
    ret += background_slope_.val() * (bin - bin_offset);
    gradients[background_slope_.x_index] = (bin - bin_offset);
  }

  if (curve_enabled_)
  {
    ret += background_curve_.val() * square(bin - bin_offset);
    gradients[background_curve_.x_index] = square(bin - bin_offset);
  }
  return ret;
}

double PolyBackground::eval_grad_at(double bin,
                                    const std::vector<double>& fit,
                                    std::vector<double>& gradients) const
{
  double ret = background_base_.val_at(fit[background_base_.x_index]);
  gradients[background_base_.x_index] = background_base_.grad_at(fit[background_base_.x_index]);
  if (slope_enabled_)
  {
    ret += background_slope_.val_at(fit[background_slope_.x_index]) * (bin - bin_offset);
    gradients[background_slope_.x_index] = (bin - bin_offset);
  }

  if (curve_enabled_)
  {
    ret += background_curve_.val_at(fit[background_curve_.x_index]) * square(bin - bin_offset);
    gradients[background_curve_.x_index] = square(bin - bin_offset);
  }
  return ret;
}

std::string PolyBackground::to_string() const
{
  std::stringstream ss;
  ss << "x=bin-" << bin_offset << "    ";
  ss << "base=" << background_base_.to_string();
  if (slope_enabled_)
    ss << "   base=" << background_slope_.to_string();
  if (curve_enabled_)
    ss << "   base=" << background_curve_.to_string();
  return ss.str();
}

void to_json(nlohmann::json& j, const PolyBackground& s)
{
  j["bin_offset"] = s.bin_offset;
  j["base"] = s.background_base_;
  if (s.slope_enabled_)
    j["slope"] = s.background_slope_;
  if (s.curve_enabled_)
    j["curve"] = s.background_curve_;
}

void from_json(const nlohmann::json& j, PolyBackground& s)
{
  s.bin_offset = j["bin_offset"];
  s.background_base_ = j["base"];
  if (j.count("slope"))
    s.background_slope_ = j["slope"];
  if (j.count("curve"))
    s.background_curve_ = j["curve"];
}

}
