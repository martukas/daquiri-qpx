#include <core/fitting/hypermet/PolyBackground.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void PolyBackground::update_indices(int32_t& i)
{
  base.x_index = i++;

  if (slope_enabled)
    slope.x_index = i++;
  else
    slope.x_index = -1;

  if (curve_enabled)
    curve.x_index = i++;
  else
    curve.x_index = -1;
}

void PolyBackground::put(std::vector<double>& fit) const
{
  base.put(fit);
  slope.put(fit);
  curve.put(fit);
}

void PolyBackground::get(const std::vector<double>& fit)
{
  base.get(fit);
  slope.get(fit);
  curve.get(fit);
}

void PolyBackground::get_uncerts(const std::vector<double>& diagonals, double chisq_norm)
{
  base.get_uncert(diagonals, chisq_norm);
  slope.get_uncert(diagonals, chisq_norm);
  curve.get_uncert(diagonals, chisq_norm);
}

double PolyBackground::eval(double bin) const
{
  double ret = base.val();
  if (slope_enabled)
    ret += slope.val() * (bin - x_offset);
  if (curve_enabled)
    ret += curve.val() * square(bin - x_offset);
  return ret;
}

double PolyBackground::eval_at(double bin, const std::vector<double>& fit) const
{
  double ret = base.val_at(fit[base.x_index]);
  if (slope_enabled)
    ret += slope.val_at(fit[slope.x_index]) * (bin - x_offset);
  if (curve_enabled)
    ret += curve.val_at(fit[curve.x_index]) * square(bin - x_offset);
  return ret;
}

double PolyBackground::eval_grad(double bin, std::vector<double>& gradients) const
{
  double ret = base.val();
  gradients[base.x_index] = base.grad();
  if (slope_enabled)
  {
    ret += slope.val() * (bin - x_offset);
    gradients[slope.x_index] = (bin - x_offset);
  }

  if (curve_enabled)
  {
    ret += curve.val() * square(bin - x_offset);
    gradients[curve.x_index] = square(bin - x_offset);
  }
  return ret;
}

double PolyBackground::eval_grad_at(double bin,
                                    const std::vector<double>& fit,
                                    std::vector<double>& gradients) const
{
  double ret = base.val_at(fit[base.x_index]);
  gradients[base.x_index] = base.grad_at(fit[base.x_index]);
  if (slope_enabled)
  {
    ret += slope.val_at(fit[slope.x_index]) * (bin - x_offset);
    gradients[slope.x_index] = (bin - x_offset);
  }

  if (curve_enabled)
  {
    ret += curve.val_at(fit[curve.x_index]) * square(bin - x_offset);
    gradients[curve.x_index] = square(bin - x_offset);
  }
  return ret;
}

double PolyBackground::eval_add(const std::vector<double>& bins, std::vector<double>& vals) const
{
  if (vals.size() != bins.size())
    vals.resize(bins.size(), 0.0);
  for (size_t i=0; i < bins.size(); ++i)
    vals[i] += eval(bins[i]);
}

std::vector<double> PolyBackground::eval(const std::vector<double>& bins) const
{
  std::vector<double> ret;
  eval_add(bins, ret);
  return ret;
}

std::string PolyBackground::to_string() const
{
  std::stringstream ss;
  ss << "x=bin-" << x_offset << "    ";
  ss << "base=" << base.to_string();
  if (slope_enabled)
    ss << "   base=" << slope.to_string();
  if (curve_enabled)
    ss << "   base=" << curve.to_string();
  return ss.str();
}

void to_json(nlohmann::json& j, const PolyBackground& s)
{
  j["bin_offset"] = s.x_offset;
  j["base"] = s.base;
  if (s.slope_enabled)
    j["slope"] = s.slope;
  if (s.curve_enabled)
    j["curve"] = s.curve;
}

void from_json(const nlohmann::json& j, PolyBackground& s)
{
  s.x_offset = j["bin_offset"];
  s.base = j["base"];
  if (j.count("slope"))
    s.slope = j["slope"];
  if (j.count("curve"))
    s.curve = j["curve"];
}

}
