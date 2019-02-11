#include <core/fitting/hypermet/PolyBackground.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

PolyBackground::PolyBackground()
{
  base.bound(0, 5000000);
  base.val(0);
  slope.bound(-100, 100);
  curve.bound(-10, 10);
}

void PolyBackground::update_indices(int32_t& i)
{
  base.update_index(i);
  slope.update_index(i);
  curve.update_index(i);
}

void PolyBackground::put(Eigen::VectorXd& fit) const
{
  base.put(fit);
  slope.put(fit);
  curve.put(fit);
}

void PolyBackground::get(const Eigen::VectorXd& fit)
{
  base.get(fit);
  slope.get(fit);
  curve.get(fit);
}

void PolyBackground::get_uncerts(const Eigen::VectorXd& diagonals, double chisq_norm)
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

double PolyBackground::eval_at(double bin, const Eigen::VectorXd& fit) const
{
  double ret = base.val_from(fit);
  if (slope_enabled)
    ret += slope.val_from(fit) * (bin - x_offset);
  if (curve_enabled)
    ret += curve.val_from(fit) * square(bin - x_offset);
  return ret;
}

double PolyBackground::eval_grad(double bin, Eigen::VectorXd& gradients) const
{
  double ret = base.val();
  gradients[base.index()] = base.grad();
  if (slope_enabled)
  {
    ret += slope.val() * (bin - x_offset);
    gradients[slope.index()] = (bin - x_offset);
  }

  if (curve_enabled)
  {
    ret += curve.val() * square(bin - x_offset);
    gradients[curve.index()] = square(bin - x_offset);
  }
  return ret;
}

double PolyBackground::eval_grad_at(double bin,
                                    const Eigen::VectorXd& fit,
                                    Eigen::VectorXd& gradients) const
{
  double ret = base.val_from(fit);
  gradients[base.index()] = base.grad_from(fit);
  if (slope_enabled)
  {
    ret += slope.val_from(fit) * (bin - x_offset);
    gradients[slope.index()] = (bin - x_offset);
  }

  if (curve_enabled)
  {
    ret += curve.val_from(fit) * square(bin - x_offset);
    gradients[curve.index()] = square(bin - x_offset);
  }
  return ret;
}

void PolyBackground::eval_add(const std::vector<double>& bins, std::vector<double>& vals) const
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

std::string PolyBackground::to_string(std::string prepend) const
{
  std::stringstream ss;
  ss << prepend << "x=bin-" << x_offset << "\n";
  ss << prepend << "base=" << base.to_string() << "\n";
  if (slope_enabled)
    ss << prepend << "slope=" << slope.to_string() << "\n";
  if (curve_enabled)
    ss << prepend << "curve=" << curve.to_string() << "\n";
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
