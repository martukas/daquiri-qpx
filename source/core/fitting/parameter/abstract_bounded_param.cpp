#include <core/fitting/parameter/abstract_bounded_param.h>
#include <core/util/more_math.h>
#include <core/util/compare.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

double BoundedValue::max() const
{
  return max_;
}

void BoundedValue::max(double new_max)
{
  max_ = new_max;
  this->val(std::min(max_, this->val()));
}

double BoundedValue::min() const
{
  return min_;
}

void BoundedValue::min(double new_min)
{
  min_ = new_min;
  this->val(std::max(min_, this->val()));
}

void BoundedValue::bound(double v1, double v2)
{
  min(std::min(v1, v2));
  max(std::max(v1, v2));
}

void BoundedValue::set(double v1, double v2, double v3)
{
  min_ = ::min(v1, v2, v3);
  max_ = ::max(v1, v2, v3);
  val(mid(v1, v2, v3));
}


bool BoundedValue::at_extremum(double min_epsilon, double max_epsilon) const
{
  return ((val() - min()) < min_epsilon) || ((max() - val()) < max_epsilon);
}

std::string BoundedValue::to_string() const
{
  auto bounds_part = fmt::format("[{:<14}-{:>14}]", min_, max_);
  return fmt::format("{} {:>30}", AbstractValue::to_string(), bounds_part);
}

void to_json(nlohmann::json& j, const BoundedValue& s)
{
  j["x"] = s.x();
  j["to_fit"] = s.to_fit;
  j["uncert_value"] = s.uncert();
  j["min"] = s.min();
  j["max"] = s.max();
}

void from_json(const nlohmann::json& j, BoundedValue& s)
{
  s.x(j["x"]);
  s.to_fit = j["to_fit"];
  s.uncert(j["uncert_value"]);
  s.min(j["min"]);
  s.max(j["max"]);
}

}
