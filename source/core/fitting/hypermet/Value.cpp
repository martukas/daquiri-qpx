#include <core/fitting/hypermet/Value.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void AbstractValue::update_index(int32_t& idx)
{
  if (idx < 0)
    throw std::runtime_error("Value cannot save negative variable index");

  if (to_fit)
    index_ = idx++;
  else
    index_ = -1;
}

int32_t AbstractValue::index() const
{
  return index_;
}

double AbstractValue::x() const
{
  return x_;
}

void AbstractValue::x(double new_x)
{
  x_ = new_x;
}

double AbstractValue::val() const
{
  return this->val_at(x_);
}

double AbstractValue::grad() const
{
  return this->grad_at(x_);
}

double AbstractValue::val_from(const Eigen::VectorXd& fit) const
{
  // \todo access without range checking once we have tests
  if (index_ >= 0)
    return this->val_at(fit(static_cast<size_t>(index_)));
  return val();
}

double AbstractValue::grad_from(const Eigen::VectorXd& fit) const
{
  // \todo access without range checking once we have tests
  if (index_ >= 0)
    return this->grad_at(fit(static_cast<size_t>(index_)));
  return grad();
}

double AbstractValue::uncert() const
{
  return val_uncert_;
}

void AbstractValue::put(Eigen::VectorXd& fit) const
{
  if (index_ >= 0)
    fit[index_] = x();
}

void AbstractValue::get(const Eigen::VectorXd& fit)
{
  if (index_ >= 0)
    x(fit[index_]);
}

void AbstractValue::get_uncert(const Eigen::VectorXd& diagonals, double chisq_norm)
{
  if (index_ >= 0)
    val_uncert_ = std::sqrt(std::abs(diagonals[index_] * this->grad() * chisq_norm));
}

std::string AbstractValue::to_string() const
{
  return fmt::format("{}\u00B1{}(x={},i={}{})",
      val(), val_uncert_, x_, index_,
      to_fit ? " fit" : "");
}

void to_json(nlohmann::json& j, const AbstractValue& s)
{
  j["x"] = s.x_;
  j["x_index"] = s.index_;
  j["to_fit"] = s.to_fit;
  j["uncert_value"] = s.val_uncert_;
}

void from_json(const nlohmann::json& j, AbstractValue& s)
{
  s.x_ = j["x"];
  s.index_ = j["x_index"];
  s.to_fit = j["to_fit"];
  s.val_uncert_ = j["uncert_value"];
}




double Value::max() const
{
  return max_;
}

void Value::max(double new_max)
{
  max_ = new_max;
  this->val(std::min(max_, this->val()));
}

double Value::min() const
{
  return min_;
}

void Value::min(double new_min)
{
  min_ = new_min;
  this->val(std::max(min_, this->val()));
}

void Value::bound(double v1, double v2)
{
  min(std::min(v1, v2));
  max(std::max(v1, v2));
}

void Value::val(double new_val)
{
  double t = (min_ + max_ - 2.0 * new_val) / (min_ - max_);
  if (std::abs(t) <= 1)
    x(std::asin((min_ + max_ - 2.0 * new_val) / (min_ - max_)));
  else if (signum(t) < 0)
    x(std::asin(-1));
  else
    x(std::asin(1));
}

double Value::val_at(double at_x) const
{
  return min_ + (max_ - min_) / 2.0 * (1.0 + std::sin(at_x));
}

double Value::grad_at(double at_x) const
{
  return (max_ - min_) * std::cos(at_x) / 2.0;
}

std::string Value::to_string() const
{
  return fmt::format("{} [{},{}]", AbstractValue::to_string(), min_, max_);
}

void to_json(nlohmann::json& j, const Value& s)
{
  j["x"] = s.x();
  j["x_index"] = s.index_;
  j["to_fit"] = s.to_fit;
  j["uncert_value"] = s.val_uncert_;
  j["min"] = s.min();
  j["max"] = s.max();
}

void from_json(const nlohmann::json& j, Value& s)
{
  s.x(j["x"]);
  s.index_ = j["x_index"];
  s.to_fit = j["to_fit"];
  s.val_uncert_ = j["uncert_value"];
  s.min(j["min"]);
  s.max(j["max"]);
}



void ValueGam::val(double new_val)
{
  x(std::sqrt(new_val));
}

double ValueGam::val_at(double at_x) const
{
  return square(at_x);
}

double ValueGam::grad_at(double at_x) const
{
  return 2.0 * at_x;
}



void ValueBkg::val(double new_val)
{
  x(new_val);
}

double ValueBkg::val_at(double at_x) const
{
  return at_x;
}

double ValueBkg::grad_at(double at_x) const
{
  (void) at_x;
  return 1.0;
}

}
