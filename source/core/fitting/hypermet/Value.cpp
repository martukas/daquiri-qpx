#include <core/fitting/hypermet/Value.h>
#include <core/util/more_math.h>

//#include <mpreal.h>

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
    reset_index();
}

void AbstractValue::reset_index()
{
  index_ = InvalidIndex;
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
  if (index_ > InvalidIndex)
    return this->val_at(fit(static_cast<size_t>(index_)));
  return val();
}

double AbstractValue::grad_from(const Eigen::VectorXd& fit) const
{
  // \todo access without range checking once we have tests
  if (index_ > InvalidIndex)
    return this->grad_at(fit(static_cast<size_t>(index_)));
  return grad();
}

double AbstractValue::uncert() const
{
  return val_uncert_;
}

void AbstractValue::put(Eigen::VectorXd& fit) const
{
  if (index_ > InvalidIndex)
    fit[index_] = x();
}

void AbstractValue::get(const Eigen::VectorXd& fit)
{
  if (index_ > InvalidIndex)
    x(fit[index_]);
}

void AbstractValue::get_uncert(const Eigen::VectorXd& diagonals, double chisq_norm)
{
  if (index_ > InvalidIndex)
    val_uncert_ = std::sqrt(std::abs(diagonals[index_] * this->grad() * chisq_norm));
}

std::string AbstractValue::to_string() const
{
  auto val_part = fmt::format("{}\u00B1{}", val(), val_uncert_);
  auto x_part = fmt::format("(x={})", x_);
  auto i_part = fmt::format("{}[{}]",
      (to_fit ? "F" : ""),
      ((index_ > InvalidIndex) ? std::to_string(index_) : "-"));
  return fmt::format("{:>20} {:<17} {:>8}",
                     val_part, x_part, i_part);
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
//  mpfr::mpreal mx = at_x;
//  mpfr::mpreal one = 1.0;
//  mpfr::mpreal two = 2.0;
//  mpfr::mpreal mmax = max_;
//  mpfr::mpreal mmin = min_;
//  auto ret = (one + mpfr::sin(mx)) * (mmax - mmin) / two + mmin;
//  return ret.toDouble();
  return (1.0 + std::sin(at_x)) * (max_ - min_) / 2.0 + min_;
}

double Value::grad_at(double at_x) const
{
//  mpfr::mpreal mx = at_x;
//  mpfr::mpreal two = 2.0;
//  mpfr::mpreal mmax = max_;
//  mpfr::mpreal mmin = min_;
//  auto ret = mpfr::cos(mx) * (mmax - mmin) / two;
//  return ret.toDouble();
  return std::cos(at_x) * (max_ - min_) / 2.0;
}

std::string Value::to_string() const
{
  auto bounds_part = fmt::format("[{:<14}-{:>14}]", min_, max_);
  return fmt::format("{} {:>30}", AbstractValue::to_string(), bounds_part);
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



void ValuePositive::val(double new_val)
{
  x(std::sqrt(new_val));
}

double ValuePositive::val_at(double at_x) const
{
  return square(at_x);
}

double ValuePositive::grad_at(double at_x) const
{
  return 2.0 * at_x;
}



void ValueSimple::val(double new_val)
{
  x(new_val);
}

double ValueSimple::val_at(double at_x) const
{
  return at_x;
}

double ValueSimple::grad_at(double at_x) const
{
  (void) at_x;
  return 1.0;
}

}
