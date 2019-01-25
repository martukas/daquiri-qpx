#include <core/fitting/hypermet/Value.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

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

void AbstractValue::put(std::vector<double>& fit) const
{
  if (x_index != -1)
    fit[x_index] = x();
}

void AbstractValue::get(const std::vector<double>& fit)
{
  if (x_index != -1)
    x(fit[x_index]);
}

void AbstractValue::get_uncert(const std::vector<double>& diagonals, double chisq_norm)
{
  if (x_index != -1)
    uncert_value = std::sqrt(std::abs(diagonals[x_index] * this->grad() * chisq_norm));
}

std::string AbstractValue::to_string() const
{
  return fmt::format("{}+-{}(x={},i={}{})",
      val(), uncert_value, x_, x_index,
      to_fit ? " fit" : "");
}

void to_json(nlohmann::json& j, const AbstractValue& s)
{
  j["x"] = s.x_;
  j["x_index"] = s.x_index;
  j["to_fit"] = s.to_fit;
  j["uncert_value"] = s.uncert_value;
}

void from_json(const nlohmann::json& j, AbstractValue& s)
{
  s.x_ = j["x"];
  s.x_index = j["x_index"];
  s.to_fit = j["to_fit"];
  s.uncert_value = j["uncert_value"];
}




double Value::max() const
{
  return max_;
}

void Value::max(double new_max)
{
  //Dim val As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(x_))
  //Value = std::min(NewMax, val)
  max_ = new_max;
  //x_ = std::Asin((_Min + _Max - 2 * val) / (-1 * _Max + _Min))
}

double Value::min() const
{
  return min_;
}

void Value::min(double new_min)
{
  //Dim val As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(x_))
  min_ = new_min;
  //temp = std::max(temp, val)
  //x_ = std::Asin((_Min + _Max - 2 * temp) / (-1 * _Max + _Min))
}

void Value::bound(double v1, double v2)
{
  min(std::min(v1, v2));
  max(std::max(v1, v2));
}

//void SetRange(double NewMin, double NewMax) {
//    if(Not NewMin < NewMax) { return; }
//    double val = _Min + (_Max - _Min) / 2 * (1 + std::sin(x_));
//    Value = std::min(NewMin, val);
//    Value = std::max(NewMax, val);
//    _Min = NewMin;
//    _Max = NewMax;
//    x_ = std::asin((_Min + _Max - 2 * val) / (-1 * _Max + _Min));
//}

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
  j["x_index"] = s.x_index;
  j["to_fit"] = s.to_fit;
  j["uncert_value"] = s.uncert_value;
  j["min"] = s.min();
  j["max"] = s.max();
}

void from_json(const nlohmann::json& j, Value& s)
{
  s.x(j["x"]);
  s.x_index = j["x_index"];
  s.to_fit = j["to_fit"];
  s.uncert_value = j["uncert_value"];
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

std::string side_to_string(const Side& s)
{
  if (s == Side::right)
    return "right";
  else
    return "left";
}

Side side_from_string(const std::string& s)
{
  if (s == "right")
    return Side::right;
  else
    return Side::left;
}


}
