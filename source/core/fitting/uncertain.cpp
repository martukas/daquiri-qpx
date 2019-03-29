#include <core/fitting/uncertain.h>

#include <iostream>
#include <limits>
#include <cmath>
#include <core/util/lexical_extensions.h>
#include <core/util/UTF_extensions.h>

UncertainDouble::UncertainDouble(double val, double sigma)
    : value_(val)
      , sigma_(std::abs(sigma)) {}

UncertainDouble UncertainDouble::from_int(int64_t val, double sigma)
{
  return UncertainDouble(static_cast<double>(val), sigma);
}

UncertainDouble UncertainDouble::from_uint(uint64_t val, double sigma)
{
  return UncertainDouble(static_cast<double>(val), sigma);
}

void UncertainDouble::set_value(double val)
{
  value_ = val;
}

void UncertainDouble::set_sigma(double sigma)
{
  sigma_ = std::abs(sigma);
}

double UncertainDouble::value() const
{
  return value_;
}

double UncertainDouble::sigma() const
{
  return sigma_;
}

double UncertainDouble::error() const
{
  if (!is_finite())
    return std::numeric_limits<double>::quiet_NaN();
  else if (value_ != 0)
    return std::abs(sigma_ / value_);
  else
    return std::numeric_limits<double>::infinity();
}

double UncertainDouble::error_percent() const
{
  if (!is_finite())
    return std::numeric_limits<double>::quiet_NaN();
  else if (value_ != 0)
    return std::abs(100.0 * sigma_ / value_);
  else
    return std::numeric_limits<double>::infinity();
}

bool UncertainDouble::is_finite() const
{
  return (std::isfinite(value_));
}

std::string UncertainDouble::debug() const
{
  std::stringstream ss;
  ss << "["
     << std::to_string(value_)
     << "\u00B1"
     << std::to_string(sigma_);
  ss << "] ==> " << to_string(false);
  return ss.str();
}

int UncertainDouble::exponent() const
{
  int orderOfValue = order_of(value_);
  int orderOfUncert = order_of(sigma_);
  int targetOrder = std::max(orderOfValue, orderOfUncert);

  if ((targetOrder > 5) || (targetOrder < -3))
    return targetOrder;
  else
    return 0;
}

std::string UncertainDouble::to_string(bool ommit_tiny) const
{
  if (!std::isfinite(value_))
    return "?";

  int orderOfValue = order_of(value_);
  int orderOfUncert = order_of(sigma_);
  int exp = exponent();

  int decimals = 0;

  std::string result;
  if (std::isinf(sigma_))
    result = "~";

  if (ommit_tiny && std::isfinite(sigma_) && ((orderOfValue - orderOfUncert) < -2))
    result = " . ";
  else
    result += to_str_decimals(value_ / pow(10.0, exp), decimals);

  if (std::isfinite(sigma_) && (sigma_ != 0.0))
  {
    std::string uncertstr;
    if (ommit_tiny && ((orderOfValue - orderOfUncert) > 8))
      uncertstr = "-";
    else if (((orderOfValue - orderOfUncert) < -2)
        || ((orderOfValue - orderOfUncert) > 6))
      uncertstr = "HUGE"; //error_percent();
    else
    {
      double unc_shifted = sigma_ / pow(10.0, exp);
      if (!decimals)
        uncertstr = to_str_decimals(unc_shifted, 0);
      else if (unc_shifted < 1.0)
        uncertstr = to_str_decimals(unc_shifted / pow(10.0, -decimals), 0);
      else
        uncertstr = to_str_decimals(unc_shifted, orderOfUncert - exp + decimals);
    }

    if (!uncertstr.empty())
      result += "(" + uncertstr + ")";
  }

  std::string times_ten("\u00D710");
  if (exp)
    result += times_ten + UTF_superscript(exp);

  return result;
}

std::string UncertainDouble::error_percent_fancy() const
{
  if (!is_finite())
    return "?";
  if ((sigma_ == 0.0) || !std::isfinite(sigma_))
    return "-";
  if (value_ == 0)
    return "inf";

  double error = std::abs(sigma_ / value_ * 100.0);

  UncertainDouble p(error, 0);
//  DBG << "perror for " << debug() << " is " << p.debug();
  if (p.exponent() != 0)
    return "(" + p.to_string(false) + ")%";
  else
    return p.to_string(false) + "%";
}

UncertainDouble& UncertainDouble::operator*=(const UncertainDouble& other)
{
  set_value(value_ * other.value_);
  if (is_finite() && other.is_finite())
    set_sigma(sqrt(value_ * value_ * other.sigma_ * other.sigma_
                       + other.value_ * other.value_ * sigma_ * sigma_));
  else
    set_sigma(std::numeric_limits<double>::quiet_NaN());
  return *this;
}

UncertainDouble& UncertainDouble::operator/=(const UncertainDouble& other)
{
  set_value(value_ / other.value_);
  if (is_finite() && other.is_finite())
    set_sigma(sqrt(value_ * value_ * other.sigma_ * other.sigma_
                       + other.value_ * other.value_ * sigma_ * sigma_));
  else
    set_sigma(std::numeric_limits<double>::quiet_NaN());
  return *this;
}

UncertainDouble& UncertainDouble::operator*=(const double& other)
{
  set_value(value_ * other);
  set_sigma(std::abs(sigma_ * other));
  return *this;
}

UncertainDouble& UncertainDouble::operator/=(const double& other)
{
  set_value(value_ / other);
  set_sigma(std::abs(sigma_ / other));
  return *this;
}

UncertainDouble& UncertainDouble::additive_uncert(const UncertainDouble& other)
{
  if (is_finite() && other.is_finite())
    set_sigma(sqrt(pow(sigma_, 2) + pow(other.sigma_, 2)));
  else
    set_sigma(std::numeric_limits<double>::quiet_NaN());
  return *this;
}

UncertainDouble& UncertainDouble::operator+=(const UncertainDouble& other)
{
  set_value(value_ + other.value_);
  additive_uncert(other);
  return *this;
}

UncertainDouble& UncertainDouble::operator-=(const UncertainDouble& other)
{
  set_value(value_ - other.value_);
  additive_uncert(other);
  return *this;
}

UncertainDouble UncertainDouble::operator+(const UncertainDouble& other) const
{
  UncertainDouble result(*this);
  result += other;
  return result;
}

UncertainDouble UncertainDouble::operator-(const UncertainDouble& other) const
{
  UncertainDouble result(*this);
  result -= other;
  return result;
}

UncertainDouble UncertainDouble::operator*(const UncertainDouble& other) const
{
  UncertainDouble result(*this);
  result *= other;
  return result;
}

UncertainDouble UncertainDouble::operator*(const double& other) const
{
  UncertainDouble result(*this);
  result *= other;
  return result;
}

UncertainDouble UncertainDouble::operator/(const UncertainDouble& other) const
{
  UncertainDouble result(*this);
  result /= other;
  return result;
}

UncertainDouble UncertainDouble::operator/(const double& other) const
{
  UncertainDouble result(*this);
  result /= other;
  return result;
}

bool UncertainDouble::operator==(const UncertainDouble& other) const
{
  return (value() == other.value()) && (sigma() == other.sigma());
}

bool UncertainDouble::operator!=(const UncertainDouble& other) const
{
  return (value() != other.value()) || (sigma() != other.sigma());
}

bool UncertainDouble::operator<(const UncertainDouble& other) const
{
  return (value() < other.value()) || (sigma() < other.sigma());
}

bool UncertainDouble::operator>(const UncertainDouble& other) const
{
  return (value() > other.value()) || (sigma() > other.sigma());
}

bool UncertainDouble::almost(const UncertainDouble& other) const
{
  if (value_ == other.value_)
    return true;
  double delta = std::abs(value_ - other.value_);
  if (std::isfinite(sigma_) && (delta <= sigma_))
    return true;
  if (std::isfinite(other.sigma_) && (delta <= other.sigma_))
    return true;
  return false;
}

UncertainDouble UncertainDouble::average(const std::list<UncertainDouble>& list)
{
  if (list.empty())
    return UncertainDouble();

  double sum = 0;
  for (auto& l : list)
    sum += l.value_;
  double avg = (sum) / list.size();
  double min = avg;
  double max = avg;
  for (auto& l : list)
  {
    if (std::isfinite(l.sigma_))
    {
      min = std::min(min, l.value_ - l.sigma_);
      max = std::max(max, l.value_ + l.sigma_);
    }
  }

  UncertainDouble ret((max + min) * 0.5, (max - min) * 0.5);
  return ret;
}

void to_json(nlohmann::json& j, const UncertainDouble& s)
{
  if (!std::isnan(s.value()))
    j["value"] = s.value();
  if (!std::isnan(s.sigma()))
    j["sigma"] = s.sigma();
}

void from_json(const nlohmann::json& j, UncertainDouble& s)
{
  double val = s.value();
  if (j.count("value") && j["value"].is_number_float())
    val = j["value"];
  double sigma = s.sigma();
  if (j.count("sigma") && j["sigma"].is_number_float())
    sigma = j["sigma"];
  s = UncertainDouble(val, sigma);
}
