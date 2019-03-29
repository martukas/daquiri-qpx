#pragma once

#include <string>
#include <list>

#include <nlohmann/json.hpp>

class UncertainDouble
{
 public:
  UncertainDouble() = default;
  UncertainDouble(double val, double sigma);

  static UncertainDouble from_int(int64_t val, double sigma);
  static UncertainDouble from_uint(uint64_t val, double sigma);

  double value() const;
  double sigma() const;
  double error() const;
  double error_percent() const;

  void set_value(double val);
  void set_sigma(double sigma);

  // \todo make static function instead?
  bool is_finite() const;

  std::string to_string(bool ommit_tiny = true) const;
  std::string error_percent_fancy() const;

  std::string debug() const;

  UncertainDouble& operator*=(const double& other);
  UncertainDouble operator*(const double& other) const;
  UncertainDouble& operator/=(const double& other);
  UncertainDouble operator/(const double& other) const;

  UncertainDouble& operator*=(const UncertainDouble& other);
  UncertainDouble operator*(const UncertainDouble& other) const;
  UncertainDouble& operator/=(const UncertainDouble& other);
  UncertainDouble operator/(const UncertainDouble& other) const;

  UncertainDouble& operator+=(const UncertainDouble& other);
  UncertainDouble operator+(const UncertainDouble& other) const;
  UncertainDouble& operator-=(const UncertainDouble& other);
  UncertainDouble operator-(const UncertainDouble& other) const;

  bool almost(const UncertainDouble& other) const;
  bool operator==(const UncertainDouble& other) const;
  bool operator!=(const UncertainDouble& other) const;
  bool operator<(const UncertainDouble& other) const;
  bool operator>(const UncertainDouble& other) const;

  static UncertainDouble average(const std::list<UncertainDouble>& list);

 private:
  double value_{std::numeric_limits<double>::quiet_NaN()};
  double sigma_{std::numeric_limits<double>::quiet_NaN()};

  UncertainDouble& additive_uncert(const UncertainDouble& other);
  UncertainDouble& multipli_uncert(const UncertainDouble& other);
  int exponent() const;
};

void to_json(nlohmann::json& j, const UncertainDouble& s);
void from_json(const nlohmann::json& j, UncertainDouble& s);
