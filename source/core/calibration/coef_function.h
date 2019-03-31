#pragma once

#include <core/calibration/parameter.h>

namespace DAQuiri
{

class CoefFunction
{
 public:
  CoefFunction() = default;
  virtual ~CoefFunction() = default;

//  void set_coeff(int degree, const Parameter& p);
//  std::map<int, Parameter> coeffs() const;

  std::vector<double> eval(const std::vector<double>& x) const;
  double inverse(double y, double e = 0.1) const;

  //TO IMPLEMENT IN CHILDREN
  virtual bool valid() const = 0;
  virtual std::string type() const = 0;
  virtual CoefFunction* clone() const = 0;
  virtual bool is_equal(CoefFunction* other) const = 0;

  virtual double operator() (double x) const = 0;
  virtual double derivative(double x) const = 0;

  virtual std::string debug() const = 0;
  virtual std::string to_UTF8(int precision, bool with_rsq) const = 0;
  virtual std::string to_markup(int precision, bool with_rsq) const = 0;
};

void to_json(nlohmann::json& j, const CoefFunction& s);
void from_json(const nlohmann::json& j, CoefFunction& s);

using CoefFunctionPtr = std::shared_ptr<CoefFunction>;

}
