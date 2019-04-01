#pragma once

#include <core/fitting/data_model/data_model.h>
#include <nlohmann/json.hpp>

namespace DAQuiri
{

class CoefFunction : public FittableRegion
{
 public:
  std::vector<double> eval_vector(const std::vector<double>& x) const;
  double inverse(double y, double e = 0.1) const;

  //TO IMPLEMENT IN CHILDREN
  virtual bool valid() const = 0;
  virtual std::string type() const = 0;
  virtual CoefFunction* clone() const = 0;
  virtual bool is_equal(CoefFunction* other) const = 0;

  virtual double d_dx(double x) const = 0;

  virtual std::string to_UTF8(int precision, bool with_rsq) const = 0;
  virtual std::string to_markup(int precision, bool with_rsq) const = 0;
};

void to_json(nlohmann::json& j, const CoefFunction& s);
void from_json(const nlohmann::json& j, CoefFunction& s);

using CoefFunctionPtr = std::shared_ptr<CoefFunction>;

}
