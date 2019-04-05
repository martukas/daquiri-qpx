#pragma once

#include <core/fitting/data_model/data_model.h>
#include <nlohmann/json.hpp>

namespace DAQuiri
{

class CalibFunction : public DataModel
{
 public:
  //This is for factory and serialization
  virtual std::string type() const = 0;
  virtual CalibFunction* clone() const = 0;
  virtual void from_json(const nlohmann::json&) {}
  virtual nlohmann::json to_json() const
  {
    nlohmann::json j;
    j["type"] = this->type();
    return j;
  };

  //This is for calibration comparison and validity
  virtual bool valid() const = 0;
  virtual bool is_equal(CalibFunction* other) const = 0;

  // This is for inverting calibration
  virtual double d_dx(double x) const = 0;

  // This is for pretty output
  virtual std::string to_UTF8(int precision) const = 0;
  virtual std::string to_markup(int precision) const = 0;
};

using CalibFunctionPtr = std::shared_ptr<CalibFunction>;

}
