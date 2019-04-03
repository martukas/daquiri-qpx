#pragma once

#include <core/calibration/calib_function.h>
#include <core/util/unique_mangle.h>

namespace DAQuiri
{

class CalibFunctionFactory
{
 public:
  static CalibFunctionFactory& singleton()
  {
    static CalibFunctionFactory singleton_instance;
    return singleton_instance;
  }

  CalibFunctionPtr create_type(std::string type) const;
  CalibFunctionPtr create_copy(CalibFunctionPtr other) const;
  CalibFunctionPtr create_from_json(const nlohmann::json&) const;

  void register_type(std::function<CalibFunction*(void)> constructor);
  std::vector<std::string> types() const;

  void clear();

 private:
  std::map<std::string, std::function<CalibFunction*(void)>> constructors_;

  //singleton assurance
  CalibFunctionFactory() {}
  CalibFunctionFactory(CalibFunctionFactory const&);
  void operator=(CalibFunctionFactory const&);
};

template<class T>
class CalibFunctionRegistrar
{
 public:
  CalibFunctionRegistrar()
  {
    CalibFunctionFactory::singleton().register_type([](void) -> CalibFunction* { return new T(); });
  }
};

#define DAQUIRI_REGISTER_COEF_FUNCTION(T) static DAQuiri::CalibFunctionRegistrar< T >\
  UNIQUE_MANGLE(MangledDAQuiriCoefFuncReg) ;


}
