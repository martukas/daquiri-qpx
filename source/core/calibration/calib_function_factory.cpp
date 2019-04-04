#include <core/calibration/calib_function_factory.h>

#include <core/util/logger.h>

namespace DAQuiri
{

void CalibFunctionFactory::register_type(std::function<CalibFunction*(void)> constructor)
{
  auto name = CalibFunctionPtr(constructor())->type();
  if (name.empty())
    WARN("<CalibFunctionFactory> failed to register nameless type");
  else if (constructors_.count(name))
    WARN("<CalibFunctionFactory> type '{}' already registered", name);
  else
  {
    constructors_[name] = constructor;
    DBG("<CalibFunctionFactory> registered '{}'", name);
  }
}

CalibFunctionPtr CalibFunctionFactory::create_type(std::string type) const
{
  auto it = constructors_.find(type);
  if (it != constructors_.end())
    return CalibFunctionPtr(it->second());
  else
    return CalibFunctionPtr();
}

CalibFunctionPtr CalibFunctionFactory::create_copy(CalibFunctionPtr other) const
{
  if (other)
    return CalibFunctionPtr(other->clone());
  return CalibFunctionPtr();
}

CalibFunctionPtr CalibFunctionFactory::create_from_json(const nlohmann::json& j) const
{
  std::string type = j["type"];
  auto ret = create_type(type);
  ret->from_json(j);
  return ret;
}

std::vector<std::string> CalibFunctionFactory::types() const
{
  std::vector<std::string> all_types;
  for (auto& q : constructors_)
    all_types.push_back(q.first);
  return all_types;
}

void CalibFunctionFactory::clear()
{
  constructors_.clear();
}


}
