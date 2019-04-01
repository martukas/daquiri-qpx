#include <core/fitting/parameter/positive_param.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void PositiveValue::val(double new_val)
{
  x(std::sqrt(new_val));
}

double PositiveValue::val_at(double at_x) const
{
  return square(at_x);
}

double PositiveValue::grad_at(double at_x) const
{
  return 2.0 * at_x;
}

}
