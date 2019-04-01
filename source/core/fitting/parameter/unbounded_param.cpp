#include <core/fitting/parameter/unbounded_param.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void UnboundedValue::val(double new_val)
{
  x(new_val);
}

double UnboundedValue::val_at(double at_x) const
{
  return at_x;
}

double UnboundedValue::grad_at(double at_x) const
{
  (void) at_x;
  return 1.0;
}

}
