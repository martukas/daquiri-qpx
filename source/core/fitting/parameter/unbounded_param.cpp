#include <core/fitting/parameter/unbounded_param.h>
#include <core/util/more_math.h>

#include <core/util/logger.h>

namespace DAQuiri
{

void UnboundedParam::val(double new_val)
{
  x(new_val);
}

double UnboundedParam::val_at(double at_x) const
{
  return at_x;
}

double UnboundedParam::grad_at(double at_x) const
{
  (void) at_x;
  return 1.0;
}

}
