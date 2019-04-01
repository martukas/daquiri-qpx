#include <core/fitting/parameter/atan_bounded_param.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

double Value2::val_at(double at_x) const
{
  return (M_PI_2 + std::atan(slope_ * at_x)) * (max() - min()) * M_1_PI + min();
}

double Value2::grad_at(double at_x) const
{
  return (1.0 / (1.0 + square(slope_ * at_x))) * slope_ * (max() - min()) * M_1_PI;
}

void Value2::val(double new_val)
{
  if (new_val >= max())
    x(std::numeric_limits<double>::max());
  else if (new_val <= min())
    x(-std::numeric_limits<double>::max());
  else
    x(std::tan(M_PI * (new_val - min()) / (max() - min()) - M_PI_2) / slope_);
}


}
