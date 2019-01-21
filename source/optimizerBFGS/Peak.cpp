#include <optimizerBFGS/Peak.h>

#include <core/util/custom_logger.h>
#include <optimizerBFGS/more_math.h>

namespace Hypermet
{

int32_t Peak::step_type() const
{
  return FEP_status_;
}

double Peak::peak_position() const
{
  return position.val();
}

double Peak::peak_position_unc() const
{
  return position.uncert_value;
}

double Peak::peak_energy(const Calibration& cal) const
{
  return cal.channel_to_energy(position.val());
}

double Peak::peak_energy_unc(const Calibration& cal) const
{
  return cal.energy_slope() * position.uncert_value;
}

bool Peak::full_energy_peak() const
{
  if (FEP_status_ == 1)
    return true;
  else if (FEP_status_ == -1)
    return false;
  return false; // false? Ask Laci
}

void Peak::full_energy_peak(bool flag)
{
  if (flag)
    FEP_status_ = 1;
  else
    FEP_status_ = -1;
}

bool Peak::operator<(const Peak& other) const
{
  return position.val() < other.position.val();
}

}
