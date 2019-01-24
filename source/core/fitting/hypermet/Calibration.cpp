#include <core/fitting/hypermet/Calibration.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace Hypermet
{

Calibration::Calibration()
{
  //EnergyCal[0].Channel = 0;
  //EnergyCal[0].val = 0;
  //EnergyCal[1].Channel = 1;
  //EnergyCal[1].val = 1;
}

uint8_t Calibration::order() const
{
  return order_;
}
void Calibration::order(uint8_t new_order)
{
  if (new_order <= 2)
  {
    order_ = new_order;
    energy_cal.resize(new_order);
    width_cal.resize(new_order);
  }
}

double Calibration::energy_const() const
{
  return -1 / (energy_cal[1].channel - energy_cal[0].channel) *
      energy_cal[0].channel * energy_cal[1].value + 1 /
      (energy_cal[1].channel - energy_cal[0].channel) *
      energy_cal[0].channel * energy_cal[0].value + energy_cal[0].value;
}

double Calibration::energy_slope() const
{
  return (1 / (energy_cal[1].channel - energy_cal[0].channel)
      * energy_cal[1].value - 1 /
      (energy_cal[1].channel - energy_cal[0].channel) * energy_cal[0].value);
}

double Calibration::channel_to_energy(double chan) const
{
  chan += nonlinearity.val(chan);
  double c0 = -1 / (energy_cal[1].channel - energy_cal[0].channel) *
      energy_cal[0].channel * energy_cal[1].value + 1 /
      (energy_cal[1].channel - energy_cal[0].channel) *
      energy_cal[0].channel * energy_cal[0].value + energy_cal[0].value;
  //(1 / (Ch2 - Ch1) * e2 - 1 / (Ch2 - Ch1) * e1)
  double c1 = (1 / (energy_cal[1].channel - energy_cal[0].channel) *
      (energy_cal[1].value - energy_cal[0].value));
  return c0 + c1 * chan;
}

double Calibration::energy_to_channel(double energy) const
{
  double c0 = -1 / (energy_cal[1].channel - energy_cal[0].channel) *
      energy_cal[0].channel * energy_cal[1].value + 1 /
      (energy_cal[1].channel - energy_cal[0].channel) * energy_cal[0].channel *
      energy_cal[0].value + energy_cal[0].value;
  double
      c1 = (1 / (energy_cal[1].channel - energy_cal[0].channel) *
      energy_cal[1].value - 1 / (energy_cal[1].channel - energy_cal[0].channel) *
      energy_cal[0].value);
  energy -= c0;
  energy /= c1;
  energy -= nonlinearity.val(energy);
  return energy;
}

double Calibration::width(double channel) const
{
  //fwhm_b = (-fwhm2 ^ 2 + fwhm1 ^ 2) / (-Energy(chfwhm2) + Energy(chfwhm1))
  //fwhm_a = -(Energy(chfwhm2) * fwhm1 ^ 2 - fwhm2 ^ 2 * Energy(chfwhm1)) / (-Energy(chfwhm2) + Energy(chfwhm1))
}

}
