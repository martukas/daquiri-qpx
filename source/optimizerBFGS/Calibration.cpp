#include <optimizerBFGS/Calibration.h>
#include <optimizerBFGS/more_math.h>

#include <core/util/custom_logger.h>

namespace Hypermet
{

CCalibration::CCalibration()
{
  //EnergyCal[0].Channel = 0;
  //EnergyCal[0].Value = 0;
  //EnergyCal[1].Channel = 1;
  //EnergyCal[1].Value = 1;
}

uint8_t CCalibration::CalOrder() const
{
  return _CalOrder;
}
void CCalibration::CalOrder(uint8_t Value)
{
  if (Value <= 2)
  {
    _CalOrder = Value;
    EnergyCal.resize(Value);
    WidthCal.resize(Value);
  }
}

double CCalibration::EnergyConst() const
{
  double c0 = -1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[1].Value + 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[0].Value + EnergyCal[0].Value;
  return c0;
}

double CCalibration::EnergySlope() const
{
  double c1 = (1 / (EnergyCal[1].Channel - EnergyCal[0].Channel)
      * EnergyCal[1].Value - 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) * EnergyCal[0].Value);
  return c1;
}

double CCalibration::ChannelToEnergy(double Channel) const
{
  Channel += Nonlinearity.Value(Channel);
  double c0 = -1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[1].Value + 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[0].Value + EnergyCal[0].Value;
  //(1 / (Ch2 - Ch1) * e2 - 1 / (Ch2 - Ch1) * e1)
  double c1 = (1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) * (EnergyCal[1].Value - EnergyCal[0].Value));
  return c0 + c1 * Channel;
}

double CCalibration::EnergyToChannel(double Energy) const
{
  double c0 = -1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[1].Value + 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) * EnergyCal[0].Channel *
      EnergyCal[0].Value + EnergyCal[0].Value;
  double
      c1 = (1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[1].Value - 1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Value);
  Energy -= c0;
  Energy /= c1;
  Energy -= Nonlinearity.Value(Energy);
  return Energy;
}

double CCalibration::Width(double Channel) const
{
  //fwhm_b = (-fwhm2 ^ 2 + fwhm1 ^ 2) / (-Energy(chfwhm2) + Energy(chfwhm1))
  //fwhm_a = -(Energy(chfwhm2) * fwhm1 ^ 2 - fwhm2 ^ 2 * Energy(chfwhm1)) / (-Energy(chfwhm2) + Energy(chfwhm1))
}

}
