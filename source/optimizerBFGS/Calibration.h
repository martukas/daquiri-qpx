#pragma once

#include <optimizerBFGS/EfficiencyCal.h>
#include <optimizerBFGS/NonlinearityCal.h>

namespace Hypermet
{

class CCalibration
{
 public:
  struct CalPoint
  {
    float Channel;
    float UncChannel;
    float Value;
    float UncValue;
  };

  std::vector<CalPoint> EnergyCal{2};
  std::vector<CalPoint> WidthCal{2};
  NonlinearityCal Nonlinearity;
  EfficiencyCal Efficiency;

 public:
  CCalibration();

  uint8_t CalOrder() const;
  void CalOrder(uint8_t Value);
  double EnergyConst() const;
  double EnergySlope() const;
  double ChannelToEnergy(double Channel) const;
  double EnergyToChannel(double Energy) const;
  double Width(double Channel) const;

 private:
  uint8_t _CalOrder;
};

}
