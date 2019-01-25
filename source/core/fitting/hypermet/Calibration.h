#pragma once

#include <core/fitting/hypermet/EfficiencyCal.h>
#include <core/fitting/hypermet/NonlinearityCal.h>

namespace DAQuiri
{

class HCalibration
{
 public:
  struct CalPoint
  {
    float channel;
    float channel_unc;
    float value;
    float value_unc;
  };

  std::vector<CalPoint> energy_cal{2};
  std::vector<CalPoint> width_cal{2};
  NonlinearityCal nonlinearity;
  EfficiencyCal efficiency;

 public:
  HCalibration();

  uint8_t order() const;
  void order(uint8_t new_order);
  double energy_const() const;
  double energy_slope() const;
  double channel_to_energy(double chan) const;
  double energy_to_channel(double energy) const;
  double width(double channel) const;

 private:
  uint8_t order_;
};

}
