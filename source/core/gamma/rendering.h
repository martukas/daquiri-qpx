#pragma once

#include <core/gamma/hypermet/gamma_region.h>
#include <core/gamma/fit_settings.h>

namespace DAQuiri
{

struct PeakRendering
{
  std::vector<double>
      peak,
      full_fit;

  void clear();
  void render(const Peak&);
};

struct RegionRendering
{
  uint8_t subdivisions{10};

  std::vector<double>
      channel,
      energy,
      background,
      back_steps,
      full_fit,
      sum4_background;

  std::map<double, PeakRendering> peaks;

  void reserve(size_t count);
  void clear();

  void render(const Region& r, const Calibration& energy_calib);
};

}
