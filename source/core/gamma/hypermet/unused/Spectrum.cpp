#include <core/fitting/hypermet/Spectrum.h>
#include <core/util/more_math.h>

#include <core/util/logger.h>

namespace DAQuiri
{

// \todo these go somewhere else

double dead_time(double real_time, double live_time)
{
  if (real_time > 0.0)
    return (real_time - live_time) / real_time * 100.0;
  return 0.0;
}

int value_quality(UncertainDouble ud, double error_threshold)
{
  if (ud.error() > error_threshold)
    return 3;
  else if (!std::isfinite(ud.sigma()) || !ud.sigma())
    return 2;
  return 1;
}

int peak_good(const Peak& h, const SUM4& s)
{
  return ((s.quality() == 1)
      && (value_quality(h.peak_position()) == 1)
      && (value_quality(h.fwhm()) == 1));
}


}
