#include <core/fitting/hypermet/Spectrum.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

// \todo these go somewhere else

double dead_time(double real_time, double live_time)
{
  if (real_time > 0.0)
    return (real_time - live_time) / real_time * 100.0;
  return 0.0;
}

double rate(double live_time, double sum_counts)
{
  if (live_time > 0.0)
    return sum_counts / live_time * 100;
  return 0.0;
}

}
