#include <core/fitting/hypermet/Spectrum.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace Hypermet
{

double CSpectrum::weight_true(size_t i) const
{
  double k0 = channels[i];
  return std::sqrt(k0);
}

double CSpectrum::weight_phillips_marlow(size_t i) const
{
  double k0 = channels[i];

  if (k0 >= 25)
    return std::sqrt(k0);
  else
  {
    double k1 = 1;
    if ((i > 0) && (i < channels.size()))
      k1 = channels[i - 1] + channels[i] + channels[i + 1] / 3.0;
    return std::max(std::sqrt(k1), 1.0);
  }
}

double CSpectrum::weight_revay_student(size_t i) const
{
  double k0 = channels[i] + 1;
  return std::sqrt(k0);
}


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

size_t mystery_function(double val)
{
  bool Ready = false;
  if ((val >= 0) && (val < pow(2, 36)))
    Ready = true;

  while (!Ready)
  {
    size_t exponent = std::log(std::abs(val)) / std::log(2);
    val = val - signum(val) * pow(2, exponent);
    if ((val >= 0) && (val < pow(2, 36)))
      Ready = true;
  }
  return val;
}

}
