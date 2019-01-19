#include <optimizerBFGS/Spectrum.h>
#include <optimizerBFGS/more_math.h>

#include <core/util/custom_logger.h>

namespace Hypermet
{

double CSpectrum::Weight(size_t i) const
{
  double k0 = Channel[i];

  if (k0 >= 25)
    return std::sqrt(k0);
  else
  {
    double k1 = 1;
    if ((i > 0) && (i < Channel.size()))
      k1 = Channel[i - 1] + Channel[i] + Channel[i + 1] / 3.0;
    return std::max(std::sqrt(k1), 1.0);
  }
}

double CSpectrum::DeadTime(double TrueTime, double LiveTime)
{
  if (TrueTime > 0.0)
    return (TrueTime - LiveTime) / TrueTime * 100.0;
  return 0.0;
}

double CSpectrum::Rate(double LiveTime, double SumCounts)
{
  if (LiveTime > 0.0)
    return SumCounts / LiveTime * 100;
  return 0.0;
}

size_t CSpectrum::mystery_function(double val)
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
