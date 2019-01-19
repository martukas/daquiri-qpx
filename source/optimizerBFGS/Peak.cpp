#include <optimizerBFGS/Peak.h>

#include <core/util/custom_logger.h>
#include <optimizerBFGS/more_math.h>

namespace Hypermet
{

double CValueDefault::Max() const
{
  return _Max;
}
void CValueDefault::Max(double NewMax)
{
  //Dim Value As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(_X))
  //Value = std::min(NewMax, Value)
  _Max = NewMax;
  //_X = std::Asin((_Min + _Max - 2 * Value) / (-1 * _Max + _Min))
}

double CValueDefault::Min()
{
  return _Min;
}
void CValueDefault::Min(double NewMin)
{
  //Dim Value As Double = _Min + (_Max - _Min) / 2 * (1 + std::Sin(_X))
  _Min = NewMin;
  //temp = std::max(temp, Value)
  //_X = std::Asin((_Min + _Max - 2 * temp) / (-1 * _Max + _Min))
}

double CValueDefault::X()
{
  return _X;
}
void CValueDefault::X(double Value)
{
  _X = Value;
}

double CValueDefault::Value() const
{
  return _Min + (_Max - _Min) / 2 * (1 + std::sin(_X));
}
void CValueDefault::Value(double val)
{
  double t = (_Min + _Max - 2 * val) / (_Min - _Max);
  if (std::abs(t) <= 1)
    _X = std::asin((_Min + _Max - 2 * val) / (_Min - _Max));
  else if (signum(t) < 0)
    _X = std::asin(-1);
  else
    _X = std::asin(1);
}

double CValueDefault::ValueAt(double atX) const
{
  try
  {
    double ret = _Min + (_Max - _Min) / 2.0 * (1.0 + std::sin(atX));
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("ValueAt failed"));
  }
}

double CValueDefault::GradAt(double atX) const
{
  try
  {
    double ret = (_Max - _Min) * std::cos(atX) / 2.0;
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("GradAt failed"));
  }
}

double CValueGam::X()
{
  return _X;
}
void CValueGam::X(double Value)
{
  _X = Value;
}

double CValueGam::Value() const
{
  try
  {
    double ret = _X * _X;
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;

  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("Value failed"));
  }
}

void CValueGam::Value(double val)
{
  try
  {
    _X = std::sqrt(val);
    if (!std::isfinite(_X))
      throw std::runtime_error("Quiet NaN overflow occurred");
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("Value failed"));
  }
}

double CValueGam::ValueAt(double atX) const
{
  try
  {
    double ret = square(atX);
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("ValueAt failed"));
  }
}

double CValueGam::GradAt(double atX) const
{
  try
  {
    double ret = 2.0 * atX;
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("GradAt failed"));
  }
}

int32_t CPeak::StepType() const
{
  return FEPStatus;
}

double CPeak::PeakPosition() const
{
  return POS.Value();
}

double CPeak::UncPeakPosition() const
{
  return POS.UncValue;
}

double CPeak::PeakEnergy(const CCalibration& cal) const
{
  return cal.ChannelToEnergy(POS.Value());
}

double CPeak::UncPeakEnergy(const CCalibration& cal) const
{
  return cal.EnergySlope() * POS.UncValue;
}

bool CPeak::IsFullEnergyPeak() const
{
  if (FEPStatus == 1)
    return true;
  else if (FEPStatus == -1)
    return false;
  return false; // false? Ask Laci
}
void CPeak::IsFullEnergyPeak(bool Value)
{
  if (Value)
    FEPStatus = 1;
  else
    FEPStatus = -1;
}

bool CPeak::operator<(const CPeak& other) const
{
  return POS.Value() < other.POS.Value();
}

double CValueBkgDefault::X() const
{
  return _X;
}

void CValueBkgDefault::X(double Value)
{
  _X = Value;
}

double CValueBkgDefault::Value() const
{
  try
  {
    double ret = square(_X);
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("Value failed"));
  }
}

void CValueBkgDefault::Value(double val)
{
  try
  {
    _X = std::sqrt(val);
    if (!std::isfinite(_X))
      throw std::runtime_error("Quiet NaN overflow occurred");
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("Value failed"));
  }
}

double CValueBkgDefault::ValueAt(double atX) const
{
  try
  {
    double ret = square(atX);
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("ValueAt failed"));
  }
}

double CValueBkgDefault::GradAt(double atX) const
{
  try
  {
    double ret = 2.0 * atX;
    if (!std::isfinite(ret))
      throw std::runtime_error("Quiet NaN overflow occurred");
    return ret;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("GradAt failed"));
  }
}

double CValueBkg::X() const
{
  return _X;
}

void CValueBkg::X(double Value)
{
  _X = Value;
}

double CValueBkg::Value() const
{
  return _X;
}

void CValueBkg::Value(double val)
{
  _X = val;
}

double CValueBkg::ValueAt(double atX) const
{
  return atX;
}

double CValueBkg::GradAt(double atX)
{
  return 1;
}

}
