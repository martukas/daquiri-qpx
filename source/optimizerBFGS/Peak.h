#pragma once

#include <optimizerBFGS/Calibration.h>
#include <math.h>

template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
  return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
  return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
  return signum(x, std::is_signed<T>());
}

template <typename T> inline constexpr
T square(T x) {
  return std::pow(x, 2);
}

template <typename T> inline constexpr
T cube(T x) {
  return std::pow(x, 3);
}

class CValueDefault {
 private:
  double _X {0};
  double _dX {0};
  double _Max {1};
  double _Min {0};
 public:
  double UncValue {0};
  int32_t XIndex {-1};

  //void SetRange(double NewMin, double NewMax) {
  //    if(Not NewMin < NewMax) { return; }
  //    double Value = _Min + (_Max - _Min) / 2 * (1 + std::sin(_X));
  //    Value = std::min(NewMin, Value);
  //    Value = std::max(NewMax, Value);
  //    _Min = NewMin;
  //    _Max = NewMax;
  //    _X = std::asin((_Min + _Max - 2 * Value) / (-1 * _Max + _Min));
  //}

  double Max() const;
  void Max(double NewMax);
  double Min();
  void Min(double NewMin);
  double X();
  void X(double Value);
  double Value() const;
  void Value (double val);
  double ValueAt(double atX);
  double GradAt(double atX);
};

class CValue: public CValueDefault {
 public:
  bool ToFit {true};
};


class CValueGam
{
 private:
  double _X{0};
  double _dX{0};
 public:
  double UncValue{0};
  int32_t XIndex{-1};

 public:
  double X();
  void X(double Value);
  double Value() const;
  void Value(double val);
  double ValueAt(double atX);
  double GradAt(double atX);
};

class CPeak
{
 private:
  int32_t FEPStatus{1};
 public:
  CValueDefault POS;
  CValueGam GAM;

  int32_t StepType() const;
  double PeakPosition() const;
  double UncPeakPosition() const;
  double PeakEnergy(const CCalibration& cal) const;
  double UncPeakEnergy(const CCalibration& cal) const;
  bool IsFullEnergyPeak() const;
  void IsFullEnergyPeak(bool Value);
  bool operator<(const CPeak& other) const;
};

class CValueBkgDefault
{
 private:
  double _X{0};
  double _dX{0};
 public:
  double UncValue{0};
  double XIndex{-1};

  double X() const;
  void X(double Value);
  double Value() const;
  void Value(double val);
  double ValueAt(double atX);
  double GradAt(double atX);
};

class CValueBkg
{
 private:
  double _X{0};
  double _dX{0};
 public:
  double UncValue{0};
  int32_t XIndex{-1};
  bool ToFit{true};

  double X() const;
  void X(double Value);

  double Value() const;

  void Value(double val);

  double ValueAt(double atX) const;

  double GradAt(double atX);
};
