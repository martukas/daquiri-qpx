#pragma once

#include <optimizerBFGS/Calibration.h>

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
  double ValueAt(double atX) const;
  double GradAt(double atX) const;
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
  double ValueAt(double atX) const;
  double GradAt(double atX) const;
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
  double ValueAt(double atX) const;
  double GradAt(double atX) const;
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
