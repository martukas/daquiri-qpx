#pragma once

#include <optimizerBFGS/Peak.h>
#include <optimizerBFGS/Spectrum.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

enum class RegTypes : uint16_t
{
  Normal = 0,
  Annihilation = 1,
  Boron = 2,
  GeTriangle = 3,
  IntPeakShape = 4
};

class Region
{
 public:
  RegTypes Type{RegTypes::Normal};
  double FirstChannel, LastChannel;
  double Chisq, UncChisq;

 protected:
  bool TailFlag{true};
  bool SlopeFlag{true};
  bool StepFlag{true};
  bool CurveFlag{true};
  bool RightTailFlag{true};

 public:
  CValueDefault DEL;

  CValue AST, BST, ALT, BLT;

  CValue ART, BRT;

  CValue SIG;

  CValueBkgDefault BLN;
  CValueBkg BSL, BCV;

  std::vector<CPeak> Peak;

  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak

 public:
  std::vector<double> Vector{4};
  std::vector<double> Gradient{4};

  Eigen::SparseMatrix<double> Hessinv{4, 4};

  long nfuncEval;

  void New(double FromChannel,
           double ToChannel);
  //, Position As Double, Min As Double, Max As Double, Optional Gamma As Double = 10)

  bool LeftTail() const;
  void LeftTail(bool Value);
  bool RightTail() const;
  void RightTail(bool Value);
  bool Slope() const;
  void Slope(bool Value);
  bool Curve();
  void Curve(bool Value);
  bool StepBkg() const;
  void StepBkg(bool Value);

  void SearchPeaks(const CSpectrum& spectrum, uint8_t Threshold = 3);
  void AddPeak(double Position, double Min, double Max, double Gamma = 10);
  void DeletePeak(size_t index);
  void SortPeak();
  virtual double PeakArea(size_t PeakIndex) const;
  virtual double UncPeakArea(size_t PeakIndex);
  virtual double PeakAreaEff(size_t PeakIndex, const CCalibration& cal);
  virtual double UncPeakAreaEff(size_t PeakIndex, const CCalibration& cal);
  virtual size_t NPeaks() const;
  virtual size_t FitVars() const;
  virtual void setupFit();
  virtual void storeFit();
  virtual void FuncValue(double E, std::vector<double>& Value);
  virtual double CalcChiSq(const CSpectrum& spectrum,
                           const std::vector<double>& XVector);
  double ChisqNorm() const;
  size_t DegreeOfFreedom() const;
  virtual void GradChiSq(const CSpectrum& spectrum,
                         const std::vector<double>& XVector,
                         std::vector<double>& XGradient, double& Chisq);

 private:
  int32_t L(int32_t i, int32_t j, int32_t m);

};
