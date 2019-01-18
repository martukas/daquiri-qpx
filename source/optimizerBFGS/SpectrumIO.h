#include <cstdint>
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

class CEff
{
  // <summary>
  // This class loads and calculates the efficiency from an existing fit, as done by Hypermet-PC.
  // See also in:
  // G.L. Molnar, Zs. Revay, T. Belgya:
  // Wide energy range efficiency calibration method for Ge detectors,
  // Nucl. Instr. Meth. A 489 (2002) 140?159
  // </summary>
  // <remarks></remarks>
 private:
  float e_maxdeg;
  double e_apol[7], e_bpol[6], e_normfact[8], e_poly_coeff[8];

  Eigen::SparseMatrix<double> e_VarMatrix{8, 8};

  double e_c0, e_c1;
  bool e_init_done{false};

 public:
  void Init(std::string flnm);
  void Close();
  bool InitDone() const;
  double Value(double Energy) const;
  double SigmaRel(double ByVal);

 private:
  double e_ortpol(size_t n, double X) const;

};

class CNonlin
{
  // <summary>
  // This class loads and calculates the nonlinearity from an existing fit, as done by Hypermet-PC.
  // See also in:
  // B. Fazekas, Zs. Revay, J. ?st?r, T. Belgya, G. Moln?r, A. Simonits:
  // A new method for determination of gamma-ray spectrometer non-linearity,
  // Nucl. Instr. Meth. A 422 (1999) 469-473
  // </summary>
  // <remarks></remarks>
 private:
  float n_maxdeg;
  double n_apol[7], n_bpol[6], n_normfact[8], n_poly_coeff[7];
  Eigen::SparseMatrix<double> n_VarMatrix{8, 8};

  double n_c0, n_c1;
  double n_bl0, n_bl1;
  bool n_init_done{false};

 public:
  void Init(std::string flnm);
  void Close();
  bool InitDone() const;
  double Value(double Position) const;
  double Sigma(double Position);
  void SetBasePoints(double ch1, double ch2);

 private:
  double n_ortpol(size_t n, double X) const;
  double nonlin1(double Position);
};

class CCalibration
{
 public:
  struct CalPoint
  {
    float Channel;
    float UncChannel;
    float Value;
    float UncValue;
  };

  std::string ChannelUnits;
  std::vector<CalPoint> EnergyCal{2};
  std::vector<CalPoint> WidthCal{2};
  std::string Source;
  CNonlin Nonlinearity;
  CEff Efficiency;

 public:
  void New();
  uint8_t CalOrder() const;
  void CalOrder(uint8_t Value);
  double EnergyConst() const;
  double EnergySlope() const;
  double ChannelToEnergy(double Channel) const;
  double EnergyToChannel(double Energy) const;
  double Width(double Channel) const;

 private:
  uint8_t _CalOrder;
};

class CSpectrum
{
 public:
  std::vector<size_t> Channel;
  CCalibration Calibration;
  double Weight(size_t i) const;

  // template type for Val;
  int8_t Sign(double Val);

  double DeadTime(double TrueTime, double LiveTime);
  double Rate(double LiveTime, double SumCounts);

 private:
  // template type for Val;
  size_t mystery_function(double Val);
};

