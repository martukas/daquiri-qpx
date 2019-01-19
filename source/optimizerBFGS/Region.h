#pragma once

#include <optimizerBFGS/Peak.h>
#include <optimizerBFGS/Spectrum.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

namespace Hypermet
{

class Region
{
 public:
  enum class Type : uint16_t
  {
    Normal = 0,
    Annihilation = 1,
    Boron = 2,
    GeTriangle = 3,
    IntPeakShape = 4
  };

  CSpectrum& spectrum;
  Type region_type{Type::Normal};
  double first_channel, last_channel;

  // background
  ValueBkgDefault background_base;
  ValueBkg background_slope, background_curve;

  // peak
  ValueDefault width;
  Value short_tail_amplitude, short_tail_slope;
  Value right_tail_amplitude, right_tail_slope;

  // step & tail
  Value long_tail_amplitude, long_tail_slope;
  Value step_amplitude;

  std::vector<Peak> peaks;
  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak

  std::vector<double> fit;
  //std::vector<double> fit_gradients; // \todo why unused?
  Eigen::SparseMatrix<double> inv_hessian;

  mutable double chi_squared {0};
  //mutable double UncChisq{0};

  Region(CSpectrum& spe, double from_channel, double to_channel);

  bool slope() const;
  void slope(bool enable);
  bool curve() const;
  void curve(bool enable);
  bool left_tail() const;
  void left_tail(bool enable);
  bool right_tail() const;
  void right_tail(bool enable);
  bool step() const;
  void step(bool enable);

  void find_peaks(uint8_t Threshold = 3);
  void add_peak(double Position, double Min, double Max, double Gamma = 10);
  void remove_peak(size_t index);
  virtual double peak_area(size_t PeakIndex) const;
  virtual double peak_area_unc(size_t PeakIndex) const;
  virtual double peak_area_eff(size_t PeakIndex, const Calibration& cal);
  virtual double peak_area_eff_unc(size_t PeakIndex, const Calibration& cal);
  virtual size_t fit_var_count() const;
  virtual void setup_fit();
  virtual void store_fit();
  virtual void eval_fit(double E, std::vector<double>& ret) const;
  virtual double calc_chi_sq(const std::vector<double>& XVector) const;
  double chi_sq_normalized() const;
  size_t degrees_of_freedom() const;
  virtual void grad_chi_sq(const std::vector<double>& XVector,
                           std::vector<double>& XGradient, double& Chisq) const;

 private:
  // \todo what does this mean?
  static int32_t L(int32_t i, int32_t j, int32_t m);

 protected:
  bool slope_enabled_{true};
  bool curve_enabled_{true};
  bool left_tail_enabled_{true};
  bool step_enabled_{true};
  bool right_tail_enabled_{true};
};

}
