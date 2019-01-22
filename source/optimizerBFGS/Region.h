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
  size_t first_channel, last_channel;

  std::vector<double> current_fit;
  //std::vector<double> fit_gradients; // \todo why unused?
  Eigen::SparseMatrix<double> inv_hessian;

  mutable double chi_squared {0};
  //mutable double UncChisq{0};

  Region(CSpectrum& spe, size_t from_channel, size_t to_channel);

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

  void find_peaks(uint8_t threshold = 3);
  void add_peak(double position, double min, double max, double amplitude = 10.0);
  void remove_peak(size_t index);
  virtual double peak_area(size_t index) const;
  virtual double peak_area_unc(size_t index) const;
  virtual double peak_area_eff(size_t index, const Calibration& cal);
  virtual double peak_area_eff_unc(size_t index, const Calibration& cal);
  virtual size_t fit_var_count() const;
  virtual void setup_fit();
  virtual void store_fit();
  virtual void eval_fit(double pos, std::vector<double>& ret) const;
  virtual double calc_chi_sq(const std::vector<double>& fit) const;
  double chi_sq_normalized() const;
  size_t degrees_of_freedom() const;
  virtual void grad_chi_sq(const std::vector<double>& fit,
                           std::vector<double>& gradients, double& Chisq) const;

  double grad_chi_sq(std::vector<double>& gradients) const;
  double calc_chi_sq() const;

 private:
  // \todo what does this mean?
  static int32_t L(int32_t i, int32_t j, int32_t m);

  // background
  ValueBkgDefault background_base_;

  bool slope_enabled_{true};
  ValueBkg background_slope_;

  bool curve_enabled_{true};
  ValueBkg background_curve_;

  // \todo why skew naming different?
  // peak
  Value width_;

  Value short_tail_amplitude_, short_tail_slope_;

  bool right_tail_enabled_{true};
  Value right_tail_amplitude_, right_tail_slope_;

  // step & tail
  bool left_tail_enabled_{true};
  Value long_tail_amplitude_, long_tail_slope_;

  bool step_enabled_{true};
  Value step_amplitude_;

  std::vector<Peak> peaks_;
  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak
};

}
