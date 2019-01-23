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

  Region(CSpectrum& spe, size_t from_channel, size_t to_channel);

  void find_peaks(uint8_t threshold = 3);
  void add_peak(double position, double min, double max, double amplitude = 10.0);
  void remove_peak(size_t index);
  double peak_area(size_t index) const;
  double peak_area_unc(size_t index) const;
  double peak_area_eff(size_t index, const Calibration& cal);
  double peak_area_eff_unc(size_t index, const Calibration& cal);
  size_t fit_var_count() const;
  double chi_sq_normalized() const;
  size_t degrees_of_freedom() const;

  void map_fit();
  void load_fit();
  void save_fit();

  void save_fit_uncerts();

  double calc_chi_sq() const;
  double grad_chi_sq(std::vector<double>& gradients) const;

  double calc_chi_sq_at(const std::vector<double>& fit) const;
  double grad_chi_sq_at(const std::vector<double>& fit,
                        std::vector<double>& gradients) const;

 private:
  // \todo what does this mean?
  static int32_t L(int32_t i, int32_t j, int32_t m);

  // background
  ValueGam background_base_;

  bool slope_enabled_{true};
  ValueBkg background_slope_;

  bool curve_enabled_{true};
  ValueBkg background_curve_;

  // \todo why skew naming different?
  // peak
  Peak default_peak_;
  std::vector<Peak> peaks_;
  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak
};

}
