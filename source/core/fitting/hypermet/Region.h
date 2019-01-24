#pragma once

#include <core/fitting/hypermet/Peak.h>
#include <core/fitting/finder.h>

#include <core/fitting/BFGS/Fittable.h>

namespace Hypermet
{

class Region : public Fittable
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
  Type region_type{Type::Normal};

  DAQuiri::Finder finder_;           // gets x & y data from fitter

  Region(const DAQuiri::Finder& finder);

  void add_peak(double position, double min, double max, double amplitude = 10.0);
  void remove_peak(size_t index);
  double peak_area(size_t index) const;
  double peak_area_unc(size_t index) const;
  double peak_area_eff(size_t index, const Calibration& cal);
  double peak_area_eff_unc(size_t index, const Calibration& cal);


  void map_fit();
  void save_fit_uncerts(const FitResult& result);

  // Fittable implementation
  std::vector<double> variables() const override;
  double degrees_of_freedom() const  override;
  double chi_sq(const std::vector<double>& fit) const override;
  double grad_chi_sq(const std::vector<double>& fit,
                     std::vector<double>& gradients) const override;

 private:
  int32_t var_count_{0};

  // background
  ValueGam background_base_;

  bool slope_enabled_{true};
  ValueBkg background_slope_;

  bool curve_enabled_{true};
  ValueBkg background_curve_;

  // peak
  Peak default_peak_;
  std::vector<Peak> peaks_;
  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak

  size_t fit_var_count() const;
  double chi_sq_normalized() const;
  void save_fit(const std::vector<double>& variables);
  double chi_sq() const;
  double grad_chi_sq(std::vector<double>& gradients) const;

};

}
