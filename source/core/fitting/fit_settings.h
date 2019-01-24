#pragma once

#include <core/calibration/calibration.h>

namespace DAQuiri
{

class FitSettings
{
 public:
  bool overriden{false};

  double finder_cutoff_kev{100};

  uint16_t KON_width{4};
  double KON_sigma_spectrum{3.0};
  double KON_sigma_resid{3.0};

  uint16_t ROI_max_peaks{10};
  double ROI_extend_peaks{3.5};
  double ROI_extend_background{0.6};
  uint16_t background_edge_samples{7};
  bool sum4_only{false};

  bool resid_auto{true};
  uint16_t resid_max_iterations{5};
  uint64_t resid_min_amplitude{5};
  double resid_too_close{0.2};

  bool small_simplify{true};
  uint64_t small_max_amplitude{500};

  bool width_common{true};
  Parameter width_common_bounds{1.0, 0.7, 1.3};
  bool width_at_511_variable{true};
  double width_at_511_tolerance{5.0};

  //hypermet
  bool gaussian_only{false};
  double lateral_slack{0.5};
  Parameter width_variable_bounds{1.0, 0.7, 4};
  //variable bounds
  bool step_enabled{true};
  Parameter step_amplitude{1.0e-10, 1.0e-10, 0.75};
  bool tail_enabled{true};
  Parameter tail_amplitude{1.0e-10, 1.0e-10, 0.015};
  Parameter tail_slope{2.75, 2.5, 50};
  bool Lskew_enabled{true};
  Parameter Lskew_amplitude{1.0e-10, 1.0e-10, 0.75};
  Parameter Lskew_slope{0.5, 0.3, 2};
  bool Rskew_enabled{true};
  Parameter Rskew_amplitude{1.0e-10, 1.0e-10, 0.75};
  Parameter Rskew_slope{0.5, 0.3, 2};
  uint16_t fitter_max_iter{3000};

  //specific to spectrum
  Calibration cali_nrg_, cali_fwhm_;
//  uint16_t bits_ {0};
  hr_duration_t real_time, live_time;

  FitSettings() = default;
  void clone(FitSettings other);
  void clear();

  double nrg_to_bin(double energy) const;
  double bin_to_nrg(double bin) const;
  double bin_to_width(double bin) const;
  double nrg_to_fwhm(double energy) const;
};

void to_json(nlohmann::json& j, const FitSettings& s);
void from_json(const nlohmann::json& j, FitSettings& s);

}