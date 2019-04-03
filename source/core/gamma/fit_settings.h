#pragma once

#include <core/calibration/calibration.h>
#include <core/gamma/hypermet/peak.h>

namespace DAQuiri
{

struct KONSettings
{
  uint16_t width{4};
  double sigma_spectrum{3.0};
  double sigma_resid{3.0};
  double edge_width_factor{3.5};
};

void to_json(nlohmann::json& j, const KONSettings& s);
void from_json(const nlohmann::json& j, KONSettings& s);

struct FCalibration
{
  //specific to spectrum
  Calibration cali_nrg_, cali_fwhm_;

  // \todo with uncerts
  double nrg_to_bin(double energy) const;
  double bin_to_nrg(double bin) const;
  double bin_to_width(double bin) const;
  double nrg_to_fwhm(double energy) const;
};


class FitSettings
{
 public:
  bool overriden{false};

  // \todo should be in KON settings? not serialized; not really used
  double finder_cutoff_kev{100};

  KONSettings kon_settings;

  uint16_t ROI_max_peaks{10};
  double ROI_extend_background{0.6};
  uint16_t background_edge_samples{7};

  bool resid_auto{true};
  uint16_t resid_max_iterations{5};
  uint64_t resid_min_amplitude{5};
  double resid_too_close{0.2};

  bool small_simplify{true};
  uint64_t small_max_amplitude{500};

  bool width_common{true};
  SineBoundedParam width_common_bounds; //{1.0, 0.7, 1.3};
  bool width_at_511_variable{true};
  double width_at_511_tolerance{5.0};

  //hypermet
  Peak default_peak;

  uint16_t fitter_max_iter{3000};

  //specific to spectrum
  FCalibration calib;
//  uint16_t bits_ {0};
  hr_duration_t real_time, live_time;

  FitSettings() = default;
  void clone(FitSettings other);
  void clear();
};

void to_json(nlohmann::json& j, const FitSettings& s);
void from_json(const nlohmann::json& j, FitSettings& s);

}