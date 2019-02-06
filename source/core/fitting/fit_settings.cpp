#include <core/fitting/fit_settings.h>

namespace DAQuiri
{


void to_json(nlohmann::json& j, const KONSettings& s)
{
  j["width"] = s.width;
  j["sigma_spectrum"] = s.sigma_spectrum;
  j["sigma_resid"] = s.sigma_resid;
  j["edge_width_factor"] = s.edge_width_factor;
}

void from_json(const nlohmann::json& j, KONSettings& s)
{
  s.width = j["width"];
  s.sigma_spectrum = j["sigma_spectrum"];
  s.sigma_resid = j["sigma_resid"];
  s.edge_width_factor = j["edge_width_factor"];
}

double FCalibration::nrg_to_bin(double energy) const
{
  return cali_nrg_.inverse(energy, 0.1);
}

double FCalibration::bin_to_nrg(double bin) const
{
  return cali_nrg_.transform(bin);
}

double FCalibration::bin_to_width(double bin) const
{
  double nrg = bin_to_nrg(bin);
  double fwhm = nrg_to_fwhm(nrg);
  return (nrg_to_bin(nrg + 0.5 * fwhm) - nrg_to_bin(nrg - 0.5 * fwhm));
}

double FCalibration::nrg_to_fwhm(double energy) const
{
  if (cali_fwhm_.valid())
    return cali_fwhm_.transform(energy);
  else
    return 1;
}

void FitSettings::clear()
{
  //bits_ = 0;
  live_time = hr_duration_t();
}
void FitSettings::clone(FitSettings other)
{
  other.calib = calib;
  //other.bits_ = bits_;
  other.live_time = live_time;
  (*this) = other;
}

void to_json(nlohmann::json& j, const FitSettings& s)
{
  j["KON"] = s.kon_settings;

  j["ROI"]["max_peaks"] = s.ROI_max_peaks;
  j["ROI"]["extend_background"] = s.ROI_extend_background;
  j["ROI"]["edge_samples"] = s.background_edge_samples;

  j["residuals"]["auto"] = s.resid_auto;
  j["residuals"]["max_iterations"] = s.resid_max_iterations;
  j["residuals"]["min_amplitude"] = s.resid_min_amplitude;
  j["residuals"]["too_close"] = s.resid_too_close;

  j["small_peaks"]["enabled"] = s.small_simplify;
  j["small_peaks"]["max_amplitude"] = s.small_max_amplitude;

  j["width_options"]["use_common"] = s.width_common;
  j["width_options"]["exception511"] = s.width_at_511_variable;
  j["width_options"]["tolerance511"] = s.width_at_511_tolerance;
  j["width_options"]["common_bounds"] = s.width_common_bounds;

  j["hypermet_peak"] = s.default_peak;
  j["fitter_max_iterations"] = s.fitter_max_iter;
}

void from_json(const nlohmann::json& j, FitSettings& s)
{
  s.kon_settings = j["KON"];

  s.ROI_max_peaks = j["ROI"]["max_peaks"];
  s.ROI_extend_background = j["ROI"]["extend_background"];
  s.background_edge_samples = j["ROI"]["edge_samples"];

  s.resid_auto = j["residuals"]["auto"];
  s.resid_max_iterations = j["residuals"]["max_iterations"];
  s.resid_min_amplitude = j["residuals"]["min_amplitude"];
  s.resid_too_close = j["residuals"]["too_close"];

  s.small_simplify = j["small_peaks"]["enabled"];
  s.small_max_amplitude = j["small_peaks"]["max_amplitude"];

  s.width_common = j["width_options"]["use_common"];
  s.width_at_511_variable = j["width_options"]["exception511"];
  s.width_at_511_tolerance = j["width_options"]["tolerance511"];
  s.width_common_bounds = j["width_options"]["common_bounds"];

  s.fitter_max_iter = j["fitter_max_iterations"];
  s.default_peak = j["hypermet_peak"];
}

}
