#include <core/fitting/fit_settings.h>

namespace DAQuiri
{

void FitSettings::clear()
{
  cali_nrg_ = Calibration();
  cali_fwhm_ = Calibration();
  //bits_ = 0;
  live_time = hr_duration_t();
}
void FitSettings::clone(FitSettings other)
{
  other.cali_nrg_ = cali_nrg_;
  other.cali_fwhm_ = cali_fwhm_;
  //other.bits_ = bits_;
  other.live_time = live_time;
  (*this) = other;
}

double FitSettings::nrg_to_bin(double energy) const
{
  return cali_nrg_.inverse(energy, 0.1);
}

double FitSettings::bin_to_nrg(double bin) const
{
  return cali_nrg_.transform(bin);
}

double FitSettings::bin_to_width(double bin) const
{
  double nrg = bin_to_nrg(bin);
  double fwhm = nrg_to_fwhm(nrg);
  return (nrg_to_bin(nrg + fwhm / 2.0) - nrg_to_bin(nrg - fwhm / 2.0));
}

double FitSettings::nrg_to_fwhm(double energy) const
{
  if (cali_fwhm_.valid())
    return cali_fwhm_.transform(energy);
  else
    return 1;
}

void to_json(nlohmann::json& j, const FitSettings& s)
{
  j["KON"]["width"] = s.KON_width;
  j["KON"]["sigma_spectrum"] = s.KON_sigma_spectrum;
  j["KON"]["sigma_resid"] = s.KON_sigma_resid;

  j["ROI"]["max_peaks"] = s.ROI_max_peaks;
  j["ROI"]["extend_peaks"] = s.ROI_extend_peaks;
  j["ROI"]["extend_background"] = s.ROI_extend_background;
  j["ROI"]["edge_samples"] = s.background_edge_samples;
  j["ROI"]["sum4_only"] = s.sum4_only;

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

  j["hypermet"]["gaussian_only"] = s.gaussian_only;
  j["hypermet"]["lateral_slack"] = s.lateral_slack;
  j["hypermet"]["fitter_max_iterations"] = s.fitter_max_iter;
  j["hypermet"]["width_variable_bounds"] = s.width_variable_bounds;
  j["hypermet"]["step_amplitude"] = s.step_amplitude;
  j["hypermet"]["tail_amplitude"] = s.tail_amplitude;
  j["hypermet"]["tail_slope"] = s.tail_slope;
  j["hypermet"]["Lskew_amplitude"] = s.Lskew_amplitude;
  j["hypermet"]["Lskew_slope"] = s.Lskew_slope;
  j["hypermet"]["Rskew_amplitude"] = s.Rskew_amplitude;
  j["hypermet"]["Rskew_slope"] = s.Rskew_slope;
}

void from_json(const nlohmann::json& j, FitSettings& s)
{
  s.KON_width = j["KON"]["width"];
  s.KON_sigma_spectrum = j["KON"]["sigma_spectrum"];
  s.KON_sigma_resid = j["KON"]["sigma_resid"];

  s.ROI_max_peaks = j["ROI"]["max_peaks"];
  s.ROI_extend_peaks = j["ROI"]["extend_peaks"];
  s.ROI_extend_background = j["ROI"]["extend_background"];
  s.background_edge_samples = j["ROI"]["edge_samples"];
  s.sum4_only = j["ROI"]["sum4_only"];

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

  s.gaussian_only = j["hypermet"]["gaussian_only"];
  s.lateral_slack = j["hypermet"]["lateral_slack"];
  s.fitter_max_iter = j["hypermet"]["fitter_max_iterations"];
  s.width_variable_bounds = j["hypermet"]["width_variable_bounds"];
  s.step_amplitude = j["hypermet"]["step_amplitude"];
  s.tail_amplitude = j["hypermet"]["tail_amplitude"];
  s.tail_slope = j["hypermet"]["tail_slope"];
  s.Lskew_amplitude = j["hypermet"]["Lskew_amplitude"];
  s.Lskew_slope = j["hypermet"]["Lskew_slope"];
  s.Rskew_amplitude = j["hypermet"]["Rskew_amplitude"];
  s.Rskew_slope = j["hypermet"]["Rskew_slope"];

}

}
