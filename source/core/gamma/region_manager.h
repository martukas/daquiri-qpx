#pragma once

#include <core/gamma/hypermet/gamma_region.h>
#include <core/gamma/fit_evaluation.h>
#include <core/gamma/fit_settings.h>

#include <core/fitting/optimizers/abstract_optimizer.h>

namespace DAQuiri {

struct FitDescription
{
  std::string description;
  size_t    peaknum {0};
  double chi_sq_norm {0};
  double sum4aggregate {0};
};

class Fit {
 public:
  Fit(const Region& r, std::string descr);

  FitDescription description;
  Region region;
};

struct PeakRendering
{
  std::vector<double>
      peak,
      full_fit;

  void clear();
  void render(const Peak&);
};

struct RegionRendering
{
  uint8_t subdivisions {10};

  std::vector<double>
      channel,
      energy,
      background,
      back_steps,
      full_fit,
      sum4_background;

  std::map<double, PeakRendering> peaks;

  void reserve(size_t count);
  void clear();

  void render(const Region& r, const Calibration& energy_calib);
};


class RegionManager {
 public:
  RegionManager() = default;
  RegionManager(const nlohmann::json& j, const FitEvaluation &finder, const FitSettings& fs);
  RegionManager(const FitSettings& fs, const FitEvaluation &parentfinder, double min, double max);

  double id() const;

  //access peaks
  size_t peak_count() const;
  bool contains(double peakID) const;
  Peak peak(double peakID) const;

  // current state
  const Region& region() const;
  void modify_region(const Region& new_region, std::string message = "");
  //access history
  size_t current_fit() const;
  std::vector<FitDescription> history() const;
  bool rollback(size_t i);

  //access other
  FitSettings fit_settings() const { return settings_; }
  const FitEvaluation &finder() const { return fit_eval_; }
  const RegionRendering &rendering() const { return rendering_; }

  //manupulation, may invoke optimizer
  bool refit(AbstractOptimizer* optimizer);

  nlohmann::json to_json(const FitEvaluation &parent_finder) const;

private:
  FitSettings settings_;
  FitEvaluation fit_eval_;           // gets x & y data from fitter
  RegionRendering rendering_;

  //history
  std::vector<Fit> fits_;
  size_t current_fit_ {0};

  //intrinsic, these are saved
  Region region_;

  void set_data(const FitEvaluation &parentfinder, double min, double max);

//  bool add_from_resid(AbstractOptimizer* optimizer);
//  void iterative_fit(AbstractOptimizer* optimizer);

  void render();
  void save_current_fit(std::string description);
};

}
