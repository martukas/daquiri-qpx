#pragma once

#include <core/fitting/peak.h>
#include <core/fitting/finder.h>
#include <atomic>
#include <set>
#include <core/fitting/BFGS/BFGS.h>
#include <core/fitting/hypermet/PolyBackground.h>

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
  Fit(const SUM4Edge &lb, const SUM4Edge &rb,
      PolyBackground bkg,
      const std::map<double, Peak> &peaks,
//      const Finder &finder,
      std::string descr);

//  const DAQuiri::Finder &finder_;

  FitDescription description;

  //FitSettings settings_;

  SUM4Edge  LB_, RB_;
  PolyBackground background;
//  Polynomial background_;


  Peak default_peak_;
  std::map<double, Peak> peaks_;
};

struct PeakRendering
{
  std::vector<double>
      peak,
      full_fit;

  void clear();
  void render(const Hypermet&);
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
  void sum4only(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const Calibration& energy_calib,
                   const Polynomial& sum4back,
                   const std::map<double, Peak>& pks);

  void with_hypermet(double start, double end,
                     const Calibration& energy_calib,
                     const Polynomial& sum4back,
                     const PolyBackground& hyp_back,
                     const std::map<double, Peak>& pks);
};


struct ROI {
  ROI() = default;
  ROI(const nlohmann::json& j,
      const Finder &finder,
      const FitSettings& fs);
  ROI(const FitSettings& fs, const Finder &parentfinder, double min, double max);

  //bounds
  double ID() const;
  double left_bin() const;
  double right_bin() const;
  double width() const;

  bool overlaps(double bin) const;
  bool overlaps(double Lbin, double Rbin) const;
  bool overlaps(const ROI& other) const;

  //access peaks
  const std::map<double, Peak> &peaks() const;

  //access other
  SUM4Edge LB() const {return LB_;}
  SUM4Edge RB() const {return RB_;}
  FitSettings fit_settings() const { return settings_; }
  const Finder &finder() const { return finder_; }

  //access history
  size_t current_fit() const;
  std::vector<FitDescription> history() const;

  //manipulation, no optimizer
  bool rollback(const Finder &parent_finder, size_t i);
  bool adjust_sum4(double &peakID, double left, double right);
  bool replace_hypermet(double &peakID, Hypermet hyp);
  bool override_energy(double peakID, double energy);

  //manupulation, may invoke optimizer
  bool auto_fit(BFGS& optimizer, std::atomic<bool>& interruptor);
  bool refit(BFGS& optimizer, std::atomic<bool>& interruptor);
  bool adjust_LB(const Finder &parentfinder, double left, double right,
                 BFGS& optimizer, std::atomic<bool>& interruptor);
  bool adjust_RB(const Finder &parentfinder, double left, double right,
                 BFGS& optimizer, std::atomic<bool>& interruptor);
  bool add_peak(const Finder &parentfinder, double left, double right,
                BFGS& optimizer, std::atomic<bool>& interruptor);
  bool remove_peaks(const std::set<double> &pks, BFGS& optimizer, std::atomic<bool>& interruptor);
  bool override_settings(const FitSettings &fs, std::atomic<bool>& interruptor);

  nlohmann::json to_json(const Finder &parent_finder) const;

private:
  FitSettings settings_;
  Finder finder_;           // gets x & y data from fitter
  Hypermet default_peak_;

  //history
  std::vector<Fit> fits_;
  size_t current_fit_ {0};

  //intrinsic, these are saved
  SUM4Edge LB_, RB_;
  PolyBackground background_;
  std::map<double, Peak> peaks_;

  RegionRendering rendering_;

  void set_data(const Finder &parentfinder, double min, double max);

  std::vector<double> remove_background();
  void init_edges();
  void init_LB();
  void init_RB();
  void init_background();

  void cull_peaks();
  bool remove_peak(double bin);

  bool add_from_resid(BFGS& optimizer, std::atomic<bool>& interruptor,
                      int32_t centroid_hint = -1);
  bool rebuild(BFGS& optimizer, std::atomic<bool>& interruptor);
  bool rebuild_as_hypermet(BFGS& optimizer, std::atomic<bool>& interruptor);
  bool rebuild_as_gaussian(BFGS& optimizer, std::atomic<bool>& interruptor);
  void iterative_fit(BFGS& optimizer, std::atomic<bool>& interruptor);

  void render();
  void save_current_fit(std::string description);
};

}
