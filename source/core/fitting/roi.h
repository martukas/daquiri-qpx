#pragma once

#include <core/fitting/peak.h>
#include <core/fitting/finder.h>
#include <core/fitting/optimizer.h>
#include <atomic>

namespace DAQuiri {

struct FitDescription
{
  std::string description;
  int    peaknum;
  double rsq;
  double sum4aggregate;

  FitDescription()
    : peaknum(0), rsq(0), sum4aggregate(0) {}
};

struct Fit {
  Fit(const SUM4Edge &lb, const SUM4Edge &rb,
      const Polynomial &backg,
      const std::map<double, Peak> &peaks,
      const Finder &finder,
      std::string descr);

  FitDescription description;

  std::map<double, Peak> peaks_;
  SUM4Edge  LB_, RB_;
  Polynomial background_;
  FitSettings settings_;
};


struct ROI {
  ROI() = default;
  ROI(const nlohmann::json& j, const Finder &finder);
  ROI(const Finder &parentfinder, double min, double max);

  //bounds
  double ID() const;
  double left_bin() const;
  double right_bin() const;
  double width() const;
  double left_nrg() const;  //may return NaN
  double right_nrg() const; //may return NaN

  bool overlaps(double bin) const;
  bool overlaps(double Lbin, double Rbin) const;
  bool overlaps(const ROI& other) const;

  //access peaks
  size_t peak_count() const;
  bool contains(double peakID) const;
  Peak peak(double peakID) const;
  const std::map<double, Peak> &peaks() const;

  //access other
  SUM4Edge LB() const {return LB_;}
  SUM4Edge RB() const {return RB_;}
  FitSettings fit_settings() const { return finder_.settings_; }
  const Finder &finder() const { return finder_; }

  //access history
  size_t current_fit() const;
  size_t history_size() const;
  std::vector<FitDescription> history() const;

  //manipulation, no optimizer
  bool rollback(const Finder &parent_finder, size_t i);
  bool adjust_sum4(double &peakID, double left, double right);
  bool replace_hypermet(double &peakID, Hypermet hyp);
  bool override_energy(double peakID, double energy);

  //manupulation, may invoke optimizer
  bool auto_fit(OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool refit(OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool adjust_LB(const Finder &parentfinder, double left, double right,
                 OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool adjust_RB(const Finder &parentfinder, double left, double right,
                 OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool add_peak(const Finder &parentfinder, double left, double right,
                OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool remove_peaks(const std::set<double> &pks, OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool override_settings(const FitSettings &fs, std::atomic<bool>& interruptor);

  nlohmann::json to_json(const Finder &parent_finder) const;

  //as rendered for graphing
  std::vector<double>
      hr_x,
      hr_x_nrg,
      hr_background,
      hr_back_steps,
      hr_fullfit,
      hr_sum4_background_;

private:
  //intrinsic, these are saved
  SUM4Edge LB_, RB_;
  Polynomial background_;
  std::map<double, Peak> peaks_;

  Finder finder_;           // gets x & y data from fitter

  //history
  std::vector<Fit> fits_;
  size_t current_fit_ {0};


  void set_data(const Finder &parentfinder, double min, double max);

  std::vector<double> remove_background();
  void init_edges();
  void init_LB();
  void init_RB();
  void init_background();

  void cull_peaks();
  bool remove_peak(double bin);

  bool add_from_resid(OptimizerPtr optimizer, std::atomic<bool>& interruptor,
                      int32_t centroid_hint = -1);
  bool rebuild(OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool rebuild_as_hypermet(OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  bool rebuild_as_gaussian(OptimizerPtr optimizer, std::atomic<bool>& interruptor);
  void iterative_fit(OptimizerPtr optimizer, std::atomic<bool>& interruptor);

  void render();
  void save_current_fit(std::string description);
};

}
