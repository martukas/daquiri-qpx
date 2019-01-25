#pragma once

#include <core/fitting/peak.h>
#include <core/fitting/finder.h>
//#include <core/fitting/optimizer.h>
#include <atomic>

#include <core/fitting/BFGS/Fittable.h>
#include <core/fitting/BFGS/BFGS.h>


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

class Fit : public Hypermet::Fittable {
 public:
  Fit(const SUM4Edge &lb, const SUM4Edge &rb,
      const std::map<double, Peak> &peaks,
      const Finder &finder,
      std::string descr);

  const DAQuiri::Finder &finder_;

  FitDescription description;

  std::map<double, Peak> peaks_;
  SUM4Edge  LB_, RB_;
//  Polynomial background_;
  FitSettings settings_;


  int32_t var_count_{0};

  // background
  Hypermet::ValueGam background_base_;
  bool slope_enabled_{true};
  Hypermet::ValueBkg background_slope_;
  bool curve_enabled_{true};
  Hypermet::ValueBkg background_curve_;

  // peak
  Peak default_peak_;

//  void add_peak(double position, double min, double max, double amplitude = 10.0);
  double peak_area(size_t index) const;
  double peak_area_unc(size_t index) const;
  double peak_area_eff(size_t index, const Hypermet::Calibration& cal);
  double peak_area_eff_unc(size_t index, const Hypermet::Calibration& cal);


  void map_fit();
  void save_fit_uncerts(const Hypermet::FitResult& result);

  // Fittable implementation
  std::vector<double> variables() const override;
  double degrees_of_freedom() const  override;
  double chi_sq(const std::vector<double>& fit) const override;
  double grad_chi_sq(const std::vector<double>& fit,
                     std::vector<double>& gradients) const override;

 private:
  size_t fit_var_count() const;
  double chi_sq_normalized() const;
  void save_fit(const std::vector<double>& variables);
  double chi_sq() const;
  double grad_chi_sq(std::vector<double>& gradients) const;

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
  bool replace_hypermet(double &peakID, Hypermet::Peak hyp);
  bool override_energy(double peakID, double energy);

  //manupulation, may invoke optimizer
  bool auto_fit(Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  bool refit(Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  bool adjust_LB(const Finder &parentfinder, double left, double right,
                 Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  bool adjust_RB(const Finder &parentfinder, double left, double right,
                 Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  bool add_peak(const Finder &parentfinder, double left, double right,
                Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  bool remove_peaks(const std::set<double> &pks, Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
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

  bool add_from_resid(Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor,
                      int32_t centroid_hint = -1);
  bool rebuild(Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  bool rebuild_as_hypermet(Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  bool rebuild_as_gaussian(Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);
  void iterative_fit(Hypermet::BFGS& optimizer, std::atomic<bool>& interruptor);

  void render();
  void save_current_fit(std::string description);
};

}
