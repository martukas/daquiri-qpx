#pragma once

#include <core/fitting/region_manager.h>
#include <core/consumer.h>

namespace DAQuiri {

class Fitter {
  
public:
  Fitter() = default;
  Fitter(ConsumerPtr spectrum)
  { setData(spectrum); }

  Fitter(const json& j, ConsumerPtr spectrum);

  FitSettings settings() const { return settings_; }
  void apply_settings(FitSettings settings);
  const FitEvaluation &finder() const { return finder_; }
//  void apply_energy_calibration(Calibration cal);
//  void apply_fwhm_calibration(Calibration cal);

  bool empty() const;

  void clear();
  void setData(ConsumerPtr spectrum);
  void find_regions();

  //access peaks
  size_t peak_count() const;
  bool contains_peak(double bin) const;
  Peak peak(double peakID) const;
  std::map<double, Peak> peaks() const;

  //access regions
  size_t region_count() const;
  bool contains_region(double bin) const;
  ROI region(double bin) const;
  const std::map<double, ROI> &regions() const;
  ROI parent_region(double peakID) const;
  std::set<double> relevant_regions(double left, double right);

  //manupulation, may invoke optimizer
  bool find_and_fit(double regionID, BFGS& optimizer);
  bool add_peak(double left, double right, BFGS& optimizer);
  bool adj_LB(double regionID, double left, double right, BFGS& optimizer);
  bool adj_RB(double regionID, double left, double right, BFGS& optimizer);
  bool merge_regions(double left, double right, BFGS& optimizer);
  bool refit_region(double regionID, BFGS& optimizer);
  bool override_ROI_settings(double regionID, const FitSettings &fs);
  bool remove_peaks(std::set<double> peakIDs, BFGS& optimizer);
  //manipulation, no optimizer
  bool adjust_sum4(double &peakID, double left, double right);
  bool replace_hypermet(double &peakID, Peak hyp);
  bool rollback_ROI(double regionID, size_t point);
  bool delete_ROI(double regionID);
  //bool override_energy(double peakID, double energy);
  void clear_all_ROIs();

  std::set<double> get_selected_peaks() const;
  void set_selected_peaks(std::set<double> selected_peaks);

  //export results
  void save_report(std::string filename);

  friend void to_json(json& j, const Fitter &s);

  //for efficiency stuff
  std::string sample_name_;
  double activity_scale_factor_ {1.0}; //should be in spectrum?

  //data from spectrum
  ConsumerMetadata metadata_;
  Detector detector_; //need this? metadata?


private:
  std::map<double, ROI> regions_;
  std::set<double> selected_peaks_;
  FitEvaluation finder_;
  FitSettings settings_;

  void render_all();
  ROI *parent_of(double peakID);

  void filter_selection();

};

typedef std::shared_ptr<Fitter> FitterPtr;

}
