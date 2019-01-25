#pragma once

#include <core/calibration/polynomial.h>
#include <core/fitting/sum4.h>

//#include <core/fitting/hypermet.h>
#include <core/fitting/hypermet/Hypermet.h>

#include <set>

namespace DAQuiri {

class Peak {
public:
  Peak(){}

  Peak(const nlohmann::json& j, const FCalibration& cal,
      const Finder &fs, const SUM4Edge& LB, const SUM4Edge& RB);

  Peak(const Hypermet &hyp, const SUM4 &s4, const FCalibration& cal);

  void reconstruct(const FCalibration& fs);

  //get rid of these
  std::vector<double> hr_peak_, hr_fullfit_;
  double intensity_theoretical_ {0.0};
  double efficiency_relative_ {0.0};

  const SUM4     &sum4() const { return sum4_;}
  const Hypermet &hypermet() const { return hypermet_;}

  const double &center() const {return center_;}
  const double &energy() const {return energy_;}
  const double &fwhm() const {return fwhm_;}
  const double &area_sum4() const {return area_sum4_;}
  const double &area_hyp() const {return area_hyp_;}
  const double &area_best() const {return area_best_;}
  const double &cps_sum4() const  {return cps_sum4_;}
  const double &cps_hyp() const {return cps_hyp_;}
  const double &cps_best() const  {return cps_best_;}

  void override_energy(double);

  int quality_energy() const;
  int quality_fwhm() const;
  bool good() const;

  bool operator<(const Peak &other) const;
  bool operator==(const Peak &other) const;
  bool operator!=(const Peak &other) const;
  bool operator>(const Peak &other) const;

  friend void to_json(nlohmann::json& j, const Peak &s);

  Hypermet hypermet_;

 private:
  SUM4 sum4_;
//  Hypermet hypermet_;

  // \todo use uncertain type
  double center_, energy_;
  double fwhm_;
  double area_sum4_, area_hyp_, area_best_;
  double cps_sum4_, cps_hyp_, cps_best_;
};

}
