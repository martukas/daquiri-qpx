#pragma once

#include <core/fitting/sum4/sum4.h>
#include <core/fitting/hypermet/Hypermet.h>

namespace DAQuiri
{

class Peak
{
 public:
  Peak() = default;

  Peak(const Hypermet& hyp, const SUM4& s4, const FCalibration& cal);

  void reconstruct(const FCalibration& fs);

  // \todo get rid of these
  double intensity_theoretical_{0.0};
  double efficiency_relative_{0.0};

  const SUM4& sum4() const { return sum4_; }
  const Hypermet& hypermet() const { return hypermet_; }

  const UncertainDouble& center() const { return center_; }
  const UncertainDouble& energy() const { return energy_; }
  const UncertainDouble& fwhm() const { return fwhm_; }
  const UncertainDouble& area_sum4() const { return area_sum4_; }
  const UncertainDouble& area_hyp() const { return area_hyp_; }
  const UncertainDouble& area_best() const { return area_best_; }
  const UncertainDouble& cps_sum4() const { return cps_sum4_; }
  const UncertainDouble& cps_hyp() const { return cps_hyp_; }
  const UncertainDouble& cps_best() const { return cps_best_; }

  void override_energy(double);

  bool operator<(const Peak& other) const;
  bool operator==(const Peak& other) const;
  bool operator!=(const Peak& other) const;
  bool operator>(const Peak& other) const;

  friend void to_json(nlohmann::json& j, const Peak& s);
  Peak(const nlohmann::json& j, const FCalibration& cal,
       const Finder& fs, const SUM4Edge& LB, const SUM4Edge& RB);

 private:
  SUM4 sum4_;
  Hypermet hypermet_;

  UncertainDouble center_, energy_;
  UncertainDouble fwhm_;
  UncertainDouble area_sum4_, area_hyp_, area_best_;
  UncertainDouble cps_sum4_, cps_hyp_, cps_best_;
};


}
