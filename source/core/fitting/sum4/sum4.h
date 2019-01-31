#pragma once

#include <core/fitting/sum4/sum4edge.h>

namespace DAQuiri {

class SUM4 {

public:
  SUM4() = default;
  SUM4(const WeightedData& d, const SUM4Edge& LB, const SUM4Edge& RB);

  // \todo consider removing these
  SUM4Edge LB() const {return LB_;}
  SUM4Edge RB() const {return RB_;}

  double left()  const {return Lchan_;}
  double right() const {return Rchan_;}

  double peak_width() const;
  int    quality() const;

  UncertainDouble gross_area()      const {return gross_area_;}
  UncertainDouble background_area() const {return background_area_;}
  UncertainDouble peak_area()       const {return peak_area_;}
  UncertainDouble centroid()        const {return centroid_;}
  UncertainDouble fwhm()            const {return fwhm_;}

  UncertainDouble peak_energy(const Calibration& cal) const;
  UncertainDouble fwhm_energy(const Calibration& cal) const;

  static int get_currie_quality_indicator(double peak_net_area, double background_variance);

  friend void to_json(nlohmann::json& j, const SUM4& s);
  friend void from_json(const nlohmann::json& j, SUM4& s);

 private:
  SUM4Edge LB_, RB_;
  double Lchan_ {0};
  double Rchan_ {0};

  UncertainDouble gross_area_;
  UncertainDouble background_area_;
  UncertainDouble peak_area_;
  UncertainDouble centroid_;
  UncertainDouble fwhm_;
};

}
