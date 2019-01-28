#pragma once

#include <core/calibration/polynomial.h>
#include <core/fitting/sum4/sum4edge.h>

namespace DAQuiri {

class SUM4 {

public:
  SUM4() = default;
  SUM4(const nlohmann::json& j, const Finder& f,
       const SUM4Edge& LB, const SUM4Edge& RB);
  SUM4(double left, double right, const Finder& f,
       const SUM4Edge& LB, const SUM4Edge& RB);

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

  static int get_currie_quality_indicator(double peak_net_area, double background_variance);

  static Polynomial sum4_background(const SUM4Edge& L,
                                    const SUM4Edge& R,
                                    const Finder& f);

  friend void to_json(nlohmann::json& j, const SUM4 &s);

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
