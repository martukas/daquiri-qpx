#pragma once

#include <core/calibration/polynomial.h>
#include <core/fitting/finder.h>
#include <cmath>

namespace DAQuiri {

class SUM4Edge {
  double Lchan_ {0};
  double Rchan_ {0};
  double min_ {0};
  double max_ {0};
  // \todo use uncertain types
  double dsum_, davg_;

public:
  SUM4Edge() = default;
  SUM4Edge(const nlohmann::json& j, const Finder& f);
  SUM4Edge(const std::vector<double> &x,
           const std::vector<double> &y,
           uint32_t left, uint32_t right);

  double left()  const {return Lchan_;}
  double right() const {return Rchan_;}
  double sum()      const {return dsum_;}
  double width()    const;
  double average()  const {return davg_;}
  double variance() const {
    return 0;
    // \todo return pow(davg_.uncertainty(),2);
  }

  double min()      const {return min_;}
  double max()      const {return max_;}
  double midpoint() const;

  friend void to_json(nlohmann::json& j, const SUM4Edge& s);
};

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

  // \todo use uncertain types
  double gross_area()      const {return gross_area_;}
  double background_area() const {return background_area_;}
  double peak_area()       const {return peak_area_;}
  double centroid()        const {return centroid_;}
  double fwhm()            const {return fwhm_;}

  static int get_currie_quality_indicator(double peak_net_area, double background_variance);

  static Polynomial sum4_background(const SUM4Edge& L,
                                    const SUM4Edge& R,
                                    const Finder& f);

  friend void to_json(nlohmann::json& j, const SUM4 &s);

private:
  SUM4Edge LB_, RB_;
  double Lchan_ {0};
  double Rchan_ {0};

  // \todo use uncertain types
  double gross_area_;
  double background_area_;
  double peak_area_;
  double centroid_;
  double fwhm_;

};

}
