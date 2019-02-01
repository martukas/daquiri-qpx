#pragma once

#include <core/fitting/weighted_data.h>
#include <core/fitting/uncertain.h>
#include <core/calibration/polynomial.h>

namespace DAQuiri {

class SUM4Edge {
  double Lchan_ {0};
  double Rchan_ {0};
  double min_ {0};
  double max_ {0};
  UncertainDouble dsum_ {0,0};
  UncertainDouble davg_ {0,0};

public:
  SUM4Edge() = default;
  SUM4Edge(const WeightedData& d);

  double left()  const {return Lchan_;}
  double right() const {return Rchan_;}
  double sum()      const {return dsum_.value();}
  double width()    const;
  double average()  const {return davg_.value();}
  double variance() const {return std::pow(davg_.sigma(), 2);}

  double min()      const {return min_;}
  double max()      const {return max_;}
  double midpoint() const;

  std::string to_string() const;

  friend void to_json(nlohmann::json& j, const SUM4Edge& s);
  friend void from_json(const nlohmann::json& j, SUM4Edge& s);

  static Polynomial sum4_background(const SUM4Edge& L,
                                    const SUM4Edge& R);
};

}
