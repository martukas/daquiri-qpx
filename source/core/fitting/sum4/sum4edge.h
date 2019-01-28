#pragma once

#include <core/fitting/finder.h>
#include <core/fitting/uncertain.h>
#include <cmath>

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
  SUM4Edge(const nlohmann::json& j, const Finder& f);
  SUM4Edge(const std::vector<double> &x,
           const std::vector<double> &y,
           uint32_t left, uint32_t right);

  double left()  const {return Lchan_;}
  double right() const {return Rchan_;}
  double sum()      const {return dsum_.value();}
  double width()    const;
  double average()  const {return davg_.value();}
  double variance() const {return std::pow(davg_.uncertainty(), 2);}

  double min()      const {return min_;}
  double max()      const {return max_;}
  double midpoint() const;

  friend void to_json(nlohmann::json& j, const SUM4Edge& s);
};

}
