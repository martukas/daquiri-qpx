#include <core/fitting/sum4/sum4edge.h>

namespace DAQuiri {

SUM4Edge::SUM4Edge(const nlohmann::json& j, const Finder& f)
  : SUM4Edge(f.x_, f.y_,
             f.find_index(j["left"]),
             f.find_index(j["right"]))
{}

SUM4Edge::SUM4Edge(const std::vector<double> &x,
                   const std::vector<double> &y,
                   uint32_t Lindex, uint32_t Rindex)
{
  dsum_ = 0.0;

  if (y.empty()
      || (y.size() != x.size())
      || (Lindex > Rindex)
      || (Lindex >= y.size())
      || (Rindex >= y.size()))
    return;

  Lchan_ = x.at(Lindex);
  Rchan_ = x.at(Rindex);

  min_ = std::numeric_limits<double>::max();
  max_ = std::numeric_limits<double>::min();

  for (size_t i=Lindex; i <= Rindex; ++i) {
    min_ = std::min(min_, y[i]);
    max_ = std::max(max_, y[i]);
    dsum_ += y[i]; // \todo  uncertaintye = sqrt(y[i]));
  }

  davg_ = dsum_ / width();
}

double SUM4Edge::width() const
{
  if (Rchan_ < Lchan_)
    return 0;
  else
    return (Rchan_ - Lchan_ + 1);
}

double SUM4Edge::midpoint() const
{
  double w = width();
  if (w <= 0)
    return std::numeric_limits<double>::quiet_NaN();
  else
    return Lchan_ + (w * 0.5);
}

void to_json(nlohmann::json& j, const SUM4Edge& s)
{
  j["left"] = s.Lchan_;
  j["right"] = s.Rchan_;
}


}
