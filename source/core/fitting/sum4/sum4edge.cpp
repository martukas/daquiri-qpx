#include <core/fitting/sum4/sum4edge.h>
#include <core/util/more_math.h>
#include <fmt/format.h>

namespace DAQuiri
{

SUM4Edge::SUM4Edge(const WeightedData& d)
{
  if (d.empty())
    throw std::runtime_error("Cannot create SUM4Edge with empty data");

  dsum_ = {0.0, 0.0};

  if (d.data.empty())
    return;

  Lchan_ = d.data.front().x;
  Rchan_ = d.data.back().x;

  min_ = std::numeric_limits<double>::max();
  max_ = std::numeric_limits<double>::min();

  for (const auto& p : d.data)
  {
    min_ = std::min(min_, p.y);
    max_ = std::max(max_, p.y);
    dsum_ += {p.y, p.weight_true};
  }

  davg_ = dsum_ / width();
}

double SUM4Edge::width() const
{
  if (!std::isfinite(Rchan_) || !std::isfinite(Lchan_) || (Rchan_ < Lchan_))
    return 0;
  else
    return (Rchan_ - Lchan_ + 1);
}

double SUM4Edge::variance() const
{
  return square(davg_.sigma());
}

double SUM4Edge::midpoint() const
{
  double w = width();
  if (w <= 0)
    return std::numeric_limits<double>::quiet_NaN();
  else
    return Lchan_ + (w * 0.5);
}

std::string SUM4Edge::to_string() const
{
  return fmt::format("x=[{},{}] cts=[{},{}] sum={} avg={}",
      Lchan_, Rchan_, min_, max_, dsum_.to_string(), davg_.to_string());
}

void to_json(nlohmann::json& j, const SUM4Edge& s)
{
  j["left"] = s.Lchan_;
  j["right"] = s.Rchan_;
}

void from_json(const nlohmann::json& j, SUM4Edge& s)
{
  s.Lchan_ = j["left"];
  s.Rchan_ = j["right"];
}

Polynomial SUM4Edge::sum4_background(const SUM4Edge& L, const SUM4Edge& R)
{
  Polynomial sum4back;
  double run = R.left() - L.right();
  auto x_offset = sum4back.x_offset();
  x_offset.constrain(L.right(), L.right());
  double s4base = L.average().value();
  double s4slope = (R.average().value() - L.average().value()) / run;
  sum4back.x_offset(x_offset);
  sum4back.set_coeff(0, {s4base, s4base, s4base});
  sum4back.set_coeff(1, {s4slope, s4slope, s4slope});
  return sum4back;
}

}
