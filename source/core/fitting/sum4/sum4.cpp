#include <core/fitting/sum4/sum4.h>

namespace DAQuiri
{

SUM4::SUM4(const nlohmann::json& j, const Finder& f, const SUM4Edge& LB, const SUM4Edge& RB)
    : SUM4(j["left"], j["right"], f, LB, RB) {}

Polynomial SUM4::sum4_background(const SUM4Edge& L, const SUM4Edge& R, const Finder& f)
{
  Polynomial sum4back;
  if (f.x_.empty())
    return sum4back;
  double run = R.left() - L.right();
  auto x_offset = sum4back.x_offset();
  x_offset.constrain(L.right(), L.right());
  double s4base = L.average();
  double s4slope = (R.average() - L.average()) / run;
  sum4back.x_offset(x_offset);
  sum4back.set_coeff(0, {s4base, s4base, s4base});
  sum4back.set_coeff(1, {s4slope, s4slope, s4slope});
  return sum4back;
}

SUM4::SUM4(double left, double right, const Finder& f,
           const SUM4Edge& LB, const SUM4Edge& RB)
{
  auto x = f.x_;
  auto y = f.y_;
  uint32_t Lindex = f.find_index(left);
  uint32_t Rindex = f.find_index(right);
  Polynomial background = sum4_background(LB, RB, f);

  if (y.empty()
      || (y.size() != x.size())
      || (Lindex > Rindex)
      || (Lindex >= y.size())
      || (Rindex >= y.size())
      || !LB.width()
      || !RB.width())
    return;

  LB_ = LB;
  RB_ = RB;
  Lchan_ = x[Lindex];
  Rchan_ = x[Rindex];

  gross_area_ = UncertainDouble::from_int(0, 0);
  for (size_t i = Lindex; i <= Rindex; ++i)
    gross_area_ += {y[i], sqrt(y[i])};

  double background_variance = pow((peak_width() / 2.0), 2) * (LB_.variance() + RB_.variance());
  background_area_ = {
      peak_width() * (background(x[Rindex]) + background(x[Lindex])) / 2.0,
      sqrt(background_variance)};

  peak_area_ = gross_area_ - background_area_;
  //peak_area_.autoSigs(1);

  double sumYnet(0), CsumYnet(0), C2sumYnet(0);
  for (size_t i = Lindex; i <= Rindex; ++i)
  {
    double yn = y[i] - background(y[i]);
    sumYnet += yn;
    CsumYnet += x[i] * yn;
    C2sumYnet += pow(x[i], 2) * yn;
  }

  double centroidval = CsumYnet / sumYnet;
  //  if ((centroidval >= 0) && (centroidval < x.size()))
  //    centroidval = x.at(static_cast<size_t>(centroidval));

  double centroid_variance = (C2sumYnet / sumYnet) - pow(centroidval, 2);
  centroid_ = {centroidval, centroid_variance};

  double fwhm_val = 2.0 * sqrt(centroid_variance * log(4));
  fwhm_ = {fwhm_val, std::numeric_limits<double>::quiet_NaN()};
}

double SUM4::peak_width() const
{
  if (Rchan_ < Lchan_)
    return 0;
  else
    return (Rchan_ - Lchan_ + 1);
}

int SUM4::quality() const
{
  return get_currie_quality_indicator(peak_area_.value(), pow(background_area_.uncertainty(), 2));
}

int SUM4::get_currie_quality_indicator(double peak_net_area, double background_variance)
{
  double currieLQ(0), currieLD(0), currieLC(0);
  currieLQ = 50 * (1 + sqrt(1 + background_variance / 12.5));
  currieLD = 2.71 + 4.65 * sqrt(background_variance);
  currieLC = 2.33 * sqrt(background_variance);

  if (peak_net_area > currieLQ)
    return 1;
  else if (peak_net_area > currieLD)
    return 2;
  else if (peak_net_area > currieLC)
    return 3;
  else if (peak_net_area > 0)
    return 4;
  else
    return 5;
}

void to_json(nlohmann::json& j, const SUM4& s)
{
  j["left"] = s.Lchan_;
  j["right"] = s.Rchan_;
}

}
