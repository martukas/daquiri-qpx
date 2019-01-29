#include <core/fitting/sum4/sum4.h>

namespace DAQuiri
{

SUM4::SUM4(const SpectrumData& d,
           const SUM4Edge& LB, const SUM4Edge& RB)
{
  Polynomial background = SUM4Edge::sum4_background(LB, RB);

  if (d.data.empty()
      || !LB.width()
      || !RB.width())
    return;

  LB_ = LB;
  RB_ = RB;
  Lchan_ = d.data.front().x;
  Rchan_ = d.data.back().x;

  gross_area_ = {0.0, 0.0};
  for (const auto& p : d.data)
    gross_area_ += {p.y, p.weight_true};

  double background_variance = pow(0.5 * peak_width(), 2) * (LB_.variance() + RB_.variance());
  background_area_ = {
      0.5 * peak_width() * (background(Rchan_) + background(Lchan_)),
      sqrt(background_variance)};

  peak_area_ = gross_area_ - background_area_;
  //peak_area_.autoSigs(1);

  double sumYnet(0), CsumYnet(0), C2sumYnet(0);
  for (const auto& p : d.data)
  {
    double yn = p.y - background(p.x);
    sumYnet += yn;
    CsumYnet += p.x * yn;
    C2sumYnet += pow(p.x, 2) * yn;
  }

  double centroidval = CsumYnet / sumYnet;
  //  if ((centroidval >= 0) && (centroidval < x.size()))
  //    centroidval = x.at(static_cast<size_t>(centroidval));

  double centroid_variance = (C2sumYnet / sumYnet) - pow(centroidval, 2);
  centroid_ = {centroidval, centroid_variance};

  double fwhm_val = 2.0 * sqrt(centroid_variance * log(4));
  fwhm_ = {fwhm_val, std::numeric_limits<double>::quiet_NaN()};
}

UncertainDouble SUM4::peak_energy(const Calibration& cal) const
{
  return {cal.transform(centroid_.value()),
          cal.function()->derivative(centroid_.value()) * centroid_.uncertainty()};
}

UncertainDouble SUM4::fwhm_energy(const Calibration& cal) const
{
  double L = centroid_.value() - 0.5 * fwhm_.value();
  double R = centroid_.value() + 0.5 * fwhm_.value();
  return {cal.transform(R) - cal.transform(L), std::numeric_limits<double>::quiet_NaN()};
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

void from_json(const nlohmann::json& j, SUM4& s)
{
  s.Lchan_ = j["left"];
  s.Rchan_ = j["right"];
}


}
