/**
 * @file sum4.cpp
 * @brief Construct for simple analysis of a background sample
 *
 * This is an abstraction of a peak area based on:
 * M. Lindstrom, Richard. (1994)
 * Sum and Mean Standard Programs for Activation Analysis.
 * Biological trace element research. 43-45. 597-603.
 * 10.1007/978-1-4757-6025-5_69.
 *
 * @author Martin Shetty
 */

#include <core/fitting/sum4/sum4.h>
#include <core/util/more_math.h>
#include <fmt/format.h>

namespace DAQuiri
{

SUM4::SUM4(const WeightedData& spectrum_data,
           const SUM4Edge& LB, const SUM4Edge& RB)
{
  if (spectrum_data.data.empty())
    throw std::runtime_error("Cannot create SUM4 with empty data");

  Polynomial background = SUM4Edge::sum4_background(LB, RB);
  double background_variance = square(0.5 * peak_width()) * (LB.variance() + RB.variance());

  Lchan_ = spectrum_data.data.front().x;
  Rchan_ = spectrum_data.data.back().x;

  for (const auto& p : spectrum_data.data)
    gross_area_ += {p.y, p.weight_true};

  background_area_ = {
      0.5 * peak_width() * (background(Rchan_) + background(Lchan_)),
      sqrt(background_variance)};

  peak_area_ = gross_area_ - background_area_;

  double sumYnet{0.0}, CsumYnet{0.0}, C2sumYnet{0.0};
  for (const auto& p : spectrum_data.data)
  {
    double yn = p.y - background(p.x);
    sumYnet += yn;
    CsumYnet += p.x * yn;
    C2sumYnet += square(p.x) * yn;
  }

  double centroidval = CsumYnet / sumYnet;
  //  if ((centroidval >= 0) && (centroidval < x.size()))
  //    centroidval = x.at(static_cast<size_t>(centroidval));
  double centroid_variance = (C2sumYnet / sumYnet) - square(centroidval);
  centroid_ = {centroidval, centroid_variance};

  double fwhm_val = 2.0 * sqrt(centroid_variance * log(4.0));
  // \todo true uncertainty?
  fwhm_ = {fwhm_val, std::numeric_limits<double>::quiet_NaN()};
}

UncertainDouble SUM4::peak_energy(const Calibration& cal) const
{
  return {cal.transform(centroid_.value()),
          cal.function()->derivative(centroid_.value()) * centroid_.sigma()};
}

UncertainDouble SUM4::fwhm_energy(const Calibration& cal) const
{
  double L = centroid_.value() - 0.5 * fwhm_.value();
  double R = centroid_.value() + 0.5 * fwhm_.value();
  return {cal.transform(R) - cal.transform(L), std::numeric_limits<double>::quiet_NaN()};
}

double SUM4::peak_width() const
{
  if (!std::isfinite(Rchan_) || !std::isfinite(Lchan_) || (Rchan_ < Lchan_))
    return 0;
  else
    return (Rchan_ - Lchan_ + 1);
}

int SUM4::quality() const
{
  return get_currie_quality_indicator(peak_area_.value(), square(background_area_.sigma()));
}

int SUM4::get_currie_quality_indicator(double peak_net_area, double background_variance)
{
  double currieLQ = 50 * (1 + sqrt(1 + background_variance / 12.5));
  double currieLD = 2.71 + 4.65 * sqrt(background_variance);
  double currieLC = 2.33 * sqrt(background_variance);

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

std::string SUM4::to_string() const
{
  return fmt::format("x=[{},{}] gross_area={} bkg_area={} peak_area={} centroid={} fwhm={}",
                     Lchan_, Rchan_,
                     gross_area_.to_string(),
                     background_area_.to_string(),
                     peak_area_.to_string(),
                     centroid_.to_string(),
                     fwhm_.to_string());
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
