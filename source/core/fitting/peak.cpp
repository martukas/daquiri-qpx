#include <core/fitting/peak.h>

namespace DAQuiri
{

Peak::Peak(const Hypermet& hyp, const SUM4& s4, const FCalibration& cal)
    : hypermet_(hyp)
      , sum4_(s4)
{
  reconstruct(cal);
}

void Peak::reconstruct(const FCalibration& fs)
{
  if (std::isfinite(hypermet_.amplitude.val()) && (hypermet_.amplitude.val() > 0))
    center_ = hypermet_.peak_position();
  else
    center_ = sum4_.centroid();

  if (!std::isfinite(center_.uncertainty()))
  {
//    DBG << "<Peak> overriding peak center uncert with sum4 for " << center_.to_string() << " with " << sum4_.centroid().to_string();
    center_.setUncertainty(sum4_.centroid().uncertainty());
  }

  energy_ = hypermet_.peak_energy(fs.cali_nrg_);

  if (std::isfinite(hypermet_.width_.val()))
    fwhm_ = hypermet_.fwhm_energy(fs.cali_nrg_);
  else
    fwhm_ = sum4_.fwhm_energy(fs.cali_nrg_);

  area_hyp_ = hypermet_.area();
  area_sum4_ = sum4_.peak_area();

//  DBG << "<Peak> hyparea = " << area_hyp.debug();

  area_best_ = area_sum4_;

  cps_best_ = cps_hyp_ = cps_sum4_ = {0.0, 0.0};
}

void Peak::override_energy(double newval)
{
  energy_.setValue(newval);
}

bool Peak::operator<(const Peak& other) const
{
  return (center_ < other.center_);
}

bool Peak::operator==(const Peak& other) const
{
  return (center_ == other.center_);
}

bool Peak::operator!=(const Peak& other) const
{
  return !(center_ == other.center_);
}

bool Peak::operator>(const Peak& other) const
{
  return (center_ > other.center_);
}

void to_json(nlohmann::json& j, const Peak& s)
{
  if (s.sum4_.peak_width())
    j["SUM4"] = s.sum4_;
  if (s.hypermet_.amplitude.val() > 0)
    j["hypermet"] = s.hypermet_;
}

Peak::Peak(
    const nlohmann::json& j,
    const FCalibration& cal,
    const Finder& f,
    const SUM4Edge& LB,
    const SUM4Edge& RB)
{
  if (j.count("hypermet"))
    hypermet_ = j["hypermet"];
  if (j.count("SUM4"))
    sum4_ = SUM4(j["SUM4"], f, LB, RB);
  reconstruct(cal);
}

}
