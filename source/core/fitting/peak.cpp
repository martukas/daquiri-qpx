#include <core/fitting/peak.h>

namespace DAQuiri
{

int value_quality(UncertainDouble ud, double error_threshold)
{
  if (ud.error() > error_threshold)
    return 3;
  else if (!std::isfinite(ud.uncertainty()) || !ud.uncertainty())
    return 2;
  return 1;
}


Peak::Peak(const Hypermet& hyp, const SUM4& s4, const FCalibration& cal)
    : hypermet_(hyp)
    , sum4_(s4)
{
  reconstruct(cal);
}

void Peak::reconstruct(const FCalibration& fs)
{
  if (std::isfinite(hypermet_.amplitude.val()) && (hypermet_.amplitude.val() > 0)) {
    center_ = UncertainDouble(hypermet_.amplitude.val(), hypermet_.amplitude.uncert_value);
  }
  else
    center_ = sum4_.centroid();

  if (!std::isfinite(center_.uncertainty())) {
//    DBG << "<Peak> overriding peak center uncert with sum4 for " << center_.to_string() << " with " << sum4_.centroid().to_string();
    center_.setUncertainty(sum4_.centroid().uncertainty());
  }


  double energyval = fs.bin_to_nrg(center_.value());
  double emin = fs.cali_nrg_.transform(center_.value() - center_.uncertainty());
  double emax = fs.cali_nrg_.transform(center_.value() + center_.uncertainty());
  energy_ = UncertainDouble::from_double(energyval, 0.5 * (emax - emin));

  if (std::isfinite(hypermet_.width_.val()))
  {
    double L = hypermet_.position.val() - hypermet_.width_.val() * sqrt(log(2.0));
    double R = hypermet_.position.val() + hypermet_.width_.val() * sqrt(log(2.0));
    double dmax = (hypermet_.width_.val() + hypermet_.width_.uncert_value) * sqrt(log(2));
    double dmin = (hypermet_.width_.val() - hypermet_.width_.uncert_value) * sqrt(log(2));
    double Lmax = hypermet_.position.val() - dmax;
    double Rmax = hypermet_.position.val() + dmax;
    double Lmin = hypermet_.position.val() - dmin;
    double Rmin = hypermet_.position.val() + dmin;
    double val = fs.cali_nrg_.transform(R) - fs.cali_nrg_.transform(L);
    double max = fs.cali_nrg_.transform(Rmax) - fs.cali_nrg_.transform(Lmax);
    double min = fs.cali_nrg_.transform(Rmin) - fs.cali_nrg_.transform(Lmin);
    double uncert = (max - min);
    fwhm_ = UncertainDouble::from_double(val, uncert);
  }
  else
  {
    double L = sum4_.centroid().value() - 0.5 * sum4_.fwhm().value();
    double R = sum4_.centroid().value() + 0.5 * sum4_.fwhm().value();
    fwhm_ = UncertainDouble::from_double(
        fs.cali_nrg_.transform(R) - fs.cali_nrg_.transform(L),
        std::numeric_limits<double>::quiet_NaN());
  }

  // \todo should be precalculated for hypermet!
  area_hyp_ = UncertainDouble::from_double(hypermet_.area(), hypermet_.area_uncert(1.0));
  area_sum4_ = sum4_.peak_area();

//  DBG << "<Peak> hyparea = " << area_hyp.debug();

  area_best_ = area_sum4_;

  cps_best_ = cps_hyp_ = cps_sum4_ = UncertainDouble::from_double(0,0);

//  double live_seconds = fs.live_time.count() * 0.001;
//
//  if (live_seconds > 0)
//  {
//    cps_hyp_ = area_hyp_ / live_seconds;
//    cps_sum4_ = area_sum4_ / live_seconds;
//    if (std::isfinite(hypermet_.amplitude.val()) && (hypermet_.amplitude.val() > 0))
//      cps_best_ = cps_hyp_;
//    else
//      cps_best_ = cps_sum4_;
//  }
}

bool Peak::good() const
{
  return ((sum4_.quality() == 1)
      && (quality_energy() == 1)
      && (quality_fwhm() == 1));
}

int Peak::quality_energy() const
{
  return value_quality(energy_, 50);
}

int Peak::quality_fwhm() const
{
  return value_quality(fwhm_, 50);
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

Peak::Peak(const nlohmann::json& j, const FCalibration& cal,
           const Finder& f, const SUM4Edge& LB, const SUM4Edge& RB)
{
  if (j.count("hypermet"))
    hypermet_ = j["hypermet"];
  if (j.count("SUM4"))
    sum4_ = SUM4(j["SUM4"], f, LB, RB);
  reconstruct(cal);
}


}
