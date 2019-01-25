#include <core/fitting/peak.h>

namespace DAQuiri {

Peak::Peak(const nlohmann::json& j, const Finder &f, const SUM4Edge &LB, const SUM4Edge &RB)
{
  if (j.count("hypermet"))
    hypermet_ = j["hypermet"];
  if (j.count("SUM4"))
    sum4_ = SUM4(j["SUM4"], f, LB, RB);
  reconstruct(f.settings_);
}

Peak::Peak(const Hypermet::Peak &hyp, const SUM4 &s4, const FitSettings &fs)
  : hypermet_(hyp)
  , sum4_(s4)
{
  reconstruct(fs);
}

void Peak::reconstruct(FitSettings fs)
{
  if (std::isfinite(hypermet_.amplitude.val()) && (hypermet_.amplitude.val() > 0)) {
    center_ = hypermet_.amplitude.val();
  }
  else
    center_ = sum4_.centroid();

//  if (!std::isfinite(center_.uncertainty())) {
////    DBG << "<Peak> overriding peak center uncert with sum4 for " << center_.to_string() << " with " << sum4_.centroid().to_string();
//    center_.setUncertainty(sum4_.centroid().uncertainty());
//  }


  double energyval = fs.cali_nrg_.transform(center_);
//  double emin = fs.cali_nrg_.transform(center_ - center_.uncertainty());
//  double emax = fs.cali_nrg_.transform(center_ + center_.uncertainty());
  energy_ = energyval; // \todo uncert = 0.5 * (emax - emin));
  //energy_.setSigFigs(center_.sigfigs());

  if (std::isfinite(hypermet_.width_.val()))
  {
    double L = hypermet_.position.val() - hypermet_.width_.val() * sqrt(log(2));
    double R = hypermet_.position.val() + hypermet_.width_.val() * sqrt(log(2));
//    double dmax = (hypermet_.width_.val() + hypermet_.width_.val().uncertainty()) * sqrt(log(2));
//    double dmin = (hypermet_.width_.val() - hypermet_.width_.val().uncertainty()) * sqrt(log(2));
//    double Lmax = hypermet_.position.val() - dmax;
//    double Rmax = hypermet_.position.val() + dmax;
//    double Lmin = hypermet_.position.val() - dmin;
//    double Rmin = hypermet_.position.val() + dmin;
    double val = fs.cali_nrg_.transform(R) - fs.cali_nrg_.transform(L);
//    double max = fs.cali_nrg_.transform(Rmax) - fs.cali_nrg_.transform(Lmax);
//    double min = fs.cali_nrg_.transform(Rmin) - fs.cali_nrg_.transform(Lmin);
    fwhm_ = val; // \todo uncert = (max - min);
  } else {
    double L = sum4_.centroid() - 0.5 * sum4_.fwhm();
    double R = sum4_.centroid() + 0.5 * sum4_.fwhm();
    fwhm_ = fs.cali_nrg_.transform(R) - fs.cali_nrg_.transform(L);
    // \todo uncert = std::numeric_limits<double>::quiet_NaN()
  }

  area_hyp_ =  hypermet_.area();
  area_sum4_ = sum4_.peak_area();

//  DBG << "<Peak> hyparea = " << area_hyp.debug();

  area_best_ = area_sum4_;

  cps_best_ = cps_hyp_ = cps_sum4_ = 0.0;

  double live_seconds = fs.live_time.count() * 0.001;

  if (live_seconds > 0) {
    cps_hyp_  = area_hyp_ / live_seconds;
    cps_sum4_ = area_sum4_ / live_seconds;
//    if (hypermet_.height_.value.finite() && (hypermet_.height_.val() > 0))
//      cps_best = cps_hyp;
//    else
      cps_best_ = cps_sum4_;
  }
}

bool Peak::good() const
{
  return ((sum4_.quality() == 1)
          && (quality_energy() == 1)
          && (quality_fwhm() == 1));
}

int Peak::quality_energy() const
{
//  if (energy_.error() > 50)
//    return 3;
//  else if (!std::isfinite(energy_.uncertainty()) || !energy_.uncertainty())
//    return 2;
//  else
    return 1;
}

int Peak::quality_fwhm() const
{
//  if (fwhm_.error() > 50)
//    return 3;
//  else if (!std::isfinite(fwhm_.uncertainty()) || !fwhm_.uncertainty())
//    return 2;
//  else
    return 1;
}

void Peak::override_energy(double newval)
{
  energy_ = newval;
}

bool Peak::operator<(const Peak &other) const
{
  return (center_ < other.center_);
}

bool Peak::operator==(const Peak &other) const
{
  return (center_ == other.center_);
}

bool Peak::operator!=(const Peak &other) const
{
  return (center_ != other.center_);
}

bool Peak::operator>(const Peak &other) const
{
  return (center_ > other.center_);
}

void to_json(nlohmann::json& j, const Peak &s)
{
  if (s.sum4_.peak_width())
    j["SUM4"] = s.sum4_;
  if (s.hypermet_.amplitude.val() > 0)
    j["hypermet"] = s.hypermet_;
}

}
