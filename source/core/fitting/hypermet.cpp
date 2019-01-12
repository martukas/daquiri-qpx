#include <core/fitting/hypermet.h>

namespace DAQuiri
{

Hypermet::Hypermet(Gaussian gauss, FitSettings settings)
    : height_(gauss.height())
      , center_(gauss.center())
      , width_(gauss.hwhm().value() / sqrt(log(2)))
      , Lskew_amplitude_(settings.Lskew_amplitude), Lskew_slope_(settings.Lskew_slope)
      , Rskew_amplitude_(settings.Rskew_amplitude), Rskew_slope_(settings.Rskew_slope)
      , tail_amplitude_(settings.tail_amplitude), tail_slope_(settings.tail_slope)
      , step_amplitude_(settings.step_amplitude)
{
  if (settings.gaussian_only)
  {
    Lskew_enabled_ = false;
    Rskew_enabled_ = false;
    tail_enabled_ = false;
    step_enabled_ = false;
  }
}

void Hypermet::set_chi2(double c2)
{
  chi2_ = c2;
}

void Hypermet::set_center(const Parameter& ncenter)
{
  center_ = ncenter;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_height(const Parameter& nheight)
{
  height_ = nheight;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_width(const Parameter& nwidth)
{
  width_ = nwidth;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_Lskew_amplitude(const Parameter& nLskew_amplitude)
{
  Lskew_amplitude_ = nLskew_amplitude;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_Lskew_slope(const Parameter& nLskew_slope)
{
  Lskew_slope_ = nLskew_slope;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_Rskew_amplitude(const Parameter& nRskew_amplitude)
{
  Rskew_amplitude_ = nRskew_amplitude;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_Rskew_slope(const Parameter& nRskew_slope)
{
  Rskew_slope_ = nRskew_slope;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_tail_amplitude(const Parameter& ntail_amplitude)
{
  tail_amplitude_ = ntail_amplitude;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_tail_slope(const Parameter& ntail_slope)
{
  tail_slope_ = ntail_slope;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::set_step_amplitude(const Parameter& nstep_amplitude)
{
  step_amplitude_ = nstep_amplitude;
  user_modified_ = true;
  chi2_ = 0;
}

void Hypermet::constrain_center(double min, double max)
{
  center_.constrain(min, max);
  user_modified_ = true;
}

void Hypermet::constrain_height(double min, double max)
{
  height_.constrain(min, max);
  user_modified_ = true;
}

void Hypermet::constrain_width(double min, double max)
{
  width_.constrain(min, max);
  user_modified_ = true;
}

std::string Hypermet::to_string() const
{
  std::string ret = "Hypermet ";
  ret += "   area=" + std::to_string(area())
      + "   rsq=" + std::to_string(chi2_) + "    where:\n";

  ret += "     " + center_.to_string() + "\n";
  ret += "     " + height_.to_string() + "\n";
  ret += "     " + width_.to_string() + "\n";
  ret += "     " + Lskew_amplitude_.to_string() + "\n";
  ret += "     " + Lskew_slope_.to_string() + "\n";
  ret += "     " + Rskew_amplitude_.to_string() + "\n";
  ret += "     " + Rskew_slope_.to_string() + "\n";
  ret += "     " + tail_amplitude_.to_string() + "\n";
  ret += "     " + tail_slope_.to_string() + "\n";
  ret += "     " + step_amplitude_.to_string();

  return ret;
}

Gaussian Hypermet::gaussian() const
{
  Gaussian ret;
  ret.set_height(height_);
  ret.set_center(center_);
  auto w = ret.hwhm();
  w.set(width_.lower() * sqrt(log(2)),
        width_.upper() * sqrt(log(2)),
        width_.value() * sqrt(log(2)));
  ret.set_hwhm(w);
  ret.set_chi2(chi2_);
  return ret;
}

double Hypermet::eval_peak(double x) const
{
  if (width_.value() == 0)
    return 0;

  double xc = x - center_.value();

  double gaussian = exp(-pow(xc / width_.value(), 2));

  double left_short = 0;
  if (Lskew_enabled_)
  {
    double lexp =
        exp(pow(0.5 * width_.value() / Lskew_slope_.value(), 2) + xc / Lskew_slope_.value());
    if ((Lskew_slope_.value() != 0) && !std::isinf(lexp))
      left_short = Lskew_amplitude_.value() * lexp
          * erfc(0.5 * width_.value() / Lskew_slope_.value() + xc / width_.value());
  }

  double right_short = 0;
  if (Rskew_enabled_)
  {
    double rexp =
        exp(pow(0.5 * width_.value() / Rskew_slope_.value(), 2) - xc / Rskew_slope_.value());
    if ((Rskew_slope_.value() != 0) && !std::isinf(rexp))
      right_short = Rskew_amplitude_.value() * rexp
          * erfc(0.5 * width_.value() / Rskew_slope_.value() - xc / width_.value());
  }

  double ret = height_.value() * (gaussian + 0.5 * (left_short + right_short));

  return ret;
}

double Hypermet::eval_step_tail(double x) const
{
  if (width_.value() == 0)
    return 0;

  double xc = x - center_.value();

  double step = 0;
  if (step_enabled_)
    step = step_amplitude_.value() * erfc(xc / width_.value());

  double tail = 0;
  if (tail_enabled_)
  {
    double lexp =
        exp(pow(0.5 * width_.value() / tail_slope_.value(), 2) + xc / tail_slope_.value());
    if ((tail_slope_.value() != 0) && !std::isinf(lexp))
      tail = tail_amplitude_.value() * lexp
          * erfc(0.5 * width_.value() / tail_slope_.value() + xc / width_.value());
  }

  return height_.value() * 0.5 * (step + tail);
}

std::vector<double> Hypermet::peak(std::vector<double> x) const
{
  std::vector<double> y;
  for (auto& q : x)
    y.push_back(eval_peak(q));
  return y;
}

std::vector<double> Hypermet::step_tail(std::vector<double> x) const
{
  std::vector<double> y;
  for (auto& q : x)
    y.push_back(eval_step_tail(q));
  return y;
}

double Hypermet::area() const
{
  return height_.value() * width_.value() * sqrt(M_PI) *
      (1.0 +
          (Lskew_enabled_ ? Lskew_amplitude_.value() * width_.value() * Lskew_slope_.value()
                          : 0.0) +
          (Rskew_enabled_ ? Rskew_amplitude_.value() * width_.value() * Rskew_slope_.value()
                          : 0.0)
      );
}

bool Hypermet::gaussian_only() const
{
  return
      !(step_enabled_
          || tail_enabled_
          || Lskew_enabled_
          || Rskew_enabled_);
}

void to_json(nlohmann::json& j, const Hypermet& s)
{
  j["chi2"] = s.chi2();
  j["center"] = s.center();
  j["height"] = s.height();
  j["width"] = s.width();
  j["Lskew_amplitude"] = s.Lskew_amplitude();
  j["Lskew_slope"] = s.Lskew_slope();
  j["Rskew_amplitude"] = s.Rskew_amplitude();
  j["Rskew_slope"] = s.Rskew_slope();
  j["tail_amplitude"] = s.tail_amplitude();
  j["tail_slope"] = s.tail_slope();
  j["step_amplitude"] = s.step_amplitude();
  j["user_modified"] = s.user_modified();
}

void from_json(const nlohmann::json& j, Hypermet& s)
{
  s.center_ = j["center"];
  s.height_ = j["height"];
  s.width_ = j["width"];
  s.Lskew_amplitude_ = j["Lskew_amplitude"];
  s.Lskew_slope_ = j["Lskew_slope"];
  s.Rskew_amplitude_ = j["Rskew_amplitude"];
  s.Rskew_slope_ = j["Rskew_slope"];
  s.tail_amplitude_ = j["tail_amplitude"];
  s.tail_slope_ = j["tail_slope"];
  s.step_amplitude_ = j["step_amplitude"];
  s.user_modified_ = j["user_modified"];
  if (j.count("chi2") && j["chi2"].is_number_float())
    s.chi2_ = j["chi2"];
}

}
