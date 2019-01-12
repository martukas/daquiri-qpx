#include <core/fitting/gaussian.h>

namespace DAQuiri
{

void Gaussian::set_center(const Parameter& ncenter)
{
  center_ = ncenter;
//  chi2_ = 0;
}

void Gaussian::set_height(const Parameter& nheight)
{
  height_ = nheight;
//  chi2_ = 0;
}

void Gaussian::set_hwhm(const Parameter& nwidth)
{
  hwhm_ = nwidth;
//  chi2_ = 0;
}

void Gaussian::constrain_center(double min, double max)
{
  center_.constrain(min, max);
}

void Gaussian::constrain_height(double min, double max)
{
  height_.constrain(min, max);
}

void Gaussian::constrain_hwhm(double min, double max)
{
  hwhm_.constrain(min, max);
}

void Gaussian::set_chi2(double c2)
{
  chi2_ = c2;
}

std::string Gaussian::to_string() const
{
  std::string ret = "gaussian ";
  ret += "   area=" + std::to_string(area())
      + "   chi2=" + std::to_string(chi2_) + "    where:\n";

  ret += "     " + center_.to_string() + "\n";
  ret += "     " + height_.to_string() + "\n";
  ret += "     " + hwhm_.to_string() + "\n";

  return ret;
}

double Gaussian::evaluate(double x)
{
  return height_.value() *
      exp(-log(2.0) * (pow(((x - center_.value()) / hwhm_.value()), 2)));
}

double Gaussian::area() const
{
  return height_.value() * hwhm_.value() * sqrt(M_PI / log(2.0));
}

std::vector<double> Gaussian::evaluate_array(std::vector<double> x)
{
  std::vector<double> y;
  for (auto& q : x)
  {
    double val = evaluate(q);
    if (val < 0)
      y.push_back(0);
    else
      y.push_back(val);
  }
  return y;
}

}