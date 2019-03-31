#include <core/calibration/coef_function.h>

namespace DAQuiri
{

std::vector<double> CoefFunction::eval_vector(const std::vector<double>& x) const
{
  std::vector<double> y;
  for (auto& q : x)
    y.push_back(this->eval(q));
  return y;
}

double CoefFunction::inverse(double y, double e) const
{
  int i = 0;
  double x0 = 0;
  double x1 = x0 + (y - this->eval(x0)) / (this->d_dx(x0));
  while (i <= 100 && std::abs(x1 - x0) > e)
  {
    x0 = x1;
    x1 = x0 + (y - this->eval(x0)) / (this->d_dx(x0));
    i++;
  }

  double x_adjusted = x1;

  if (std::abs(x1 - x0) <= e)
    return x_adjusted;

  else
  {
//    WARN("<" << this->type() << "> Maximum iteration reached in CoefFunction inverse evaluation";
    return nan("");
  }
}

void to_json(nlohmann::json& j, const CoefFunction& s)
{
  j["type"] = s.type();
//  for (auto c : s.coeffs())
//  {
//    nlohmann::json cc;
//    cc["degree"] = c.first;
//    cc["coefficient"] = c.second;
//    j["coefficients"].push_back(cc);
//  }
//  j["chi2"] = s.chi2();
}

void from_json(const nlohmann::json& j, CoefFunction& s)
{
//  if (j.count("coefficients"))
//  {
//    auto o = j["coefficients"];
//    for (nlohmann::json::iterator it = o.begin(); it != o.end(); ++it)
//    {
//      int d = it.value()["degree"];
//      Parameter p = it.value()["coefficient"];
//      s.set_coeff(d, p);
//    }
//  }
//  s.chi2(j["chi2"]);
}

}
