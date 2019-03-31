#include <core/calibration/polynomial.h>

#include <core/util/UTF_extensions.h>
#include <core/util/lexical_extensions.h>
#include <core/util/custom_logger.h>

namespace DAQuiri
{

Polynomial::Polynomial(const std::vector<double>& coeffs, double uncert)
{
//  size_t i{0};
  coeffs_.resize(coeffs.size());
  for (size_t i=0; i < coeffs.size(); ++i)
    coeffs_[i].val(coeffs[i]);
//  for (const auto& c : coeffs)
//    coeffs_[i++].val(c);
  // = Parameter(c - uncert, c, c + uncert);
}

bool Polynomial::valid() const
{
  return !coeffs_.empty();
}

bool Polynomial::is_equal(CoefFunction* other) const
{
  if (this->type() != other->type())
    return false;
  auto* pother = dynamic_cast<Polynomial*>(other);
  if (coeffs_.size() != pother->coeffs_.size())
    return false;
  for (size_t i=0; i < coeffs_.size(); ++i)
    if (coeffs_[i].x() != pother->coeffs_[i].x())
      return false;
  return  true;
//  return (coeffs_ == dynamic_cast<Polynomial*>(other)->coeffs_);
}

std::string Polynomial::debug() const
{
  std::string ret = type() + " = ";
  std::string vars;
  int i = 0;
  for (auto& c : coeffs_)
  {
    if (i > 0)
      ret += " + ";
    ret += "p" + std::to_string(i);
    if (i > 0)
      ret += "*(x - x_offset)";
    if (i > 1)
      ret += "^" + std::to_string(i);
    i++;
    vars += "     p" + std::to_string(i) + "=" + c.to_string() + "\n";
  }

  return ret;
}

std::string Polynomial::to_UTF8(int precision, bool with_rsq) const
{
  std::string x_str = "x";

  std::string calib_eqn;
  int i = 0;
  for (auto& c : coeffs_)
  {
    if (i > 0)
      calib_eqn += " + ";
    calib_eqn += to_str_precision(c.val(), precision);
    if (i > 0)
      calib_eqn += x_str;
    if (i > 1)
      calib_eqn += UTF_superscript(i);
    i++;
  }

//  if (with_rsq)
//    calib_eqn += std::string("   r")
//        + UTF_superscript(2)
//        + std::string("=")
//        + to_str_precision(chi2_, precision);

  return calib_eqn;
}

std::string Polynomial::to_markup(int precision, bool with_rsq) const
{
  std::string x_str = "x";

  std::string calib_eqn;
  int i = 0;
  for (auto& c : coeffs_)
  {
    if (i > 0)
      calib_eqn += " + ";
    calib_eqn += to_str_precision(c.val(), precision);
    if (i > 0)
      calib_eqn += x_str;
    if (i > 1)
      calib_eqn += "<sup>" + std::to_string(i) + "</sup>";
    i++;
  }

//  if (with_rsq)
//    calib_eqn += "   r<sup>2</sup>"
//        + std::string("=")
//        + to_str_precision(chi2_, precision);

  return calib_eqn;
}

double Polynomial::operator()(double x) const
{
  double x_adjusted = x;
  double result = 0.0;
  int i = 0;
  for (auto& c : coeffs_)
  {
    result += c.val() * pow(x_adjusted, i);
    i++;
  }
  return result;
}

double Polynomial::derivative(double x) const
{
  Polynomial new_poly;  // derivative not true if offset != 0
  if (coeffs_.size() > 1)
  {
    new_poly.coeffs_.resize(coeffs_.size() - 1);
    for (size_t i=1; i < coeffs_.size(); ++i)
    {
        new_poly.coeffs_[i - 1].val(coeffs_[i].val() * i);
//          = {c.second.lower() * c.first,
//             c.second.upper() * c.first,
//             c.second.value() * c.first};
    }

  }


//  INFO("Poly deriv {} -> {}", to_UTF8(6, false), new_poly.to_UTF8(6, false));

  return new_poly(x);
}

}