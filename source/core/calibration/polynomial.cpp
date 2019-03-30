#include <core/calibration/polynomial.h>

#include <core/util/UTF_extensions.h>
#include <core/util/lexical_extensions.h>
#include <core/util/custom_logger.h>

namespace DAQuiri
{

Polynomial::Polynomial(const std::vector<double>& coeffs, double uncert)
{
  size_t i{0};
  for (const auto& c : coeffs)
    coeffs_[i++] = Parameter(c - uncert, c, c + uncert);
}

bool Polynomial::valid() const
{
  return !coeffs_.empty();
}

bool Polynomial::is_equal(CoefFunction* other) const
{
  if (this->type() != other->type())
    return false;
  return (coeffs_ == dynamic_cast<Polynomial*>(other)->coeffs_);
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
    ret += "p" + std::to_string(c.first);
    if (c.first > 0)
      ret += "*(x - x_offset)";
    if (c.first > 1)
      ret += "^" + std::to_string(c.first);
    i++;
    vars += "     p" + std::to_string(c.first) + "=" + c.second.to_string() + "\n";
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
    calib_eqn += to_str_precision(c.second.value(), precision);
    if (c.first > 0)
      calib_eqn += x_str;
    if (c.first > 1)
      calib_eqn += UTF_superscript(c.first);
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
    calib_eqn += to_str_precision(c.second.value(), precision);
    if (c.first > 0)
      calib_eqn += x_str;
    if (c.first > 1)
      calib_eqn += "<sup>" + std::to_string(c.first) + "</sup>";
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
  for (auto& c : coeffs_)
    result += c.second.value() * pow(x_adjusted, c.first);
  return result;
}

double Polynomial::derivative(double x) const
{
  Polynomial new_poly;  // derivative not true if offset != 0

  for (auto& c : coeffs_)
  {
    if (c.first != 0)
    {
      new_poly.coeffs_[c.first - 1]
          = {c.second.lower() * c.first,
             c.second.upper() * c.first,
             c.second.value() * c.first};
    }
  }

//  INFO("Poly deriv {} -> {}", to_UTF8(6, false), new_poly.to_UTF8(6, false));

  return new_poly(x);
}

}