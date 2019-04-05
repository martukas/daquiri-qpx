#include <core/calibration/polynomial.h>

#include <core/util/UTF_extensions.h>
#include <core/util/lexical_extensions.h>
#include <core/util/logger.h>

namespace DAQuiri
{

double Polynomial::eval(double x) const
{
  double result{0.};
  size_t i{0};
  for (auto& c : coefficients)
    result += c.val() * pow(x, i++);
  return result;
}

double Polynomial::eval_at(double x, const Eigen::VectorXd& fit) const
{
  double result{0.};
  size_t i{0};
  for (auto& c : coefficients)
    result += c.val_from(fit) * pow(x, i++);
  return result;
}

double Polynomial::eval_grad_at(double x, const Eigen::VectorXd& fit,
                                Eigen::VectorXd& grads) const
{
  double result{0.};
  size_t i {0};
  for (auto& c : coefficients)
  {
    double p = pow(x, i++);
    result += c.val_from(fit) * p;
    grads[c.index()] += c.grad_from(fit) * p;
  }
  return result;
}

double Polynomial::d_dx(double x) const
{
  Polynomial new_poly;
  if (coefficients.size() > 1)
  {
    new_poly.coefficients.resize(coefficients.size() - 1);
    for (size_t i = 1; i < coefficients.size(); ++i)
      new_poly.coefficients[i - 1].val(coefficients[i].val() * i);
  }
  return new_poly.eval(x);
}

std::string Polynomial::to_string(std::string prepend) const
{
  std::string ret = type() + " = ";
  int i = 0;
  for (auto& c : coefficients)
  {
    if (i > 0)
      ret += " + ";
    ret += "p" + std::to_string(i);
    if (i > 0)
      ret += "*x";
    if (i > 1)
      ret += "^" + std::to_string(i);
    i++;
  }

  return ret + "\n" + AbstractPolynomial::to_string(prepend);
}

std::string Polynomial::to_UTF8(int precision) const
{
  std::string ret;
  int i = 0;
  for (auto& c : coefficients)
  {
    if (i > 0)
      ret += "+";
    ret += to_str_precision(c.val(), precision);
    if (i > 0)
      ret += "x";
    if (i > 1)
      ret += UTF_superscript(i);
    i++;
  }
  return ret;
}

std::string Polynomial::to_markup(int precision) const
{
  std::string ret;
  int i = 0;
  for (auto& c : coefficients)
  {
    if (i > 0)
      ret += "+";
    ret += to_str_precision(c.val(), precision);
    if (i > 0)
      ret += "x";
    if (i > 1)
      ret += "<sup>" + std::to_string(i) + "</sup>";
    i++;
  }
  return ret;
}

}