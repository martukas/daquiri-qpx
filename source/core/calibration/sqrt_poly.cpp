#include <core/calibration/sqrt_poly.h>

#include <core/util/UTF_extensions.h>
#include <core/util/lexical_extensions.h>
#include <core/util/logger.h>

namespace DAQuiri
{

double SqrtPoly::eval(double x) const
{
  return std::sqrt(Polynomial::eval(x));
}

double SqrtPoly::eval_at(double x, const Eigen::VectorXd& fit) const
{
  return std::sqrt(Polynomial::eval_at(x, fit));
}

double SqrtPoly::eval_grad_at(double x, const Eigen::VectorXd& fit,
                              Eigen::VectorXd& grads) const
{
  auto ret = SqrtPoly::eval_at(x, fit);

  size_t i{0};
  for (auto& c : coefficients)
    grads[c.index()] += c.grad_from(fit) * std::pow(x, i++) / (2.0 * ret);

  return ret;
}

double SqrtPoly::d_dx(double x) const
{
  return Polynomial::d_dx(x) / (2.0 * eval(x));
}

std::string SqrtPoly::to_string(std::string prepend) const
{
  std::string ret = type() + " = sqrt(";
  int i = 0;
  for (auto& c : coefficients)
  {
    if (i > 0)
      ret += "+";
    ret += "p" + std::to_string(i);
    if (i > 0)
      ret += "*x";
    if (i > 1)
      ret += "^" + std::to_string(i);
    i++;
  }

  return ret + ")\n" + AbstractPolynomial::to_string(prepend);
}

std::string SqrtPoly::to_UTF8(int precision) const
{
  return "\u221A(" + Polynomial::to_UTF8(precision) + ")";
}

std::string SqrtPoly::to_markup(int precision) const
{
  return "&radic;<span style=\"text-decoration:overline;\">"
      + Polynomial::to_markup(precision) + "</span>";
}

}