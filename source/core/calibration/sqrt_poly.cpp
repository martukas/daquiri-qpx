#include <core/calibration/sqrt_poly.h>

#include <core/util/UTF_extensions.h>
#include <core/util/lexical_extensions.h>
#include <core/util/logger.h>

namespace DAQuiri
{

SqrtPoly::SqrtPoly(const std::vector<double>& coeffs)
{
  coeffs_.resize(coeffs.size());
  for (size_t i = 0; i < coeffs.size(); ++i)
    coeffs_[i].val(coeffs[i]);
}

bool SqrtPoly::valid() const
{
  return !coeffs_.empty();
}

bool SqrtPoly::is_equal(CalibFunction* other) const
{
  if (this->type() != other->type())
    return false;
  auto* pother = dynamic_cast<SqrtPoly*>(other);
  if (coeffs_.size() != pother->coeffs_.size())
    return false;
  for (size_t i = 0; i < coeffs_.size(); ++i)
    if (coeffs_[i].x() != pother->coeffs_[i].x())
      return false;
  return true;
//  return (coeffs_ == dynamic_cast<SqrtPoly*>(other)->coeffs_);
}

void SqrtPoly::update_indices()
{
  variable_count = 0;
  for (auto& c : coeffs_)
    c.update_index(variable_count);
}

Eigen::VectorXd SqrtPoly::variables() const
{
  Eigen::VectorXd ret;
  ret.setConstant(variable_count, 0.0);
  for (const auto& c : coeffs_)
    c.put(ret);
  return ret;
}

double SqrtPoly::eval(double x) const
{
  double result{0.};
  int i = 0;
  for (auto& c : coeffs_)
  {
    result += c.val() * pow(x, i);
    i++;
  }
  return std::sqrt(result);
}

double SqrtPoly::eval_grad_at(double x, const Eigen::VectorXd& fit,
                                Eigen::VectorXd& grads) const
{
  double result{0.};
  int i = 0;
  for (auto& c : coeffs_)
    result += c.val_from(fit) * pow(x, i++);

  i = 0;
  for (auto& c : coeffs_)
    grads[c.index()] += c.grad_from(fit) * std::pow(x, i++) / (2.0 * result);

  return result;
}

void SqrtPoly::save_fit(const DAQuiri::FitResult& result)
{
  for (auto& c : coeffs_)
    c.get(result.variables);

  if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
    return;

  double dof = degrees_of_freedom();

  Eigen::VectorXd diags = result.inv_hessian.diagonal();
  diags *= dof;

  double chisq_norm = this->chi_sq(result.variables) / dof;

  for (auto& c : coeffs_)
    c.get_uncert(diags, chisq_norm);
}

double SqrtPoly::d_dx(double x) const
{
  SqrtPoly new_poly;  // derivative not true if offset != 0
  if (coeffs_.size() > 1)
  {
    new_poly.coeffs_.resize(coeffs_.size() - 1);
    for (size_t i = 1; i < coeffs_.size(); ++i)
    {
      new_poly.coeffs_[i - 1].val(coeffs_[i].val() * i);
//          = {c.second.lower() * c.first,
//             c.second.upper() * c.first,
//             c.second.value() * c.first};
    }

  }

  return new_poly.eval(x) / (2.0 * eval(x));
}


std::string SqrtPoly::to_string(std::string prepend) const
{
  std::string ret = type() + " = sqrt(";
  std::string vars;
  int i = 0;
  for (auto& c : coeffs_)
  {
    if (i > 0)
      ret += "+";
    ret += "p" + std::to_string(i);
    if (i > 0)
      ret += "*x";
    if (i > 1)
      ret += "^" + std::to_string(i);
    i++;
    vars += prepend + "  p" + std::to_string(i) + "=" + c.to_string() + "\n";
  }

  return ret + ")\n" + vars;
}

std::string SqrtPoly::to_UTF8(int precision) const
{
  std::string ret = "\u221A(";
  int i = 0;
  for (auto& c : coeffs_)
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
  ret += ")";

  return ret;
}

std::string SqrtPoly::to_markup(int precision) const
{
  std::string ret = "&radic;<span style=\"text-decoration:overline;\">";
  int i = 0;
  for (auto& c : coeffs_)
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
  ret += "</span>";

  return ret;
}


}