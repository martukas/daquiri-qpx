#include <core/calibration/polynomial.h>

#include <core/util/UTF_extensions.h>
#include <core/util/lexical_extensions.h>
#include <core/util/custom_logger.h>

namespace DAQuiri
{

Polynomial::Polynomial(const std::vector<double>& coeffs, double uncert)
{
  coeffs_.resize(coeffs.size());
  for (size_t i = 0; i < coeffs.size(); ++i)
    coeffs_[i].val(coeffs[i]);
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
  for (size_t i = 0; i < coeffs_.size(); ++i)
    if (coeffs_[i].x() != pother->coeffs_[i].x())
      return false;
  return true;
//  return (coeffs_ == dynamic_cast<Polynomial*>(other)->coeffs_);
}

void Polynomial::update_indices()
{
  variable_count = 0;
  for (auto& c : coeffs_)
    c.update_index(variable_count);
}

Eigen::VectorXd Polynomial::variables() const
{
  Eigen::VectorXd ret;
  ret.setConstant(variable_count, 0.0);
  for (const auto& c : coeffs_)
    c.put(ret);
  return ret;
}

double Polynomial::eval(double x) const
{
  double result{0.};
  int i = 0;
  for (auto& c : coeffs_)
  {
    result += c.val() * pow(x, i);
    i++;
  }
  return result;
}

double Polynomial::eval_grad_at(double chan, const Eigen::VectorXd& fit,
                                Eigen::VectorXd& grads) const
{
  double result{0.};
  int i = 0;
  for (auto& c : coeffs_)
  {
    result += c.val_from(fit) * pow(chan, i);
    grads[c.index()] += c.grad_from(fit) * std::pow(chan, i);
    i++;
  }
  return result;
}


void Polynomial::save_fit(const DAQuiri::FitResult& result)
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

double Polynomial::d_dx(double x) const
{
  Polynomial new_poly;  // derivative not true if offset != 0
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

  return new_poly.eval(x);
}


std::string Polynomial::to_string(std::string prepend) const
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
      ret += "*x";
    if (i > 1)
      ret += "^" + std::to_string(i);
    i++;
    vars += prepend + "  p" + std::to_string(i) + "=" + c.to_string() + "\n";
  }

  return ret + "\n" + vars;
}

std::string Polynomial::to_UTF8(int precision, bool with_rsq) const
{
  std::string calib_eqn;
  int i = 0;
  for (auto& c : coeffs_)
  {
    if (i > 0)
      calib_eqn += " + ";
    calib_eqn += to_str_precision(c.val(), precision);
    if (i > 0)
      calib_eqn += "x";
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
  std::string calib_eqn;
  int i = 0;
  for (auto& c : coeffs_)
  {
    if (i > 0)
      calib_eqn += " + ";
    calib_eqn += to_str_precision(c.val(), precision);
    if (i > 0)
      calib_eqn += "x";
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


}