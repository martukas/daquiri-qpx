#include <core/calibration/abstract_poly.h>

namespace DAQuiri
{

AbstractPolynomial::AbstractPolynomial(const std::vector<double>& coeffs)
{
  coefficients.resize(coeffs.size());
  for (size_t i = 0; i < coeffs.size(); ++i)
    coefficients[i].val(coeffs[i]);
}

bool AbstractPolynomial::valid() const
{
  return !coefficients.empty();
}

bool AbstractPolynomial::is_equal(CalibFunction* other) const
{
  if (this->type() != other->type())
    return false;
  auto* pother = dynamic_cast<AbstractPolynomial*>(other);
  if (coefficients.size() != pother->coefficients.size())
    return false;
  for (size_t i = 0; i < coefficients.size(); ++i)
    if (coefficients[i].x() != pother->coefficients[i].x())
      return false;
  return true;
}

void AbstractPolynomial::update_indices()
{
  variable_count = 0;
  for (auto& c : coefficients)
    c.update_index(variable_count);
}

Eigen::VectorXd AbstractPolynomial::variables() const
{
  Eigen::VectorXd ret;
  ret.setConstant(variable_count, 0.0);
  for (const auto& c : coefficients)
    c.put(ret);
  return ret;
}

void AbstractPolynomial::save_fit(const DAQuiri::FitResult& result)
{
  for (auto& c : coefficients)
    c.get(result.variables);

  if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
    return;

  double dof = degrees_of_freedom();

  Eigen::VectorXd diags = result.inv_hessian.diagonal();
  diags *= dof;

  double chisq_norm = this->chi_sq(result.variables) / dof;

  for (auto& c : coefficients)
    c.get_uncert(diags, chisq_norm);
}

nlohmann::json AbstractPolynomial::to_json() const
{
  nlohmann::json j;
  j["type"] = this->type();

  nlohmann::json coeffs;
  for (const auto&c : coefficients)
    coeffs.push_back(nlohmann::json(c));
  j["coefficients"] = coeffs;

  return j;
}

void AbstractPolynomial::from_json(const nlohmann::json& j)
{
  std::string t = j["type"];
  if (t != this->type())
    return;

  // \todo maybe throw instead?
  nlohmann::json coeffs = j["coefficients"];

  coefficients.clear();
  for (const nlohmann::json& c : coeffs)
    coefficients.emplace_back(UnboundedParam(c));
}

std::string AbstractPolynomial::to_string(std::string prepend) const
{
  std::string vars;
  size_t i {0};
  for (auto& c : coefficients)
    vars += prepend + " p" + std::to_string(i++) + "=" + c.to_string() + "\n";
  return vars;
}

}