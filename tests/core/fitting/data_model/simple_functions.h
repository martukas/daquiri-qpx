#include <core/fitting/data_model/data_model.h>
#include <core/fitting/data_model/value.h>

#include <core/util/UTF_extensions.h>

template<typename T>
class ConstFunction : public DAQuiri::FittableRegion
{
 public:
  T val;

  void update_indices() override
  {
    variable_count = 0;
    val.update_index(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    val.put(ret);
    return ret;
  }

  double eval(double chan) const override
  {
    (void) chan;
    return val.val();
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    (void) chan;
    return val.val_from(fit);
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    double ret = eval_at(chan, fit);
    if (val.valid_index())
      grads[val.index()] = val.grad_from(fit);
    return ret;
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    val.get(result.variables);

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();
    val.get_uncert(diags, chi_sq());
  }

  std::string to_string(std::string prepend = "") const override
  {
    Eigen::VectorXd grads;
    auto x = variables();
    chi_sq_gradient(x, grads);
    std::stringstream ss;
    ss << grads.transpose();
    return prepend + val.to_string()
        + "  chi" + UTF_superscript(2) + "=" + std::to_string(chi_sq())
        + "  grads=" + ss.str();
  }
};

template<typename T>
class LinearFunction : public DAQuiri::FittableRegion
{
 public:
  T val;

  void update_indices() override
  {
    variable_count = 0;
    val.update_index(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    val.put(ret);
    return ret;
  }

  double eval(double chan) const override
  {
    return val.val() * chan;
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return val.val_from(fit) * chan;
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    double ret = eval_at(chan, fit);
    if (val.valid_index())
      grads[val.index()] = val.grad_from(fit) * chan;
    return ret;
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    val.get(result.variables);

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();
    val.get_uncert(diags, chi_sq());
  }

  std::string to_string(std::string prepend = "") const override
  {
    Eigen::VectorXd grads;
    auto x = variables();
    chi_sq_gradient(x, grads);
    std::stringstream ss;
    ss << grads.transpose();
    return prepend + val.to_string()
        + "  chi" + UTF_superscript(2) + "=" + std::to_string(chi_sq())
        + "  grads=" + ss.str();
  }
};

template<typename T>
class QuadraticFunction : public DAQuiri::FittableRegion
{
 public:
  T val;

  void update_indices() override
  {
    variable_count = 0;
    val.update_index(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    val.put(ret);
    return ret;
  }

  double eval(double chan) const override
  {
    return val.val() * chan * chan;
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return val.val_from(fit) * chan * chan;
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    double ret = eval_at(chan, fit);
    if (val.valid_index())
      grads[val.index()] = val.grad_from(fit) * chan * chan;
    return ret;
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    val.get(result.variables);

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();
    val.get_uncert(diags, chi_sq());
  }

  std::string to_string(std::string prepend = "") const override
  {
    Eigen::VectorXd grads;
    auto x = variables();
    chi_sq_gradient(x, grads);
    std::stringstream ss;
    ss << grads.transpose();
    return prepend + val.to_string()
        + "  chi" + UTF_superscript(2) + "=" + std::to_string(chi_sq())
        + "  grads=" + ss.str();
  }
};
