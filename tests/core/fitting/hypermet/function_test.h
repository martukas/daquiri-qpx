#pragma once

#include "gtest_color_print.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Value.h>
#include <core/fitting/weighted_data.h>

#include <core/fitting/BFGS/Fittable.h>

class TestFittable : public DAQuiri::Fittable
{
 public:
  int32_t var_count{0};
  DAQuiri::WeightedData data;

  virtual double eval(double chan) const = 0;

  virtual double eval_at(double chan, const Eigen::VectorXd& fit) const = 0;

  virtual double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const = 0;

  DAQuiri::WeightedData generate_data(size_t bins) const
  {
    std::vector<double> channels;
    std::vector<double> y;
    for (size_t i = 0; i < bins; ++i)
    {
      channels.push_back(i);
      y.push_back(this->eval(i));
    }
    return DAQuiri::WeightedData(channels, y);
  }

  double chi_sq(const Eigen::VectorXd& fit) const override
  {
    double ChiSq = 0;
    double dof = data.data.size() - var_count;
    for (const auto& data : data.data)
    {
      double FTotal = this->eval_at(data.x, fit);
      ChiSq += square((data.y - FTotal) / data.weight_phillips_marlow);
    }
    return ChiSq / dof;
  }

  double operator ()(const Eigen::VectorXd& fit,
                     Eigen::VectorXd& gradients) const override
  {
    gradients.setConstant(fit.size(), 0.0);
    Eigen::VectorXd chan_gradients;
    chan_gradients.setConstant(fit.size(), 0.0);
    double Chisq{0.0};
    double dof = data.data.size() - var_count;
    for (const auto& data : data.data)
    {
      chan_gradients.setConstant(fit.size(), 0.0);

      double FTotal = this->eval_grad_at(data.x, fit, chan_gradients);

      double t3 = -2.0 * (data.y - FTotal) / square(data.weight_phillips_marlow);
      for (size_t var = 0; var < static_cast<size_t>(var_count); ++var)
        gradients[var] += chan_gradients[var] * t3;
      Chisq += square((data.y - FTotal) / data.weight_phillips_marlow);
    }
    return Chisq / dof;
  }
};


class FunctionTest : public TestBase
{
 protected:
  std::vector<double> val_proxy;
  std::vector<double> val_val;
  std::vector<double> chi_sq_norm;
  std::vector<double> gradient;

  void survey_grad(const TestFittable* fittable,
      DAQuiri::Value& variable,
      double step_size = 0.1)
  {
    size_t chosen_var_idx = variable.index();

    val_proxy.clear();
    val_val.clear();
    chi_sq_norm.clear();
    gradient.clear();

    Eigen::VectorXd variables = fittable->variables();
    Eigen::VectorXd gradients;
    for (double proxy = -M_PI; proxy < M_PI; proxy += step_size)
    {
      variables[chosen_var_idx] = proxy;

      val_proxy.push_back(proxy);
      val_val.push_back(variable.val_at(proxy));
      chi_sq_norm.push_back(fittable->operator()(variables, gradients));
      gradient.push_back(gradients[chosen_var_idx]);
    }
  }

  double check_chi_sq(bool print) const
  {
    auto min_chi = std::min_element(chi_sq_norm.begin(), chi_sq_norm.end());
    auto min_chi_i = std::distance(chi_sq_norm.begin(), min_chi);

    if (print)
    {
//      MESSAGE() << "chi_sq(proxy):\n" << visualize(val_proxy, chi_sq_norm, 100) << "\n";
      MESSAGE() << "chi_sq(val):\n" << visualize_all(val_val, chi_sq_norm, 100) << "\n";
    }
    MESSAGE() << "min(chi_sq)=" << chi_sq_norm[min_chi_i] << " at val=" << val_val[min_chi_i] << "\n";

    return val_val[min_chi_i];
  }

  double check_gradients(bool print) const
  {
    double min_abs = std::numeric_limits<double>::max();
    size_t grad_i = 0;
    for (size_t i=0; i < gradient.size(); ++i)
    {
      if (std::abs(gradient[i]) < min_abs)
      {
        grad_i = i;
        min_abs = std::abs(gradient[i]);
      }
    }

    if (print)
    {
//      MESSAGE() << "gradient(proxy):\n" << visualize(val_proxy, gradient, 100) << "\n";
      MESSAGE() << "gradient(val):\n" << visualize_all(val_val, gradient, 100) << "\n";
    }
    MESSAGE() << "min(abs(grad))=" << gradient[grad_i] << " at val=" << val_val[grad_i] << "\n";

    return val_val[grad_i];
  }


};
