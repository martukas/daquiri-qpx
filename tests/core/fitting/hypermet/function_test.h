#pragma once

#include "gtest_color_print.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/fittable_region.h>
#include <core/fitting/hypermet/Value.h>

class FunctionTest : public TestBase
{
 protected:
  std::vector<double> val_proxy;
  std::vector<double> val_val;
  std::vector<double> chi_sq_norm;
  std::vector<double> gradient;

  DAQuiri::WeightedData generate_data(const DAQuiri::FittableRegion* fittable, size_t bins) const
  {
    std::vector<double> channels;
    std::vector<double> y;
    for (size_t i = 0; i < bins; ++i)
    {
      channels.push_back(i);
      y.push_back(fittable->eval(i));
    }
    return DAQuiri::WeightedData(channels, y);
  }

  void survey_grad(const DAQuiri::FittableRegion* fittable,
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
      chi_sq_norm.push_back(fittable->chi_sq_gradient(variables, gradients));
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
