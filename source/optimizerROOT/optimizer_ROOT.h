#pragma once

#include <core/fitting/optimizer.h>

namespace DAQuiri {

class OptimizerROOT : public Optimizer
{
public:

  void fit(std::shared_ptr<CoefFunction> func,
          const std::vector<double> &x,
          const std::vector<double> &y,
          const std::vector<double> &x_sigma,
          const std::vector<double> &y_sigma) override;

  void fit(Gaussian& gaussian,
           const std::vector<double> &x,
           const std::vector<double> &y) override;

  std::vector<Gaussian> fit_multiplet(const std::vector<double> &x,
                                      const std::vector<double> &y,
                                      std::vector<Gaussian> old,
                                      Polynomial &background,
                                      FitSettings settings) override;

  std::vector<Hypermet> fit_multiplet(const std::vector<double> &x,
                                      const std::vector<double> &y,
                                      std::vector<Hypermet> old,
                                      Polynomial &background,
                                      FitSettings settings) override;
};

}
