#pragma once

#include <core/calibration/coef_function.h>
#include <core/fitting/gaussian.h>
#include <core/fitting/hypermet.h>
#include <memory>

namespace DAQuiri
{

class Optimizer
{
 public:

  virtual void fit(std::shared_ptr<CoefFunction> func,
                   const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& x_sigma,
                   const std::vector<double>& y_sigma) = 0;

  virtual void fit(Gaussian& gaussian,
                   const std::vector<double>& x,
                   const std::vector<double>& y) = 0;

  virtual std::vector<Gaussian> fit_multiplet(const std::vector<double>& x,
                                              const std::vector<double>& y,
                                              std::vector<Gaussian> old,
                                              Polynomial& background,
                                              FitSettings settings) = 0;

  virtual std::vector<Hypermet> fit_multiplet(const std::vector<double>& x,
                                              const std::vector<double>& y,
                                              std::vector<Hypermet> old,
                                              Polynomial& background,
                                              FitSettings settings) = 0;

  static void sanity_check(Gaussian& gaussian,
                           double xmin, double xmax,
                           double ymin, double ymax);

  static void constrain_center(Gaussian& gaussian, double slack);

  static void sanity_check(Hypermet& hyp,
                           double xmin, double xmax,
                           double ymin, double ymax);

  static void constrain_center(Hypermet& gaussian, double slack);

  static void initial_sanity(Gaussian& gaussian,
                             double xmin, double xmax,
                             double ymin, double ymax);

};

using OptimizerPtr = std::shared_ptr<Optimizer>;

class OptimizerFactory
{
 public:
  static OptimizerFactory& getInstance()
  {
    static OptimizerFactory singleton_instance;
    return singleton_instance;
  }

  OptimizerPtr create_any() const;
  OptimizerPtr create_type(std::string type) const;
  void register_type(std::string name, std::function<Optimizer*(void)> typeConstructor);
  std::vector<std::string> types() const;

 private:
  std::map<std::string, std::function<Optimizer*(void)>> constructors;

  //singleton assurance
  OptimizerFactory() {}
  OptimizerFactory(OptimizerFactory const&);
  void operator=(OptimizerFactory const&);
};

template<class T>
class OptimizerRegistrar
{
 public:
  OptimizerRegistrar(std::string name)
  {
    OptimizerFactory::getInstance().register_type(name,
                                                  [](void) -> Optimizer* { return new T(); });
  }
};

}
