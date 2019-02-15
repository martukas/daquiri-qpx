#pragma once

#include <core/fitting/abstract_optimizer.h>

#include <dlib/optimization.h>


namespace DAQuiri
{

// \todo add locks
class DLib : public AbstractOptimizer
{
 public:
  FitResult BFGSMin(Fittable* fittable, double tolf) override;

 private:
  Fittable* function_;

  using fitter_vector = dlib::matrix<double, 0, 1>;
  using fitter_matrix = dlib::matrix<double>;

  // dlib implementation
  double eval(const fitter_vector& m) const;
  fitter_vector derivative (const fitter_vector& m) const;
  fitter_matrix hessian(const fitter_vector& m) const;

};

}
