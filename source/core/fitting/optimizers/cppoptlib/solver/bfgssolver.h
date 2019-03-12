#pragma once

// CppNumericalSolver
#include <iostream>
#include <iomanip>
#include <Eigen/LU>
#include "isolver.h"
#include "../linesearch/armijo.h"

namespace cppoptlib
{

template<typename ProblemType>
class BfgsSolver : public ISolver<ProblemType, 1>
{
 public:
  using Superclass = ISolver<ProblemType, 1>;
  using typename Superclass::Scalar;
  using typename Superclass::TVector;
  using typename Superclass::THessian;

  void minimize(ProblemType& objFunc, TVector& x0)
  {
    const size_t DIM = x0.rows();
    THessian H = THessian::Identity(DIM, DIM);
    TVector grad(DIM);
    TVector x_old;
    this->m_current.reset();
    this->m_current.fDelta = std::numeric_limits<Scalar>::infinity();
    objFunc.gradient(x0, grad);

    if (Superclass::m_debug >= DebugLevel::High) {
      std::cout << "init x=" << x0.transpose() << std::endl;
      std::cout << "init g=" << grad.transpose() << std::endl;
    }

    do
    {
      x_old = x0;

      TVector searchDir = -1 * H * grad;
      if (Superclass::m_debug >= DebugLevel::High)
        std::cout << "searchDir=" << searchDir.transpose() << std::endl;

      // check "positive definite"
      Scalar phi = grad.dot(searchDir);

      // positive definit ?
      if ((phi > 0) || (phi != phi))
      {
        // no, we reset the hessian approximation
        if (Superclass::m_debug >= DebugLevel::High)
          std::cout << "resetting Hessian approximation, phi=" << phi << std::endl;
        H = THessian::Identity(DIM, DIM);
        searchDir = -1 * grad;
      }

      const Scalar rate = Armijo<ProblemType, 1>::linesearch(x0, searchDir, objFunc);
      x0 = x0 + rate * searchDir;

      if (Superclass::m_debug >= DebugLevel::High)
        std::cout << "rate=" << rate << std::endl;

      TVector grad_old = grad;
      objFunc.gradient(x0, grad);
      TVector s = rate * searchDir;
      TVector y = grad - grad_old;

      const Scalar rho = 1.0 / y.dot(s);
      H = H - rho * (s * (y.transpose() * H) + (H * y) * s.transpose())
          + rho * (rho * y.dot(H * y) + 1.0) * (s * s.transpose());

      if (this->m_current.iterations)
        this->m_current.fDelta = objFunc.value(x0) - objFunc.value(x_old);
      this->m_current.xDelta = (x_old - x0).norm();
      this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
      this->m_status = checkConvergence(this->m_stop, this->m_current);

      if (Superclass::m_debug >= DebugLevel::Low)
      {
        std::cout << this->m_current;
        std::cout << "     f=" << std::right << std::setw(12) << objFunc.value(x0);
        if (Superclass::m_debug >= DebugLevel::High)
        {
          std::cout << "\n     x=" << x0.transpose();
          std::cout << "\n     g=" << grad.transpose();
        }
        std::cout << std::endl;
      }
      ++this->m_current.iterations;

      if (!objFunc.callback(this->m_current, x0))
        this->m_status = Status::UserDefined;

      if (!std::isfinite(H.diagonal().maxCoeff() / H.diagonal().minCoeff()))
        this->m_status = Status::Condition;
    }
    while (this->m_status == Status::Continue);
  }
};

}
/* namespace cppoptlib */
