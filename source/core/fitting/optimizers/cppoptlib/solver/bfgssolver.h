#pragma once

// CppNumericalSolver
#include <iostream>
#include <iomanip>
#include <Eigen/LU>
#include "isolver.h"
#include "../linesearch/armijo.h"
#include "../linesearch/morethuente.h"

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

  void minimize(ProblemType& objFunc, TVector& x0, std::ostream* os = &std::cout) override
  {
    const size_t DIM = x0.rows();
    THessian H = THessian::Identity(DIM, DIM);
    TVector grad(DIM);
    TVector x_old;
    this->m_current.reset();
    this->m_current.fDelta = std::numeric_limits<Scalar>::infinity();
    objFunc.gradient(x0, grad);

    if (Superclass::m_debug >= DebugLevel::High)
    {
      (*os) << "init x=" << x0.transpose() << std::endl;
      (*os) << "init g=" << grad.transpose() << std::endl;
    }

    do
    {
      x_old = x0;

      TVector searchDir = -1 * H * grad;
      if (Superclass::m_debug >= DebugLevel::High)
        (*os) << "   searchDir=" << searchDir.transpose() << std::endl;

      // check "positive definite"
      Scalar phi = grad.dot(searchDir);

      // positive definit ?
      if ((phi > 0) || (phi != phi))
      {
        // no, we reset the hessian approximation
        H = THessian::Identity(DIM, DIM);
        searchDir = -1 * grad;
        if (Superclass::m_debug >= DebugLevel::High)
          (*os) << "   resetting Hessian approximation after phi=" << phi
                << "   New searchDir=" << searchDir << std::endl;
      }

      const Scalar rate =
          Armijo<ProblemType, 1>::linesearch(
              x0, searchDir, objFunc, 1.0,
              ((Superclass::m_debug >= DebugLevel::High) ? os
                                                         : nullptr));

      if (Superclass::m_debug >= DebugLevel::High)
      {
        (*os) << "   rate=" << rate << std::endl;
        (*os) << "   rate*dir=" << (rate * searchDir) << std::endl;
      }

      x0 = x0 + rate * searchDir;


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
        (*os) << this->m_current;
        (*os) << "     f=" << std::right << std::setw(12) << objFunc.value(x0);
        if (Superclass::m_debug >= DebugLevel::High)
        {
          (*os) << "\n     x=" << x0.transpose();
          (*os) << "\n     g=" << grad.transpose();
        }
        (*os) << std::endl;
      }
      ++this->m_current.iterations;

      if (!objFunc.callback(this->m_current, x0))
        this->m_status = Status::UserDefined;

//      if (!std::isfinite(searchDir))
//        this->m_status = Status::Condition;

      if (!std::isfinite(H.diagonal().maxCoeff() / H.diagonal().minCoeff()))
        this->m_status = Status::Condition;
    }
    while (this->m_status == Status::Continue);
  }
};

}
/* namespace cppoptlib */
