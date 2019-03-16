#pragma once

// CppNumericalSolver
#include <iostream>
#include <iomanip>
#include <core/util/eigen_fix.h>
#include "isolver.h"
#include "../linesearch/armijo.h"
#include "../linesearch/morethuente.h"
#include "../linesearch/brent.h"

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
  using Superclass::verbosity;
  using Superclass::os;

  void minimize(ProblemType& objFunc, TVector& x0) override
  {
    Brent<ProblemType, 1> brent;
    brent.verbosity = (verbosity > 1) ? (verbosity - 2) : 0;
    brent.os = os;

    const size_t DIM = x0.rows();
    THessian H = THessian::Identity(DIM, DIM);
    TVector grad(DIM);
    TVector x_old;
    this->m_current.reset();
    this->m_current.fDelta = std::numeric_limits<Scalar>::infinity();
    objFunc.gradient(x0, grad);

    if (verbosity >= 2)
    {
      (*os) << "init x=" << x0.transpose() << "\n";
      (*os) << "init g=" << grad.transpose() << "\n";
    }

    do
    {
      x_old = x0;

      TVector searchDir = -1 * H * grad;
      if (verbosity >= 2)
        (*os) << " searchDir=" << searchDir.transpose() << "\n";

      // check "positive definite"
      Scalar phi = grad.dot(searchDir);

      // positive definit ?
      if ((phi > 0) || (phi != phi))
      {
        // no, we reset the hessian approximation
        H = THessian::Identity(DIM, DIM);
        searchDir = -1 * grad;
        if (verbosity >= 2)
          (*os) << " resetting Hessian approximation after phi=" << phi
                << "   New searchDir=" << searchDir << "\n";
      }

//      const Scalar rate =
//          Armijo<ProblemType, 1>::linesearch(
//              x0, searchDir, objFunc, 1.0,
//              ((verbosity >= 2) ? os
//                                                         : nullptr));

      const Scalar rate = brent.linesearch(x0, searchDir, objFunc);

      if (verbosity >= 2)
      {
        (*os) << " rate=" << rate << "\n";
        (*os) << " rate*dir=" << (rate * searchDir).transpose() << "\n";
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

      if (verbosity >= 1)
      {
        (*os) << this->m_current;
        (*os) << "     f=" << std::right << std::setw(12) << objFunc.value(x0);
        if (verbosity >= 2)
        {
          (*os) << "\n  x=" << x0.transpose();
          (*os) << "\n  g=" << grad.transpose();
        }
        (*os) << "\n";
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
