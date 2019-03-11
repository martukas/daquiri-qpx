// CppNumericalSolver
#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>

namespace cppoptlib
{

/*template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
*/

enum class DebugLevel { None = 0, Low, High };
enum class Status
{
  NotStarted = -1,
  Continue = 0,
  IterationLimit,
  XDeltaTolerance,
  FDeltaTolerance,
  GradNormTolerance,
  Condition,
  UserDefined
};

enum class SimplexOp
{
  Place,
  Reflect,
  Expand,
  ContractIn,
  ContractOut,
  Shrink
};

inline std::ostream& operator<<(std::ostream& os, const SimplexOp& op)
{
  switch (op)
  {
    case SimplexOp::Place:
      os << "place";
      break;
    case SimplexOp::Reflect:
      os << "reflect";
      break;
    case SimplexOp::Expand:
      os << "expand";
      break;
    case SimplexOp::ContractIn:
      os << "contract-in";
      break;
    case SimplexOp::ContractOut:
      os << "contract-out";
      break;
    case SimplexOp::Shrink:
      os << "shrink";
      break;
  }
  return os;
}

inline std::string op_to_string(SimplexOp op)
{
  switch (op)
  {
    case SimplexOp::Place:
      return "place";
    case SimplexOp::Expand:
      return "expand";
    case SimplexOp::Reflect:
      return "reflect";
    case SimplexOp::ContractIn:
      return "contract-in";
    case SimplexOp::ContractOut:
      return "contract-out";
    case SimplexOp::Shrink:
      return "shrink";
  }
  return "unknown";
}

template<typename T>
class Criteria
{
 public:
  size_t iterations; //!< Maximum number of iterations
  T xDelta;          //!< Minimum change in parameter vector
  T fDelta;          //!< Minimum change in cost function
  T gradNorm;        //!< Minimum norm of gradient vector
  T condition;       //!< Maximum condition number of Hessian

  Criteria()
  {
    reset();
  }

  static Criteria defaults()
  {
    Criteria d;
    d.iterations = 10000;
    d.xDelta = 0;
    d.fDelta = 0;
    d.gradNorm = 1e-4;
    d.condition = 0;
    return d;
  }

  void reset()
  {
    iterations = 0;
    xDelta = 0;
    fDelta = 0;
    gradNorm = 0;
    condition = 0;
  }

  void print(std::ostream& os) const
  {
    os << "iter=" << std::right << std::setw(8) << iterations;
    os << "     \u0394x=" << std::right << std::setw(12) << xDelta;
    os << "     \u0394f=" << std::right << std::setw(12) << fDelta;
    os << "     ||g||=:   " << std::right << std::setw(12) << gradNorm;
    os << "     condition=" << std::right << std::setw(12) << condition;
  }
};

template<typename T>
Status checkConvergence(const Criteria<T>& stop, const Criteria<T>& current)
{

  if ((stop.iterations > 0) && (current.iterations > stop.iterations))
  {
    return Status::IterationLimit;
  }
  if ((stop.xDelta > 0) && (current.xDelta < stop.xDelta))
  {
    return Status::XDeltaTolerance;
  }
  if ((stop.fDelta > 0) && (std::abs(current.fDelta) < stop.fDelta))
  {
    return Status::FDeltaTolerance;
  }
  if ((stop.gradNorm > 0) && (current.gradNorm < stop.gradNorm))
  {
    return Status::GradNormTolerance;
  }
  if ((stop.condition > 0) && (current.condition > stop.condition))
  {
    return Status::Condition;
  }
  return Status::Continue;
}

inline std::ostream& operator<<(std::ostream& os, const Status& s)
{
  switch (s)
  {
    case Status::NotStarted:
      os << "Solver not started.";
      break;
    case Status::Continue:
      os << "Convergence criteria not reached.";
      break;
    case Status::IterationLimit:
      os << "Iteration limit reached.";
      break;
    case Status::XDeltaTolerance:
      os << "Change in parameter vector too small.";
      break;
    case Status::FDeltaTolerance:
      os << "Change in cost function value too small.";
      break;
    case Status::GradNormTolerance:
      os << "Gradient vector norm too small.";
      break;
    case Status::Condition:
      os << "Condition of Hessian/Covariance matrix too large.";
      break;
    case Status::UserDefined:
      os << "User-defined condition reached.";
      break;
  }
  return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Criteria<T>& c)
{
  c.print(os);
  return os;
}

} // End namespace cppoptlib

