#pragma once

#include "../meta.h"

#include <boost/math/tools/minima.hpp>
//#include <boost/math/tools/toms748_solve.hpp>

// For is_close_at_tolerance and is_small
//#include <boost/test/tools/floating_point_comparison.hpp>

#include <type_traits>
#include <typeinfo>

//! Test if two values are close within a given tolerance.
//template<typename FPT>
//inline bool
//is_close_to(FPT left, FPT right, FPT tolerance)
//{
//  return boost::math::fpc::close_at_tolerance<FPT>(tolerance)(left, right);
//}

//! Compare if value got is close to expected,
//! checking first if expected is very small
//! (to avoid divide by tiny or zero during comparison)
//! before comparing expect with value got.
//template<class T>
//bool is_close(T expect, T got, T tolerance)
//{
//  using boost::math::fpc::close_at_tolerance;
//  using boost::math::fpc::is_small;
//  using boost::math::fpc::FPC_STRONG;
//
//  if (is_small<T>(expect, tolerance))
//  {
//    return is_small<T>(got, tolerance);
//  }
//
//  return close_at_tolerance<T>(tolerance, FPC_STRONG)(expect, got);
//}

namespace cppoptlib
{

template<typename ProblemType, int Ord>
class Toms748
{
 public:
  using Scalar = typename ProblemType::Scalar;
  using TVector = typename ProblemType::TVector;

  std::ostream* os{nullptr};
  bool log_bracket{false};
  bool log_brent{false};
  size_t verbosity{0};

  size_t brent_maximum_iterations{500};
  Scalar brent_zeps{1e-10};

  Scalar bracket_glimit{100.0};
  Scalar bracket_tiny{1e-20};

  Scalar linmin_tolerance{0.0001};

  Scalar golden_ratio;

  Toms748()
  {
    golden_ratio = (1. + std::sqrt(5.)) / 2.;
  }

  struct StepEval
  {
    ProblemType* fittable_;
    const Eigen::VectorXd* variables_;
    const Eigen::VectorXd* search_direction_;

    Scalar size;
    Scalar f;
    Scalar dot;

    StepEval(ProblemType* fittable,
             const Eigen::VectorXd& variables,
             const Eigen::VectorXd& search_direction,
             Scalar lambda = 1.0)
        : fittable_(fittable)
          , variables_(&variables)
          , search_direction_(&search_direction)
          , size(lambda) {}

    StepEval(const StepEval& other) = default;

    StepEval& operator=(const StepEval& other)
    {
      fittable_ = other.fittable_;
      variables_ = other.variables_;
      search_direction_ = other.search_direction_;
      size = other.size;
      f = other.f;
      dot = other.dot;
      return *this;
    }

    void recalc_f(Scalar lambda)
    {
      size = lambda;
      f = fittable_->value((*variables_) + lambda * (*search_direction_));
    }

    void recalc_df(Scalar lambda)
    {
      recalc_f(lambda);
      Eigen::VectorXd gradient(variables_->size());
      fittable_->gradient((*variables_) + lambda * (*search_direction_), gradient);
      dot = gradient.dot(*search_direction_);
    }

    void recalc_f()
    {
      recalc_f(size);
    }

    void recalc_df()
    {
      recalc_df(size);
    }

    friend std::ostream& operator<<(std::ostream& stream, const StepEval& se)
    {
      stream << "(size:" << se.size << " f:" << se.f << " dot:" << se.dot << ")";
      return stream;
    }

    // for boost
    Scalar operator()(Scalar const& x)
    {
      return fittable_->value((*variables_) + x * (*search_direction_));
    }
  };

  inline static Scalar Sign(Scalar a, Scalar b)
  {
    if (b >= 0)
      return std::abs(a);
    else
      return -std::abs(a);
  }

  void bracket(StepEval& a, StepEval& b, StepEval& c)
  {
    // Just setting it up, value doesn't matter
    StepEval u = b;

    a.recalc_f();
    b.recalc_f();

    if (b.f > a.f)
      std::swap(a, b);

    c.recalc_f(b.size + 0.5 * (b.size - a.size));

    if (os && (verbosity > 1))
      (*os) << "     bracket starting a="
            << a << " b=" << b << " c=" << c << "\n";

    while (b.f > c.f)
    {
      Scalar r = (b.size - a.size) * (b.f - c.f);
      Scalar q = (b.size - c.size) * (b.f - a.f);
      Scalar n{true};
      u.size = std::abs(q - r);
      if (bracket_tiny > u.size)
        u.size = bracket_tiny;
      if (r > q)
        u.size = -u.size;
      u.size = b.size - ((b.size - c.size) * q
          - (b.size - a.size) * r) / (2 * u.size);
      Scalar ulim = b.size + bracket_glimit * (c.size - b.size);

      if (os && (verbosity > 1))
        (*os) << "      r=" << r << " q=" << q
              << " u.size=" << u.size << " ulim=" << ulim << "\n";

      if (os && (verbosity > 1))
        (*os) << "      ";

      if ((b.size - u.size) * (u.size - c.size) > 0)
      {
        u.recalc_f();

        if (u.f < c.f)
        {
          a = b;
          b = u;
          n = false;
          if (os && (verbosity > 2))
            (*os) << " case 1A ";
        }
        else if (u.f > b.f)
        {
          c = u;
          n = false;
          if (os && (verbosity > 2))
            (*os) << " case 1B ";
        }
        else
        {
          u.recalc_f(c.size + golden_ratio * (c.size - b.size));
          if (os && (verbosity > 2))
            (*os) << " case 1C ";
        }
      }
      else if ((c.size - u.size) * (u.size - ulim) > 0)
      {
        u.recalc_f();
        if (u.f < c.f)
        {
          b = c;
          c = u;
          u.recalc_f(c.size + golden_ratio * (c.size - b.size));
          if (os && (verbosity > 2))
            (*os) << " case 2A ";
        }
        else if (os && (verbosity > 2))
          (*os) << " case 2B ";

      }
      else if ((u.size - ulim) * (ulim - c.size) >= 0)
      {
        u.recalc_f(ulim);
        if (os && (verbosity > 2))
          (*os) << " case 3 ";
      }
      else
      {
        u.recalc_f(c.size + golden_ratio * (c.size - b.size));
        if (os && (verbosity > 2))
          (*os) << " case 4 ";
      }

      if (n)
      {
        a = b;
        b = c;
        c = u;
        if (os && (verbosity > 2))
          (*os) << "(N) ";
      }

      if (os && (verbosity > 1))
        (*os) << "a=" << a << " b=" << b << " c=" << c << " u=" << u << "\n";
    }
  }

  Scalar linesearch(const TVector& x, const TVector& searchDir, ProblemType& objFunc,
                    const Scalar alpha_init = 1.0)
  {
    StepEval step_min(&objFunc, x, searchDir, 0.0),
        step_init(&objFunc, x, searchDir, alpha_init),
        step_max(&objFunc, x, searchDir, 2.0 * alpha_init);

    if (os && (verbosity > 0))
      (*os) << "   line search\n"
            << "    min =" << step_min << "\n"
            << "    init=" << step_init << "\n"
            << "    max =" << step_max << "\n";

    bracket(step_min, step_max, step_init);

    if (step_min.size > step_max.size)
      std::swap(step_min, step_max);

    if (os && (verbosity > 0))
      (*os) << "   bracketed\n"
            << "    min =" << step_min << "\n"
            << "    init=" << step_init << "\n"
            << "    max =" << step_max << "\n";

    std::pair<Scalar, Scalar> ret{0., 0.};
    const boost::uintmax_t maxit = brent_maximum_iterations;
    boost::uintmax_t it = maxit;

    try
    { // Always use try'n'catch blocks with Boost.Math to ensure you get any error messages.

      int bits = std::numeric_limits<Scalar>::digits / 2; // Maximum is digits/2;

//      using tol_type = boost::math::tools::eps_tolerance<Scalar>;
//      tol_type tol(bits);
//      ret = boost::math::tools::toms748_solve<StepEval, Scalar, tol_type>(
//          step_init, step_min.size, step_max.size,
//          step_min.f, step_max.f, tol, it);

      ret = boost::math::tools::brent_find_minima<StepEval, Scalar>(
          step_init, step_min.size, step_max.size, bits, it);

    }
    catch (const std::exception& e)
    {
      // Always useful to include try & catch blocks because default policies
      // are to throw exceptions on arguments that cause errors like underflow, overflow.
      // Lacking try & catch blocks, the program will abort without a message below,
      // which may give some helpful clues as to the cause of the exception.

      if (os)
        (*os) << "Message from thrown boost exception was:\n   "
              << e.what() << "\n";
    }

    // x = ret.first;
    // f(x) = ret.second

    if (os && (verbosity > 0))
    {
      if (it >= maxit)
        (*os) << "   Brent did NOT meet after " << it << " iterations!\n";
      else
        (*os) << "   Brent met after " << it << " iterations, f("
              << ret.first << ")=" << ret.second << "\n";
    }

    return ret.first;
  }

};

}
