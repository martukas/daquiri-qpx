#pragma once

#include "../meta.h"

namespace cppoptlib
{

template<typename ProblemType, int Ord>
class Brent
{
 public:
  using Scalar = typename ProblemType::Scalar;
  using TVector = typename ProblemType::TVector;

  std::ostream* os {nullptr};
  bool log_bracket {false};
  bool log_brent {false};
  size_t verbosity {0};

  size_t brent_maximum_iterations{500};
  Scalar brent_zeps{1e-10};

  Scalar bracket_glimit{100.0};
  Scalar bracket_tiny{1e-20};

  Scalar linmin_tolerance{0.0001};

  Scalar golden_ratio;

  Brent()
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

    friend std::ostream& operator<<(std::ostream& stream,
                                    const StepEval& se)
    {
      stream << "(size:" << se.size << " f:" << se.f << " dot:" << se.dot << ")";
      return stream;
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

    c.recalc_f(b.size + golden_ratio * (b.size - a.size));

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
        (*os) << "a=" << a << " b=" << b << " c=" << c << " u=" << u  << "\n";
    }
  }

  StepEval brent_search(StepEval x, Scalar lambda1, Scalar lambda2)
  {
    x.recalc_df();
    StepEval v = x;
    StepEval w = x;

    StepEval u = x;
    u.recalc_df(0.);

    Scalar lambda_min = std::min(lambda1, lambda2);
    Scalar lambda_max = std::max(lambda1, lambda2);

    bool done{false};
    Scalar e{0.0};
    Scalar d{0.0}; // \todo is this really the default value?

    if (os && (verbosity > 1))
      (*os) << "     Brent starting x=" << x
            << "  u=" << u
            << "  lambda_min=" << lambda_min
            << "  lambda_max=" << lambda_max << "\n";

    // \todo check for cancel
    for (size_t iteration = 0; iteration < brent_maximum_iterations; ++iteration)
    {
      Scalar lambda_mid = 0.5 * (lambda_min + lambda_max);
      Scalar tol1 = linmin_tolerance * std::abs(x.size) + brent_zeps;
      Scalar tol2 = 2. * tol1;

      done = (std::abs(x.size - lambda_mid) <= (tol2 - 0.5 * (lambda_max - lambda_min)));

      if (os && (verbosity > 1))
        (*os) << "     i=" << iteration
              << " lambda_mid=" << lambda_mid
              << " tol1=" << tol1
              << " tol2=" << tol2
              << " done=" << done
              << "\n";

      if (!done)
      {
        bool ok1 = false;
        if (std::abs(e) > tol1)
        {
          Scalar d1 = 2 * (lambda_max - lambda_min);
          Scalar d2 = d1;
          if (w.dot != x.dot)
            d1 = (w.size - x.size) * x.dot / (x.dot - w.dot);
          if (v.dot != x.dot)
            d2 = (v.size - x.size) * x.dot / (x.dot - v.dot);
          Scalar u1 = x.size + d1;
          Scalar u2 = x.size + d2;
          ok1 = ((lambda_min - u1) * (u1 - lambda_max) > 0) && (x.dot * d1 <= 0);
          bool ok2 = ((lambda_min - u2) * (u2 - lambda_max) > 0) && (x.dot * d2 <= 0);

          if (os && (verbosity > 2))
            (*os) << "     d1=" << d1 << " d2=" << d2
                  << " u1=" << u1 << " u2=" << u2
                  << " ok1=" << ok1 << " ok2=" << ok2 << "\n";

          Scalar olde = e;
          e = d;
          if (ok1 && ok2)
            d = (std::abs(d1) < std::abs(d2)) ? d1 : d2;
          else if (ok1 && !ok2)
            d = d1;
          else if (!ok1 && ok2)
          {
            d = d2;
            ok1 = true;
          }

          if (std::abs(d) > std::abs(0.5 * olde))
            ok1 = false;

          if (ok1)
          {
            u.size = x.size + d;
            if (((u.size - lambda_min) < tol2) || ((lambda_max - u.size) < tol2))
              d = Sign(tol1, lambda_mid - x.size);
          }

          if (os && (verbosity > 2))
            (*os) << "     olde=" << olde << " e=" << e << " d=" << d
                  << " ok1=" << ok1 << " u.size=" << u.size << "\n";
        }

        if (!ok1)
        {
          e = (x.dot > 0) ?
              (lambda_min - x.size) : (lambda_max - x.size);
          d = 0.5 * e;
        }

        if (std::abs(d) >= tol1)
          u.recalc_df(x.size + d);
        else
        {
          u.recalc_df(x.size + Sign(tol1, d));
          done = (u.f > x.f);
        }

        if (os && (verbosity > 2))
          (*os) << "     e=" << e << " d=" << d << " u=" << u
                << " done=" << done << "\n";

        if (!done)
        {
          if (u.f < x.f)
          {
            if (u.size >= x.size)
              lambda_min = x.size;
            else
              lambda_max = x.size;
            v = w;
            w = x;
            x = u;
            if (os && (verbosity > 2))
              (*os) << "     Case1  u=" << u << "  v=" << v
                    << "  w=" << w << "  x=" << x << "\n";
          }
          else
          {
            if (u.size < x.size)
              lambda_min = u.size;
            else
              lambda_max = u.size;

            if ((u.f <= w.f) || (v.size == x.size))
            {
              v = w;
              w = u;
            }
            else if ((u.f < v.f) ||
                (v.size == x.size) ||
                (v.size == w.size))
              v = u;
            if (os && (verbosity > 2))
              (*os) << "     Case2  u=" << u << "  v=" << v
                    << "  w=" << w << "  x=" << x << "\n";
          }
        }
      }

      if (done)
        break;
    }

    if (!done && os && (verbosity > 0))
      (*os) <<  "     Brent maximum number of iterations reached\n";

    // \todo how to indicate failure?
//    if (!done)
//      x.recalc_df(0.);

    return x;
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

    bracket(step_min, step_init, step_max);

    if (os && (verbosity > 0))
      (*os) << "   bracketed\n"
            << "    min =" << step_min << "\n"
            << "    init=" << step_init << "\n"
            << "    max =" << step_max << "\n";

    StepEval step = brent_search(step_init, step_min.size, step_max.size);

    if (os && (verbosity > 0))
      (*os) << "   result min =" << step << "\n";

    return step.size;
  }

};

}
