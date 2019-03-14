#include <core/fitting/optimizers/BFGS.h>
#include <core/util/more_math.h>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/tools/minima.hpp>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

double BudapestOptimizer::Sign(double a, double b)
{
  if (b >= 0)
    return std::abs(a);
  else
    return -std::abs(a);
}

double BudapestOptimizer::BrentDeriv(FittableFunction* fittable,
                        double a,
                        double b,
                        double c,
                        double tol,
                        double& xmin,
                        const Eigen::VectorXd& variables,
                        const Eigen::VectorXd& search_direction)
{
  double d{0.0}; // \todo is this really the default value?
  double d1, d2, du, dv, dw, dx;
  double fu, fv, fw, fx;
  size_t iteration;
  bool ok1, ok2, done;
  double olde, tol1, tol2, u, u1, u2, v, w, x, xm;

  double sa = (a < c) ? a : c;
  double sb = (a > c) ? a : c;

  w = x = v = b;
  fw = fv = fx = fgv(fittable, x, variables, search_direction);
  dw = dv = dx = dfgv(fittable, x, variables, search_direction);

  double e{0.0};
  // \todo check for cancel
  for (iteration = 0; iteration < brent_maximum_iterations; ++iteration)
  {
    xm = 0.5 * (sa + sb);
    tol1 = tol * std::abs(x) + brent_zeps;
    tol2 = 2 * tol1;
    done = (std::abs(x - xm) <= (tol2 - 0.5 * (sb - sa)));
    if (!done)
    {
      ok1 = false;
      if (std::abs(e) > tol1)
      {
        d2 = d1 = 2 * (sb - sa);
        if (dw != dx)
          d1 = (w - x) * dx / (dx - dw);
        if (dv != dx)
          d2 = (v - x) * dx / (dx - dv);
        u1 = x + d1;
        u2 = x + d2;
        ok1 = ((sa - u1) * (u1 - sb) > 0) && (dx * d1 <= 0);
        ok2 = ((sa - u2) * (u2 - sb) > 0) && (dx * d2 <= 0);
        olde = e;
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
          u = x + d;
          if (((u - sa) < tol2) || ((sb - u) < tol2))
            d = Sign(tol1, xm - x);
        }
      }

      if (!ok1)
      {
        e = (dx > 0) ? (sa - x) : (sb - x);
        d = 0.5 * e;
      }

      if (std::abs(d) >= tol1)
      {
        u = x + d;
        fu = fgv(fittable, u, variables, search_direction);
      }
      else
      {
        u = x + Sign(tol1, d);
        fu = fgv(fittable, u, variables, search_direction);
        done = (fu > fx);
      }

      if (!done)
      {
        du = dfgv(fittable, u, variables, search_direction);
        if (fu < fx)
        {
          if (u >= x)
            sa = x;
          else
            sb = x;
          v = w;
          fv = fw;
          dv = dw;
          w = x;
          fw = fx;
          dw = dx;
          x = u;
          fx = fu;
          dx = du;
        }
        else
        {
          if (u < x)
            sa = u;
          else
            sb = u;

          if ((fu <= fw) || (v == x))
          {
            v = w;
            fv = fw;
            dv = dw;
            w = u;
            fw = fu;
            dw = du;
          }
          else if ((fu < fv) || (v == x) || (v == w))
          {
            v = u;
            fv = fu;
            dv = du;
          }
        }
      }
    }

    if (done)
      break;
  }
  if (!done && verbosity)
    WARN("Warning: The maximum number of iterations ({}) reached in Brent line search", iteration);

  xmin = x;
  return fx;
}

void BudapestOptimizer::Bracket(FittableFunction* fittable,
                   double& a, double& b, double& c,
                   double& fa, double& fb, double& fc,
                   const Eigen::VectorXd& variables, const Eigen::VectorXd& search_direction)
{
  double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;

  double ulim, u, r, q, fu, dum;
  bool n;

  fa = fgv(fittable, a, variables, search_direction);
  fb = fgv(fittable, b, variables, search_direction);

  if (fb > fa)
  {
    dum = a;
    a = b;
    b = dum;
    dum = fb;
    fb = fa;
    fa = dum;
  }

  c = b + golden_ratio * (b - a);
  fc = fgv(fittable, c, variables, search_direction);

  while (fb > fc)
  {
    r = (b - a) * (fb - fc);
    q = (b - c) * (fb - fa);
    n = true;
    u = std::abs(q - r);
    if (bracket_tiny > u)
      u = bracket_tiny;
    if (r > q)
      u = -1 * u;
    u = b - ((b - c) * q - (b - a) * r) / (2 * u);
    ulim = b + bracket_glimit * (c - b);
    if ((b - u) * (u - c) > 0)
    {
      fu = fgv(fittable, u, variables, search_direction);
      if (fu < fc)
      {
        a = b;
        fa = fb;
        b = u;
        fb = fu;
        n = false;
      }
      else if (fu > fb)
      {
        c = u;
        fc = fu;
        n = false;
      }
      else
      {
        u = c + golden_ratio * (c - b);
        fu = fgv(fittable, u, variables, search_direction);
      }
    }
    else if ((c - u) * (u - ulim) > 0)
    {
      fu = fgv(fittable, u, variables, search_direction);
      if (fu < fc)
      {
        b = c;
        c = u;
        u = c + golden_ratio * (c - b);
        fb = fc;
        fc = fu;
        fu = fgv(fittable, u, variables, search_direction);
      }
    }
    else if ((u - ulim) * (ulim - c) >= 0)
    {
      u = ulim;
      fu = fgv(fittable, u, variables, search_direction);
    }
    else
    {
      u = c + golden_ratio * (c - b);
      fu = fgv(fittable, u, variables, search_direction);
    }

    if (n)
    {
      a = b;
      b = c;
      c = u;
      fa = fb;
      fb = fc;
      fc = fu;
    }
  }

}

double BudapestOptimizer::fgv(FittableFunction* fittable,
                 double lambda,
                 Eigen::VectorXd variables,
                 Eigen::VectorXd search_direction)
{
  return fittable->chi_sq(variables + lambda * search_direction);
}

double BudapestOptimizer::dfgv(FittableFunction* fittable,
                  double lambda,
                  Eigen::VectorXd variables,
                  Eigen::VectorXd search_direction)
{
  Eigen::VectorXd dflocal(variables.size());
  fittable->chi_sq_gradient(variables + lambda * search_direction, dflocal);
  return dflocal.dot(search_direction);
}

struct BoostFuncWrapper
{
  FittableFunction* func;
};

double BudapestOptimizer::LinMin(FittableFunction* fittable, Eigen::VectorXd& variables,
    Eigen::VectorXd search_direction)
{
  double lambdak, xk, fxk, fa, fb, a, b;
  a = 0.0;
  xk = 1.0;
  b = 2.0;
  Bracket(fittable, a, xk, b, fa, fxk, fb, variables, search_direction);
  double fmin = BrentDeriv(fittable, a, xk, b, linmin_tol, lambdak, variables, search_direction);

//  boost::math::tools::brent_find_minima<>()

  if (verbosity)
    DBG("lambda={}", lambdak);
  variables += lambdak * search_direction;
  return fmin;
}

FitResult BudapestOptimizer::minimize(FittableFunction* fittable)
{
  FitResult ret;

  ret.variables = fittable->variables();
  auto var_count = static_cast<size_t>(ret.variables.size());
  Eigen::VectorXd
      search_direction(var_count),
      gradients(var_count),
      prev_val(var_count),
      Adg(var_count);

  double f = fittable->chi_sq_gradient(ret.variables, gradients);

  ret.inv_hessian.resize(var_count, var_count);
  ret.inv_hessian.setIdentity();

  for (; ret.iterations < maximum_iterations; ++ret.iterations)
  {
    double fmin = LinMin(fittable, ret.variables, search_direction);
    //if (std::abs(f - fmin) < 0.000001) { done = true; }
    ret.converged = (2 * std::abs(fmin - f)) <= (tolerance * (std::abs(fmin) + std::abs(f) + eps));
    f = fmin;

    for (size_t i = 0; i < var_count; ++i)
      prev_val[i] = gradients[i];

    fmin = fittable->chi_sq_gradient(ret.variables, gradients);
    if (verbosity)
      INFO("Fitting... Iteration = {}, Chisq = {}", ret.iterations, fmin);

    for (size_t i = 0; i < var_count; ++i)
      prev_val[i] = gradients[i] - prev_val[i];

    for (size_t i = 0; i < var_count; ++i)
    {
      auto& Adgi = Adg[i];
      Adgi = 0.0;
      for (size_t j = 0; j < var_count; ++j)
        Adgi += ret.inv_hessian.coeff(i, j) * prev_val[j];
    }

    double s{0.0};
    double s1{0.0};
    for (size_t i = 0; i < var_count; ++i)
    {
      s1 += prev_val[i] * search_direction[i];
      s += prev_val[i] * Adg[i];
    }

    if (s1 != 0.0)
      s1 = 1.0 / s1;

    double s2{0.0};
    if (s != 0)
      s2 = 1.0 / s;

    for (size_t i = 0; i < var_count; ++i)
      prev_val[i] = s1 * search_direction[i] - s2 * Adg[i];

    for (size_t i = 0; i < var_count; ++i)
      for (size_t j = 0; j < var_count; ++j)
        ret.inv_hessian.coeffRef(i, j) +=
            s1 * square(search_direction[i])
                - s2 * square(Adg[i])
                + s * square(prev_val[i]);

    for (size_t i = 0; i < var_count; ++i)
    {
      search_direction[i] = 0;
      for (size_t j = 0; j < var_count; ++j)
        search_direction[i] -= ret.inv_hessian.coeff(i, j) * gradients[j];
    }

    if (ret.converged || cancel.load())
      break;
  }

  return ret;
}

}
