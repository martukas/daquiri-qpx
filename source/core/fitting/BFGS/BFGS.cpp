#include <core/fitting/BFGS/BFGS.h>

#include <core/util/more_math.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

#include <core/util/custom_logger.h>

namespace DAQuiri
{

double BFGS::Sign(double a, double b)
{
  if (b >= 0)
    return std::abs(a);
  else
    return -std::abs(a);
}

double BFGS::BrentDeriv(Fittable* fittable,
                        double a,
                        double b,
                        double c,
                        double tol,
                        double& xmin,
                        const Eigen::VectorXd& variables,
                        const Eigen::VectorXd& hessian)
{
  static constexpr int32_t brent_iter{500};
  static constexpr double zeps{0.0000000001};

  double d{0.0}; // \todo is this really the default value?
  double d1, d2, du, dv, dw, dx;
  double fu, fv, fw, fx;
  int32_t iteration;
  bool ok1, ok2, done;
  double olde, tol1, tol2, u, u1, u2, v, w, x, xm;

  double sa = (a < c) ? a : c;
  double sb = (a > c) ? a : c;

  w = x = v = b;
  fw = fv = fx = fgv(fittable, x, variables, hessian);
  dw = dv = dx = dfgv(fittable, x, variables, hessian);

  double e{0.0};
  // \todo check for cancel
  for (iteration = 0; iteration < brent_iter; ++iteration)
  {
    xm = 0.5 * (sa + sb);
    tol1 = tol * std::abs(x) + zeps;
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
        fu = fgv(fittable, u, variables, hessian);
      }
      else
      {
        u = x + Sign(tol1, d);
        fu = fgv(fittable, u, variables, hessian);
        done = (fu > fx);
      }

      if (!done)
      {
        du = dfgv(fittable, u, variables, hessian);
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
  if (!done)
    WARN("Warning: The maximum number of iterations ({}) reached in Brent line search", iteration);

  xmin = x;
  return fx;
}

void BFGS::Bracket(Fittable* fittable,
                   double& a, double& b, double& c, double& fa, double& fb, double& fc,
                   const Eigen::VectorXd& variables, const Eigen::VectorXd& hessian)
{
  static constexpr double glimit{100.0};
  static constexpr double tiny{1.0E-20};

  double gold = (1.0 + std::sqrt(5.0)) / 2.0;

  double ulim, u, r, q, fu, dum;
  bool n;

  fa = fgv(fittable, a, variables, hessian);
  fb = fgv(fittable, b, variables, hessian);

  if (fb > fa)
  {
    dum = a;
    a = b;
    b = dum;
    dum = fb;
    fb = fa;
    fa = dum;
  }

  c = b + gold * (b - a);
  fc = fgv(fittable, c, variables, hessian);

  while (fb > fc)
  {
    r = (b - a) * (fb - fc);
    q = (b - c) * (fb - fa);
    n = true;
    u = std::abs(q - r);
    if (tiny > u)
      u = tiny;
    if (r > q)
      u = -1 * u;
    u = b - ((b - c) * q - (b - a) * r) / (2 * u);
    ulim = b + glimit * (c - b);
    if ((b - u) * (u - c) > 0)
    {
      fu = fgv(fittable, u, variables, hessian);
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
        u = c + gold * (c - b);
        fu = fgv(fittable, u, variables, hessian);
      }
    }
    else if ((c - u) * (u - ulim) > 0)
    {
      fu = fgv(fittable, u, variables, hessian);
      if (fu < fc)
      {
        b = c;
        c = u;
        u = c + gold * (c - b);
        fb = fc;
        fc = fu;
        fu = fgv(fittable, u, variables, hessian);
      }
    }
    else if ((u - ulim) * (ulim - c) >= 0)
    {
      u = ulim;
      fu = fgv(fittable, u, variables, hessian);
    }
    else
    {
      u = c + gold * (c - b);
      fu = fgv(fittable, u, variables, hessian);
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

double BFGS::fgv(Fittable* fittable,
                 double lambda,
                 Eigen::VectorXd variables,
                 Eigen::VectorXd hessian)
{
  auto n = variables.size();
  Eigen::VectorXd xlocal(n);
  for (size_t i = 0; i < n; ++i)
    xlocal[i] = variables[i] + lambda * hessian[i];
  return fittable->chi_sq(xlocal);
}

double BFGS::dfgv(Fittable* fittable,
                  double lambda,
                  Eigen::VectorXd variables,
                  Eigen::VectorXd hessian)
{
  auto n = variables.size();
  Eigen::VectorXd xlocal(n);
  Eigen::VectorXd dflocal(n);
  for (size_t i = 0; i < n; ++i)
    xlocal[i] = variables[i] + lambda * hessian[i];
  (*fittable)(xlocal, dflocal);
  double s = 0;
  for (size_t i = 0; i < n; ++i)
    s += dflocal[i] * hessian[i];
  return s;
}

double BFGS::LinMin(Fittable* fittable, Eigen::VectorXd& variables,
    Eigen::VectorXd hessian)
{
  static constexpr float linmin_tol{0.0001};
  double lambdak, xk, fxk, fa, fb, a, b;
  auto n = variables.size();
  a = 0.0;
  xk = 1.0;
  b = 2.0;
  Bracket(fittable, a, xk, b, fa, fxk, fb, variables, hessian);
  double fmin = BrentDeriv(fittable, a, xk, b, linmin_tol, lambdak, variables, hessian);
  DBG("lambda={}", lambdak);
  for (size_t j = 0; j < n; ++j)
  {

    hessian[j] *= lambdak;
    variables[j] += hessian[j];
  }
  return fmin;
}

FitResult BFGS::BFGSMin(Fittable* fittable, double tolf)
{
  static constexpr double eps{0.0000000001};
  static constexpr size_t maxit{500};

  FitResult ret;

  ret.variables = fittable->variables();
  auto var_count = ret.variables.size();
  double free_vars = fittable->degrees_of_freedom();
  Eigen::VectorXd
      hessian(var_count),
      gradients(var_count),
      prev_val(var_count),
      Adg(var_count);

  double f = (*fittable)(ret.variables, gradients);

  ret.inv_hessian.resize(var_count, var_count);
  ret.inv_hessian.setIdentity();

  for (; ret.iterations < maxit; ++ret.iterations)
  {
    double fmin = LinMin(fittable, ret.variables, hessian);
    //if (std::abs(f - fmin) < 0.000001) { done = true; }
    ret.converged = (2 * std::abs(fmin - f)) <= (tolf * (std::abs(fmin) + std::abs(f) + eps));
    f = fmin;

    for (size_t i = 0; i < var_count; ++i)
      prev_val[i] = gradients[i];

    fmin = (*fittable)(ret.variables, gradients);
    INFO("Fitting... Iteration = {}, Chisq = {}", ret.iterations, fmin / free_vars);

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
      s1 += prev_val[i] * hessian[i];
      s += prev_val[i] * Adg[i];
    }

    if (s1 != 0.0)
      s1 = 1.0 / s1;

    double s2{0.0};
    if (s != 0)
      s2 = 1.0 / s;

    for (size_t i = 0; i < var_count; ++i)
      prev_val[i] = s1 * hessian[i] - s2 * Adg[i];

    for (size_t i = 0; i < var_count; ++i)
      for (size_t j = 0; j < var_count; ++j)
        ret.inv_hessian.coeffRef(i, j) +=
            s1 * square(hessian[i])
                - s2 * square(Adg[i])
                + s * square(prev_val[i]);

    for (size_t i = 0; i < var_count; ++i)
    {
      hessian[i] = 0;
      for (size_t j = 0; j < var_count; ++j)
        hessian[i] -= ret.inv_hessian.coeff(i, j) * gradients[j];
    }

    if (ret.converged || cancel.load())
      break;
  }

  return ret;
}

}
