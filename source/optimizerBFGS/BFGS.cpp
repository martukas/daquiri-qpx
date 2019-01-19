#include <optimizerBFGS/BFGS.h>

#include <optimizerBFGS/more_math.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

#include <core/util/custom_logger.h>

namespace Hypermet
{

double BFGS::Sign(double a, double b)
{
  if (b >= 0)
    return std::abs(a);
  else
    return -1 * std::abs(a);
}

static constexpr int32_t brent_iter{500};
static constexpr double zeps{0.0000000001};

double BFGS::BrentDeriv(const Region& region,
                        double a,
                        double b,
                        double c,
                        double tol,
                        double& xmin,
                        const std::vector<double>& gx,
                        const std::vector<double>& gh)
{
  double sa, sb, d, d1, d2, du, dv, dw, dx;
  double e, fu, fv, fw, fx;
  int32_t k;
  bool ok1, ok2, done;
  double olde, tol1, tol2, u, u1, u2, v, w, x, xm;
  try
  {
    if (a < c)
      sa = a;
    else
      sa = c;

    if (a > c)
      sb = a;
    else
      sb = c;

    v = b;
    w = v;
    x = v;
    e = 0;
    fx = fgv(region, x, gx, gh);
    fv = fx;
    fw = fx;
    dx = dfgv(region, x, gx, gh);
    dv = dx;
    dw = dx;

    for (k = 0; k <= brent_iter; ++k)
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
          d1 = 2 * (sb - sa);
          d2 = d1;
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
          {
            if (std::abs(d1) < std::abs(d2))
              d = d1;
            else
              d = d2;
          }
          if (ok1 && !ok2)
            d = d1;
          if (!ok1 && ok2)
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
            {
              d = Sign(tol1, xm - x);
            }
          }
        }

        if (!ok1)
        {
          if (dx > 0)
            e = sa - x;
          else
            e = sb - x;
          d = 0.5 * e;
        }

        if (std::abs(d) >= tol1)
        {
          u = x + d;
          fu = fgv(region, u, gx, gh);
        }
        else
        {
          u = x + Sign(tol1, d);
          fu = fgv(region, u, gx, gh);
          done = (fu > fx);
        }
        if (!done)
        {
          du = dfgv(region, u, gx, gh);
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

            if ((fu <= fw) || (v = x))
            {
              v = w;
              fv = fw;
              dv = dw;
              w = u;
              fw = fu;
              dw = du;
            }
            else if ((fu < fv) || (v = x) || (v = w))
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
      WARN("Warning: The maximum number of iterations reached in Brent line search");

    xmin = x;
    return fx;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error(""));
  }
}

static constexpr double gold{1.618034};
static constexpr double glimit{100.0};
static constexpr double tiny{1.0E-20};

void BFGS::Bracket(const Region& region, double& a, double& b, double& c, double& fa, double& fb, double& fc,
                   const std::vector<double>& gx, const std::vector<double>& gh)
{
  double ulim, u, r, q, fu, dum;
  bool n;
  try
  {
    fa = fgv(region, a, gx, gh);
    fb = fgv(region, b, gx, gh);

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
    fc = fgv(region, c, gx, gh);

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
        fu = fgv(region, u, gx, gh);
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
          fu = fgv(region, u, gx, gh);
        }
      }
      else if ((c - u) * (u - ulim) > 0)
      {
        fu = fgv(region, u, gx, gh);
        if (fu < fc)
        {
          b = c;
          c = u;
          u = c + gold * (c - b);
          fb = fc;
          fc = fu;
          fu = fgv(region, u, gx, gh);
        }
      }
      else if ((u - ulim) * (ulim - c) >= 0)
      {
        u = ulim;
        fu = fgv(region, u, gx, gh);
      }
      else
      {
        u = c + gold * (c - b);
        fu = fgv(region, u, gx, gh);
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
  catch (...)
  {
    std::throw_with_nested(std::runtime_error(""));
  }
}

double BFGS::fgv(const Region& region, double lambda, std::vector<double> gx, std::vector<double> gh)
{
  try
  {
    auto n = gx.size();
    std::vector<double> xlocal(n);
    for (size_t i = 0; i < n; ++i)
      xlocal[i] = gx[i] + lambda * gh[i];
    return region.calc_chi_sq(xlocal);
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error(""));
  }
}

double BFGS::dfgv(const Region& region, double lambda, std::vector<double> gx, std::vector<double> gh)
{
  auto n = gx.size();
  std::vector<double> xlocal(n);
  std::vector<double> dflocal(n);
  double s, junk;
  try
  {
    for (size_t i = 0; i < n; ++i)
      xlocal[i] = gx[i] + lambda * gh[i];
    region.grad_chi_sq(xlocal, dflocal, junk);
    s = 0;
    for (size_t i = 0; i < n; ++i)
      s += dflocal[i] * gh[i];
    return s;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error(""));
  }
}

static constexpr float linmin_tol{0.0001};

void BFGS::LinMin(const Region& region, std::vector<double>& x, std::vector<double> h, double& fmin)
{
  double lambdak, xk, fxk, fa, fb, a, b;
  auto n = x.size();
  try
  {
    a = 0.0;
    xk = 1.0;
    b = 2.0;
    Bracket(region, a, xk, b, fa, fxk, fb, x, h);
    fmin = BrentDeriv(region, a, xk, b, linmin_tol, lambdak, x, h);
    DBG("lambda={}", lambdak);
    for (size_t j = 0; j < n; ++j)
    {

      h[j] *= lambdak;
      x[j] += h[j];
    }
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error(""));
  }
}

static constexpr double eps{0.0000000001};
static constexpr size_t maxit{500};

void BFGS::BFGSMin(Region& region, double tolf, size_t& iter)
{
  try
  {
    auto x = region.fit;
    auto n = x.size();
    double f, fmin, s, s1, s2;
    double fv = region.degrees_of_freedom();
    std::vector<double> h(n), g(n), u(n), Adg(n);
    bool done{false};

    region.grad_chi_sq(x, g, f);

    Eigen::SparseMatrix<double> A(n, n);
    A.setIdentity();

    size_t k;
    for (k = 0; k <= maxit; ++k)
    {
      LinMin(region, x, h, fmin);
      //if (std::abs(f - fmin) < 0.000001) { done = true; }
      done = (2 * std::abs(fmin - f)) <= (tolf * (std::abs(fmin) + std::abs(f) + eps));
      f = fmin;

      for (size_t i = 0; i < n; ++i)
        u[i] = g[i];

      region.grad_chi_sq(x, g, fmin);
      INFO("Fitting... Iteration = {}, Chisq = {}", k, fmin / fv);

      for (size_t i = 0; i < n; ++i)
        u[i] = g[i] - u[i];

      for (size_t i = 0; i < n; ++i)
      {
        s = 0;
        for (size_t j = 0; j < n; ++j)
          s += A.coeff(i, j) * u[j];
        Adg[i] = s;
      }
      s1 = 0;
      s = 0;

      for (size_t i = 0; i < n; ++i)
      {
        s1 += u[i] * h[i];
        s += u[i] * Adg[i];
      }

      if (s1 != 0)
        s1 = 1 / s1;
      if (s != 0)
        s2 = 1 / s;
      else
        s2 = 0;

      for (size_t i = 0; i < n; ++i)
        u[i] = s1 * h[i] - s2 * Adg[i];

      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
          A.coeffRef(i, j) += s1 * h[i] * h[j] - s2 * Adg[i] * Adg[j] + s * u[i] * u[j];

      for (size_t i = 0; i < n; ++i)
      {
        h[i] = 0;
        for (size_t j = 0; j < n; ++j)
          h[i] -= A.coeff(i, j) * g[j];
      }

      if (done || Cancelfit)
        break;
    }
    if (!done && !Cancelfit)
    {
      WARN("Warning: The number of iterations reached the limit without convergence");
    }
    else
    {
      if (!done && Cancelfit)
        WARN("Warning: The fit was interrupted");
      region.fit = x;
      iter = k;
    }
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error(""));
  }
}

}
