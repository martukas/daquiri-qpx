#include <core/fitting/optimizers/brent.h>

namespace DAQuiri
{

inline double Sign(double a, double b)
{
  if (b >= 0)
    return std::abs(a);
  else
    return -std::abs(a);
}

Brent::StepEval::StepEval(FittableFunction* fittable,
    const Eigen::VectorXd& variables,
    const Eigen::VectorXd& search_direction, double lambda)
    : fittable_(fittable)
    , variables_(&variables)
    , search_direction_(&search_direction)
    , size (lambda)
{}

Brent::StepEval& Brent::StepEval::operator=(const StepEval& other)
{
  fittable_ = other.fittable_;
  variables_ = other.variables_;
  search_direction_ = other.search_direction_;
  size = other.size;
  f = other.f;
  dot = other.dot;
  return *this;
}

void Brent::StepEval::recalc_f(double lambda)
{
  size = lambda;
  f = fittable_->chi_sq((*variables_) + lambda * (*search_direction_));
}

void Brent::StepEval::recalc_df(double lambda)
{
  size = lambda;
  Eigen::VectorXd gradient(variables_->size());
  f = fittable_->chi_sq_gradient((*variables_) + lambda * (*search_direction_), gradient);
  dot = gradient.dot(*search_direction_);
}

void Brent::StepEval::recalc_f()
{
  recalc_f(size);
}

void Brent::StepEval::recalc_df()
{
  recalc_df(size);
}

Brent::StepEval Brent::brent_search(StepEval step_x, double lambda1, double lambda2)
{
  step_x.recalc_df();
  StepEval step_u = step_x;
  StepEval step_w = step_x;
  StepEval step_v = step_x;

  step_u.recalc_df(0);

  double lambda_min = std::min(lambda1, lambda2);
  double lambda_max = std::max(lambda1, lambda2);

  bool done {false};
  double e{0.0};
  double d{0.0}; // \todo is this really the default value?

  // \todo check for cancel
  for (size_t iteration = 0; iteration < brent_maximum_iterations; ++iteration)
  {
    double lambda_mid = 0.5 * (lambda_min + lambda_max);
    double tol1 = linmin_tolerance * std::abs(step_x.size) + brent_zeps;
    double tol2 = 2. * tol1;

    done = (std::abs(step_x.size - lambda_mid) <= (tol2 - 0.5 * (lambda_max - lambda_min)));

    if (!done)
    {
      bool ok1 = false;
      if (std::abs(e) > tol1)
      {
        double d1 = 2 * (lambda_max - lambda_min);
        double d2 = d1;
        if (step_w.dot != step_x.dot)
          d1 = (step_w.size - step_x.size) * step_x.dot / (step_x.dot - step_w.dot);
        if (step_v.dot != step_x.dot)
          d2 = (step_v.size - step_x.size) * step_x.dot / (step_x.dot - step_v.dot);
        double u1 = step_x.size + d1;
        double u2 = step_x.size + d2;
        ok1 = ((lambda_min - u1) * (u1 - lambda_max) > 0) && (step_x.dot * d1 <= 0);
        bool ok2 = ((lambda_min - u2) * (u2 - lambda_max) > 0) && (step_x.dot * d2 <= 0);

        double olde = e;
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
          step_u.size = step_x.size + d;
          if (((step_u.size - lambda_min) < tol2) || ((lambda_max - step_u.size) < tol2))
            d = Sign(tol1, lambda_mid - step_x.size);
        }
      }

      if (!ok1)
      {
        e = (step_x.dot > 0) ? (lambda_min - step_x.size) : (lambda_max - step_x.size);
        d = 0.5 * e;
      }

      if (std::abs(d) >= tol1)
        step_u.recalc_df(step_x.size + d);
      else
      {
        step_u.recalc_df(step_x.size + Sign(tol1, d));
        done = (step_u.f > step_x.f);
      }

      if (!done)
      {
        if (step_u.f < step_x.f)
        {
          if (step_u.size  >= step_x.size)
            lambda_min = step_x.size;
          else
            lambda_max = step_x.size;
          step_v = step_w;
          step_w = step_x;
          step_x = step_u;
        }
        else
        {
          if (step_u.size < step_x.size)
            lambda_min = step_u.size;
          else
            lambda_max = step_u.size;

          if ((step_u.f <= step_w.f) || (step_v.size == step_x.size))
          {
            step_v = step_w;
            step_w = step_u;
          }
          else if ((step_u.f < step_v.f) || (step_v.size == step_x.size) || (step_v.size == step_w.size))
            step_v = step_u;
        }
      }
    }

    if (done)
      break;
  }

//  if (!done && verbosity)
//    WARN("Warning: The maximum number of iterations reached in Brent line search");

  // \todo how to indicate failure?

  return step_x;
}

void Brent::bracket(StepEval& a_step, StepEval& b_step, StepEval& c_step)
{
  double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;

  StepEval u_step = b_step;

  double ulim, r, q;
  bool n;

  a_step.recalc_f();
  b_step.recalc_f();

  if (b_step.f > a_step.f)
    std::swap(a_step, b_step);

  c_step.recalc_f(b_step.size + golden_ratio * (b_step.size - a_step.size));

  while (b_step.f > c_step.f)
  {
    r = (b_step.size - a_step.size) * (b_step.f - c_step.f);
    q = (b_step.size - c_step.size) * (b_step.f - a_step.f);
    n = true;
    u_step.size = std::abs(q - r);
    if (bracket_tiny > u_step.size)
      u_step.size = bracket_tiny;
    if (r > q)
      u_step.size = -u_step.size;
    u_step.size = b_step.size - ((b_step.size - c_step.size) * q - (b_step.size - a_step.size) * r) / (2 * u_step.size);
    ulim = b_step.size + bracket_glimit * (c_step.size - b_step.size);
    if ((b_step.size - u_step.size) * (u_step.size - c_step.size) > 0)
    {
      u_step.recalc_f();

      if (u_step.f < c_step.f)
      {
        a_step = b_step;
        b_step = u_step;
        n = false;
      }
      else if (u_step.f > b_step.f)
      {
        c_step = u_step;
        n = false;
      }
      else
        u_step.recalc_f(c_step.size + golden_ratio * (c_step.size - b_step.size));
    }
    else if ((c_step.size - u_step.size) * (u_step.size - ulim) > 0)
    {
      u_step.recalc_f();
      if (u_step.f < c_step.f)
      {
        b_step = c_step;
        c_step = u_step;
        u_step.recalc_f(c_step.size + golden_ratio * (c_step.size - b_step.size));
      }
    }
    else if ((u_step.size - ulim) * (ulim - c_step.size) >= 0)
      u_step.recalc_f(ulim);
    else
      u_step.recalc_f(c_step.size + golden_ratio * (c_step.size - b_step.size));

    if (n)
    {
      a_step = b_step;
      b_step = c_step;
      c_step = u_step;
    }
  }

}

double Brent::line_search(FittableFunction* fittable,
    const Eigen::VectorXd& variables,
    const Eigen::VectorXd& search_direction)
{
  StepEval step_min(fittable, variables, search_direction, 0.0),
           step_init(fittable, variables, search_direction, 1.0),
           step_max(fittable, variables, search_direction, 2.0);

  bracket(step_min, step_init, step_max);

  StepEval step = brent_search(step_init, step_min.size, step_max.size);

  return step.size;
}

}
