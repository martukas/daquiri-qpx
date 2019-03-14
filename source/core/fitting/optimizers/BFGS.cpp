#include <core/fitting/optimizers/BFGS.h>
#include <core/util/more_math.h>

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

BudapestOptimizer::StepEval::StepEval(FittableFunction* fittable,
    const Eigen::VectorXd& variables,
    const Eigen::VectorXd& search_direction)
    : fittable_(fittable)
    , variables_(&variables)
    , search_direction_(&search_direction)
{}

BudapestOptimizer::StepEval& BudapestOptimizer::StepEval::operator=(const StepEval& other)
{
  fittable_ = other.fittable_;
  variables_ = other.variables_;
  search_direction_ = other.search_direction_;
  step = other.step;
  f = other.f;
  dot = other.dot;
  return *this;
}

void BudapestOptimizer::StepEval::recalc_df(double lambda)
{
  step = lambda;
  Eigen::VectorXd dflocal(variables_->size());
  f = fittable_->chi_sq_gradient((*variables_) + lambda * (*search_direction_), dflocal);
  dot = dflocal.dot(*search_direction_);
}

void BudapestOptimizer::StepEval::recalc_f(double lambda)
{
  step = lambda;
  f = fittable_->chi_sq((*variables_) + lambda * (*search_direction_));
}

BudapestOptimizer::StepEval BudapestOptimizer::BrentDeriv(FittableFunction *fittable,
                                                          double lambda1,
                                                          double initial_lambda,
                                                          double lambda2,
                                                          const Eigen::VectorXd &variables,
                                                          const Eigen::VectorXd &search_direction)
{
  StepEval step_x(fittable, variables, search_direction),
           step_v(fittable, variables, search_direction),
           step_w(fittable, variables, search_direction),
           step_u(fittable, variables, search_direction);
  step_x.recalc_df(initial_lambda);
  step_w = step_v = step_x;

  double lambda_min = std::min(lambda1, lambda2);
  double lambda_max = std::max(lambda1, lambda2);

  bool done {false};
  double e{0.0};
  double d{0.0}; // \todo is this really the default value?

  // \todo check for cancel
  for (size_t iteration = 0; iteration < brent_maximum_iterations; ++iteration)
  {
    double lambda_mid = 0.5 * (lambda_min + lambda_max);
    double tol1 = linmin_tolerance * std::abs(step_x.step) + brent_zeps;
    double tol2 = 2. * tol1;

    done = (std::abs(step_x.step - lambda_mid) <= (tol2 - 0.5 * (lambda_max - lambda_min)));

    if (!done)
    {
      bool ok1 = false;
      if (std::abs(e) > tol1)
      {
        double d1 = 2 * (lambda_max - lambda_min);
        double d2 = d1;
        if (step_w.dot != step_x.dot)
          d1 = (step_w.step - step_x.step) * step_x.dot / (step_x.dot - step_w.dot);
        if (step_v.dot != step_x.dot)
          d2 = (step_v.step - step_x.step) * step_x.dot / (step_x.dot - step_v.dot);
        double u1 = step_x.step + d1;
        double u2 = step_x.step + d2;
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
          step_u.step = step_x.step + d;
          if (((step_u.step - lambda_min) < tol2) || ((lambda_max - step_u.step) < tol2))
            d = Sign(tol1, lambda_mid - step_x.step);
        }
      }

      if (!ok1)
      {
        e = (step_x.dot > 0) ? (lambda_min - step_x.step) : (lambda_max - step_x.step);
        d = 0.5 * e;
      }

      if (std::abs(d) >= tol1)
        step_u.recalc_df(step_x.step + d);
      else
      {
        step_u.recalc_df(step_x.step + Sign(tol1, d));
        done = (step_u.f > step_x.f);
      }

      if (!done)
      {
        if (step_u.f < step_x.f)
        {
          if (step_u.step  >= step_x.step)
            lambda_min = step_x.step;
          else
            lambda_max = step_x.step;
          step_v = step_w;
          step_w = step_x;
          step_x = step_u;
        }
        else
        {
          if (step_u.step < step_x.step)
            lambda_min = step_u.step;
          else
            lambda_max = step_u.step;

          if ((step_u.f <= step_w.f) || (step_v.step == step_x.step))
          {
            step_v = step_w;
            step_w = step_u;
          }
          else if ((step_u.f < step_v.f) || (step_v.step == step_x.step) || (step_v.step == step_w.step))
          {
            step_v = step_u;
          }
        }
      }
    }

    if (done)
      break;
  }

  if (!done && verbosity)
    WARN("Warning: The maximum number of iterations reached in Brent line search");

  // \todo how to indicate failure?

  return step_x;
}

void BudapestOptimizer::Bracket(FittableFunction* fittable,
                                StepEval& a_step, StepEval& b_step, StepEval& c_step,
                                const Eigen::VectorXd& variables, const Eigen::VectorXd& search_direction)
{
  double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;

  StepEval u_step(fittable, variables, search_direction);

  double ulim, r, q;
  bool n;

  a_step.recalc_f(a_step.step);
  b_step.recalc_f(b_step.step);

  if (b_step.f > a_step.f)
    std::swap(a_step, b_step);

  c_step.recalc_f(b_step.step + golden_ratio * (b_step.step - a_step.step));

  while (b_step.f > c_step.f)
  {
    r = (b_step.step - a_step.step) * (b_step.f - c_step.f);
    q = (b_step.step - c_step.step) * (b_step.f - a_step.f);
    n = true;
    u_step.step = std::abs(q - r);
    if (bracket_tiny > u_step.step)
      u_step.step = bracket_tiny;
    if (r > q)
      u_step.step = -u_step.step;
    u_step.step = b_step.step - ((b_step.step - c_step.step) * q - (b_step.step - a_step.step) * r) / (2 * u_step.step);
    ulim = b_step.step + bracket_glimit * (c_step.step - b_step.step);
    if ((b_step.step - u_step.step) * (u_step.step - c_step.step) > 0)
    {
      u_step.recalc_f(u_step.step);

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
      {
        u_step.recalc_f(c_step.step + golden_ratio * (c_step.step - b_step.step));

      }
    }
    else if ((c_step.step - u_step.step) * (u_step.step - ulim) > 0)
    {
      u_step.recalc_f(u_step.step);
      if (u_step.f < c_step.f)
      {
        b_step = c_step;
        c_step = u_step;
        u_step.recalc_f(c_step.step + golden_ratio * (c_step.step - b_step.step));
      }
    }
    else if ((u_step.step - ulim) * (ulim - c_step.step) >= 0)
    {
      u_step.recalc_f(ulim);
    }
    else
    {
      u_step.recalc_f(c_step.step + golden_ratio * (c_step.step - b_step.step));
    }

    if (n)
    {
      a_step = b_step;
      b_step = c_step;
      c_step = u_step;
    }
  }

}

double BudapestOptimizer::LinMin(FittableFunction* fittable, Eigen::VectorXd& variables,
                                 const Eigen::VectorXd& search_direction)
{
  StepEval lambda_min(fittable, variables, search_direction),
           lambda_max(fittable, variables, search_direction),
           lambda_init(fittable, variables, search_direction);
  lambda_min.step = 0.0;
  lambda_init.step = 1.0;
  lambda_max.step = 2.0;
  Bracket(fittable, lambda_min, lambda_init, lambda_max, variables, search_direction);
  StepEval step = BrentDeriv(fittable, lambda_min.step, lambda_init.step, lambda_max.step, variables, search_direction);

  if (verbosity)
    DBG("lambda={}", step.step);
  variables += step.step * search_direction;
  return step.f;
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
