#include <core/fitting/hypermet/Region.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

Region::Region(const SpectrumData& data)
  : data_(data)
{
  if (data_.data.empty())
    throw std::runtime_error("Attempting to construct Region from empty sample");

  background.x_offset = data_.data.front().x;

  default_peak_.width_.bound(0.8, 4.0);

  default_peak_.short_tail.amplitude.bound(0.02, 1.5);
  default_peak_.short_tail.amplitude.to_fit = true;
  default_peak_.short_tail.slope.bound(0.2, 0.5);
  default_peak_.short_tail.slope.to_fit = true;

  default_peak_.right_tail.amplitude.bound(0.01, 0.9);
  default_peak_.right_tail.amplitude.to_fit = true;
  default_peak_.right_tail.slope.bound(0.3, 1.5);
  default_peak_.right_tail.slope.to_fit = true;

  default_peak_.long_tail.amplitude.bound(0.0001, 0.15);
  default_peak_.long_tail.amplitude.to_fit = true;
  default_peak_.long_tail.slope.bound(2.5, 50);
  default_peak_.long_tail.slope.to_fit = true;

  default_peak_.step.amplitude.bound(0.000001, 0.05);
  default_peak_.step.amplitude.to_fit = true;
}

void Region::map_fit()
{
  size_t unique_widths{0};
  size_t unique_short_tails{0};
  size_t unique_right_tails{0};
  size_t unique_long_tails{0};
  size_t unique_steps{0};
  for (auto& p : peaks_)
  {
    if (p.width_override)
      unique_widths++;
    if (p.short_tail.override)
      unique_short_tails++;
    if (p.right_tail.override)
      unique_right_tails++;
    if (p.long_tail.override)
      unique_long_tails++;
    if (p.step.override)
      unique_steps++;
  }

  var_count_ = 0;
  background.update_indices(var_count_);

  if (unique_widths < peaks_.size())
    default_peak_.width_.x_index = var_count_++;
  else
    default_peak_.width_.x_index = -1;

  if (default_peak_.short_tail.enabled &&
      (unique_short_tails < peaks_.size()))
    default_peak_.short_tail.update_indices(var_count_);

  if (default_peak_.right_tail.enabled &&
      (unique_right_tails < peaks_.size()))
    default_peak_.right_tail.update_indices(var_count_);

  if (default_peak_.long_tail.enabled &&
      (unique_long_tails < peaks_.size()))
    default_peak_.long_tail.update_indices(var_count_);

  if (default_peak_.step.enabled &&
      (unique_steps < peaks_.size()))
    default_peak_.step.update_indices(var_count_);

  for (auto& p : peaks_)
  {
    p.apply_defaults(default_peak_);
    p.update_indices(var_count_);
  }
}

size_t Region::fit_var_count() const
{
  return static_cast<size_t>(var_count_);
}

std::vector<double> Region::variables() const
{
  std::vector<double> ret;
  ret.resize(fit_var_count(), 0.0);

  background.put(ret);
  default_peak_.put(ret);

  for (auto& p : peaks_)
    p.put(ret);

  return ret;
}

void Region::save_fit(const std::vector<double>& variables)
{
  background.get(variables);
  default_peak_.get(variables);

  for (auto& p : peaks_)
    p.get(variables);
}

void Region::save_fit_uncerts(const FitResult& result)
{
  save_fit(result.variables);

  std::vector<double> diagonals;
  diagonals.reserve(result.variables.size());

  double df = degrees_of_freedom();
  for (size_t i = 0; i < result.variables.size(); ++i)
    diagonals.push_back(result.inv_hessian.coeff(i, i) * df);

  double chisq_norm = std::max(chi_sq_normalized(), 1.0) * 0.5;

  background.get_uncerts(diagonals, chisq_norm);
  default_peak_.get_uncerts(diagonals, chisq_norm);

  for (auto& p : peaks_)
    p.get_uncerts(diagonals, chisq_norm);
}

double Region::chi_sq_normalized() const
{
  return chi_sq() / degrees_of_freedom();
}

double Region::degrees_of_freedom() const
{
  // \todo what if channel range is < fit_var_count?
  return (data_.data.size() - fit_var_count());
}

//Calculates the Chi-square over a region
double Region::chi_sq() const
{
  double ChiSq = 0;
  for (const auto& data : data_.data)
  {
    double FTotal = background.eval(data.x);
    for (auto& p : peaks_)
      FTotal += p.eval(data.x).all();
    ChiSq += square((data.y - FTotal) / data.weight_true);
  }
  return ChiSq;
}


//Calculates the Chi-square and its gradient
double Region::grad_chi_sq(std::vector<double>& gradients) const
{
  gradients.assign(gradients.size(), 0.0);
  auto chan_gradients = gradients;

  double Chisq = 0;

  for (const auto& data : data_.data)
  {
    chan_gradients.assign(chan_gradients.size(), 0.0);

    double FTotal = background.eval_grad(data.x, chan_gradients);
    for (auto& p : peaks_)
      FTotal += p.eval_grad(data.x, chan_gradients).all();

    double t3 = -2.0 * (data.y - FTotal) / square(data.weight_true);
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((data.y - FTotal) / data.weight_true);
  }
  //Chisq /= df

  return Chisq;
}

//Calculates the Chi-square over a region
double Region::chi_sq(const std::vector<double>& fit) const
{
  double ChiSq = 0;
  for (const auto& data : data_.data)
  {
    double FTotal = background.eval_at(data.x, fit);
    for (auto& p : peaks_)
      FTotal += p.eval_at(data.x, fit).all();
    ChiSq += square((data.y - FTotal) / data.weight_phillips_marlow);
  }
  return ChiSq;
}

//Calculates the Chi-square and its gradient
double Region::grad_chi_sq(const std::vector<double>& fit,
                           std::vector<double>& gradients) const
{
  gradients.assign(gradients.size(), 0.0);
  auto chan_gradients = gradients;

  double Chisq = 0;

  for (const auto& data : data_.data)
  {
    chan_gradients.assign(chan_gradients.size(), 0.0);

    double FTotal = background.eval_grad_at(data.x, fit, chan_gradients);
    for (auto& p : peaks_)
      FTotal += p.eval_grad_at(data.x, fit, chan_gradients).all();

    double t3 = -2.0 * (data.y - FTotal) / square(data.weight_phillips_marlow);
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((data.y - FTotal) / data.weight_phillips_marlow);
  }
  //Chisq /= df

  return Chisq;
}

}
