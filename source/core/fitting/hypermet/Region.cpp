#include <core/fitting/hypermet/Region.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

Region::Region()
{
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

Region::Region(const SpectrumData& data, uint16_t background_samples)
  : Region()
{
  if (data.data.empty())
    throw std::runtime_error("Attempting to construct Region from empty sample");

  data_ = data;
  LB_ = SUM4Edge(data_.left(background_samples));
  RB_ = SUM4Edge(data_.right(background_samples));
  init_background();
}

void Region::replace_data(const SpectrumData& data)
{
  replace_data(data, LB_.width(), RB_.width());
}

void Region::replace_data(const SpectrumData& data, uint16_t left_samples, uint16_t right_samples)
{
  replace_data(data, SUM4Edge(data.left(left_samples)), SUM4Edge(data.right(right_samples)));
}

void Region::replace_data(const SpectrumData& data, const SUM4Edge& lb, const SUM4Edge& rb)
{
  if (data.data.empty())
    throw std::runtime_error("Attempting to construct Region from empty sample");

  if ((lb.left() < data.data.front().x) || (data.data.back().x < rb.right()))
    throw std::runtime_error("Sum4 edges outside of region");

  LB_ = lb;
  RB_ = rb;
  data_ = data;

  // \todo only if edge values changed
  init_background();
  cull_peaks();
  dirty = true;
}

void Region::adjust_LB(const SUM4Edge& lb)
{
  if ((lb.left() < left()) || (lb.right() > right()))
    return;
  LB_ = lb;
  init_background();
  dirty = true;
}

void Region::adjust_RB(const SUM4Edge& rb)
{
  if ((rb.left() < left()) || (rb.right() > right()))
    return;
  RB_ = rb;
  init_background();
  dirty = true;
}

double Region::left() const
{
  return data_.data.front().x;
}

double Region::right() const
{
  return data_.data.back().x;
}

bool Region::adjust_sum4(double peakID, double left, double right)
{
  if (!peaks_.count(peakID))
    return false;
  auto subregion = data_.subset(left, right);
  if (subregion.data.empty())
    return false;

  peaks_[peakID].sum4 = SUM4(subregion, LB_, RB_);
  return true;
}

bool Region::replace_hypermet(double peakID, Peak hyp)
{
  if (!peaks_.count(peakID))
    return false;

  // \todo sanity check of new peak definition?

  // \todo should this happen?
  hyp.sum4 = peaks_[peakID].sum4;
  peaks_[hyp.id()] = hyp;
  reindex_peaks();
  dirty = true;
  return true;
}

bool Region::remove_peak(double peakID)
{
  if (!peaks_.count(peakID))
    return false;
  peaks_.erase(peakID);
  dirty = true;
  return true;
}

bool Region::remove_peaks(const std::set<double>& ids)
{
  bool found = false;
  for (auto& q : ids)
    if (remove_peak(q))
      found = true;
  return found;
}

void Region::cull_peaks()
{
  std::map<double, Peak> peaks;
  for (const auto& p : peaks_)
    if (p.second.sanity_check(left(), right()))
      peaks[p.first] = p.second;
  peaks_ = peaks;
}

void Region::init_background()
{
  background = PolyBackground();
  background.x_offset = data_.data.front().x;

  // \todo make this more rigorous

  //by default, linear
  double run = RB_.left() - LB_.right();

  double minslope = 0, maxslope = 0;
  double ymin, ymax, yav;
  if (LB_.average() < RB_.average())
  {
    run = RB_.right() - LB_.right();
    minslope = (RB_.min() - LB_.max()) / (RB_.right() - LB_.left());
    maxslope = (RB_.max() - LB_.min()) / (RB_.left() - LB_.right());
    ymin = LB_.min();
    ymax = RB_.max();
    yav = LB_.average();
  }
  else if (RB_.average() < LB_.average())
  {
    run = RB_.left() - LB_.left();
    minslope = (RB_.min() - LB_.max()) / (RB_.left() - LB_.right());
    maxslope = (RB_.max() - LB_.min()) / (RB_.right() - LB_.left());
    ymin = RB_.min();
    ymax = LB_.max();
    yav = RB_.average();
  }

  double slope = (RB_.average() - LB_.average()) / run;
  double maxcurve = (square(run) - std::min(LB_.min(), RB_.min())) / std::max(LB_.max(), RB_.max());

  // \todo bounds for polynomial
//  background.base. set_coeff(0, {ymin, ymax, yav});
//  background.slope. set_coeff(1, {0.5 * minslope, 2 * maxslope, slope});
//  background.curve. set_coeff(2, {-maxcurve, maxcurve, 0});

  background.base.val(yav);
  background.slope.val(slope);
  background.curve.val(0.0);

  // \todo invalidate uncertanties
}

void Region::reindex_peaks()
{
  std::map<double, Peak> new_peaks;
  for (const auto& p : peaks_)
    new_peaks[p.second.id()] = p.second;
  peaks_ = new_peaks;
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
    if (p.second.width_override)
      unique_widths++;
    if (p.second.short_tail.override)
      unique_short_tails++;
    if (p.second.right_tail.override)
      unique_right_tails++;
    if (p.second.long_tail.override)
      unique_long_tails++;
    if (p.second.step.override)
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
    p.second.apply_defaults(default_peak_);
    p.second.update_indices(var_count_);
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
    p.second.put(ret);

  return ret;
}

void Region::save_fit(const std::vector<double>& variables)
{
  background.get(variables);
  default_peak_.get(variables);

  for (auto& p : peaks_)
    p.second.get(variables);

  reindex_peaks();
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
    p.second.get_uncerts(diagonals, chisq_norm);
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
      FTotal += p.second.eval(data.x).all();
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
      FTotal += p.second.eval_grad(data.x, chan_gradients).all();

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
      FTotal += p.second.eval_at(data.x, fit).all();
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
      FTotal += p.second.eval_grad_at(data.x, fit, chan_gradients).all();

    double t3 = -2.0 * (data.y - FTotal) / square(data.weight_phillips_marlow);
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((data.y - FTotal) / data.weight_phillips_marlow);
  }
  //Chisq /= df

  return Chisq;
}

void to_json(nlohmann::json& j, const Region& s)
{

}

void from_json(const nlohmann::json& j, Region& s)
{

}


}
