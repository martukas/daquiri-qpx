#include <core/fitting/hypermet/Region.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

#include <dlib/matrix/matrix_mat.h>

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

Region::Region(const WeightedData& data, uint16_t background_samples)
    : Region()
{
  if (data.data.empty())
    throw std::runtime_error("Attempting to construct Region from empty sample");

  data_ = data;
  LB_ = SUM4Edge(data_.left(background_samples));
  RB_ = SUM4Edge(data_.right(background_samples));
  init_background();
}

void Region::replace_data(const WeightedData& data)
{
  // \todo preserve position, not just width
  replace_data(data, LB_.width(), RB_.width());
}

void Region::replace_data(const WeightedData& data, uint16_t left_samples, uint16_t right_samples)
{
  replace_data(data, SUM4Edge(data.left(left_samples)), SUM4Edge(data.right(right_samples)));
}

void Region::replace_data(const WeightedData& data, const SUM4Edge& lb, const SUM4Edge& rb)
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
  dirty_ = true;
}

void Region::adjust_LB(const SUM4Edge& lb)
{
  if ((lb.left() < left()) || (lb.right() > right()))
    return;
  LB_ = lb;
  init_background();
  dirty_ = true;
}

void Region::adjust_RB(const SUM4Edge& rb)
{
  if ((rb.left() < left()) || (rb.right() > right()))
    return;
  RB_ = rb;
  init_background();
  dirty_ = true;
}

double Region::left() const
{
  if (!data_.empty())
    return data_.data.front().x;
  return std::numeric_limits<double>::quiet_NaN();
}

double Region::right() const
{
  if (!data_.empty())
    return data_.data.back().x;
  return std::numeric_limits<double>::quiet_NaN();
}

bool Region::empty() const
{
  return peaks_.empty();
}

bool Region::dirty() const
{
  return dirty_;
}

bool Region::add_peak(double l, double r, double amp_hint)
{
  if ((l < left()) || (right() < r))
    throw std::runtime_error("Attempting to add peak outside of region bounds");
  Peak p = default_peak_;
  p.position.bound(l, r);
  p.position.val(0.5 * (l + r));

  data_.subset(l, r);
  double max_val{0.0};
  double min_bkg{std::numeric_limits<double>::max()};
  for (const auto& v : data_.data)
  {
    max_val = std::max(max_val, v.y);
    min_bkg = std::min(min_bkg, background.eval(v.x));
  }
  double amp_max = max_val - min_bkg;

  p.amplitude.bound(0, 1.1 * amp_max);
  if (amp_hint > 0.0)
    p.amplitude.val(amp_hint);
  else
    p.amplitude.val(0.9 * amp_max);
  peaks_[p.id()] = p;
  dirty_ = true;
  return true;
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

void Region::auto_sum4()
{
  for (auto& p : peaks_)
  {
    if (p.second.sum4.peak_width())
      continue;
    // \todo use const from settings?
    // \todo do we really need to multiply with sqrt(log(2))?
    double edge = p.second.fwhm().value() * sqrt(log(2)) * 3;
    double left = p.second.peak_position().value() - edge;
    double right = p.second.peak_position().value() + edge;
    adjust_sum4(p.second.id(), left, right);
  }
}

bool Region::replace_hypermet(double peakID, Peak hyp)
{
  if (!peaks_.count(peakID))
    return false;

  // \todo sanity check of new peak definition?

  // \todo should this happen?
  hyp.sum4 = peaks_[peakID].sum4;
  peaks_.erase(peakID);
  peaks_[hyp.id()] = hyp;
  reindex_peaks();
  dirty_ = true;
  return true;
}

bool Region::remove_peak(double peakID)
{
  if (!peaks_.count(peakID))
    return false;
  peaks_.erase(peakID);
  dirty_ = true;
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
  double global_min = data_.data.front().y;
  for (const auto& d : data_.data)
    global_min = std::min(global_min, d.y);
  global_min = std::floor(global_min);

  //by default, linear
  double run = RB_.right() - LB_.left();

  INFO("run={}", run);

//  // ascending slope
//  if (LB_.average() < RB_.average())
//    run = RB_.right() - LB_.right();
//  else
//    run = RB_.left() - LB_.left();

  double ymax = std::max(LB_.max(), RB_.max());
  double ymin = std::min(LB_.min(), RB_.min());

  INFO("ymax={} ymin={}", ymax, ymin);

  double maxcurve = std::abs((ymax - ymin) / square(run));

  INFO("maxcurve={}", maxcurve);

  double minslope{0.0}, maxslope{0.0};
  double yav;

  // ascending slope
  if (LB_.average() < RB_.average())
  {
    //minslope = (RB_.min() - LB_.max()) / (RB_.right() - LB_.left());
    maxslope = (RB_.max() - LB_.min()) / (RB_.left() - LB_.right());
    yav = LB_.average().value();
  }
  else
  {
    minslope = (RB_.min() - LB_.max()) / (RB_.left() - LB_.right());
    //maxslope = (RB_.max() - LB_.min()) / (RB_.right() - LB_.left());
    yav = RB_.average().value();
  }

  background = PolyBackground();
  background.x_offset = data_.data.front().x;

  background.base.bound(global_min, ymax);
  background.slope.bound(minslope, maxslope);
  background.curve.bound(-maxcurve, maxcurve);

  background.base.val(yav);
  background.slope.val((RB_.average().value() - LB_.average().value()) / run);
  background.curve.val(0.0);

  // \todo invalidate uncertanties

  dirty_ = true;
}

void Region::reindex_peaks()
{
  std::map<double, Peak> new_peaks;
  for (const auto& p : peaks_)
    new_peaks[p.second.id()] = p.second;
  peaks_ = new_peaks;
}

size_t Region::variable_count() const
{
  return static_cast<size_t>(variable_count_);
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

  variable_count_ = 0;
  background.update_indices(variable_count_);

  if (!peaks_.empty())
  {
    if (unique_widths < peaks_.size())
      default_peak_.width_.x_index = variable_count_++;
    else
      default_peak_.width_.x_index = -1;

    if (default_peak_.short_tail.enabled &&
        (unique_short_tails < peaks_.size()))
      default_peak_.short_tail.update_indices(variable_count_);

    if (default_peak_.right_tail.enabled &&
        (unique_right_tails < peaks_.size()))
      default_peak_.right_tail.update_indices(variable_count_);

    if (default_peak_.long_tail.enabled &&
        (unique_long_tails < peaks_.size()))
      default_peak_.long_tail.update_indices(variable_count_);

    if (default_peak_.step.enabled &&
        (unique_steps < peaks_.size()))
      default_peak_.step.update_indices(variable_count_);
  }

  for (auto& p : peaks_)
  {
    p.second.apply_defaults(default_peak_);
    p.second.update_indices(variable_count_);
  }
}

size_t Region::fit_var_count() const
{
  return static_cast<size_t>(variable_count_);
}

Eigen::VectorXd Region::variables() const
{
  Eigen::VectorXd ret;
  ret.setConstant(fit_var_count(), 0.0);

  background.put(ret);

  if (!peaks_.empty())
    default_peak_.put(ret);

  for (auto& p : peaks_)
    p.second.put(ret);

  return ret;
}

void Region::save_fit(const Eigen::VectorXd& variables)
{
  background.get(variables);

  if (!peaks_.empty())
    default_peak_.get(variables);

  for (auto& p : peaks_)
    p.second.get(variables);

  reindex_peaks();
}

void Region::save_fit_uncerts(const FitResult& result)
{
  save_fit(result.variables);

  Eigen::VectorXd diagonals;
  diagonals.resize(result.variables.size());

  double df = degrees_of_freedom();
  for (size_t i = 0; i < static_cast<size_t>(result.variables.size()); ++i)
    diagonals[i] = result.inv_hessian.coeff(i, i) * df;

  double chisq_norm = std::max(chi_sq_normalized(), 1.0) * 0.5;

  background.get_uncerts(diagonals, chisq_norm);
  if (!peaks_.empty())
    default_peak_.get_uncerts(diagonals, chisq_norm);

  for (auto& p : peaks_)
    p.second.get_uncerts(diagonals, chisq_norm);

  dirty_ = false;
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
  return ChiSq / degrees_of_freedom();
}

//Calculates the Chi-square and its gradient
double Region::grad_chi_sq(Eigen::VectorXd& gradients) const
{
  gradients.setConstant(gradients.size(), 0.0);
  Eigen::VectorXd chan_gradients;
  chan_gradients.setConstant(gradients.size(), 0.0);
  double Chisq = 0;
  for (const auto& data : data_.data)
  {
    chan_gradients.setConstant(gradients.size(), 0.0);

    double FTotal = background.eval_grad(data.x, chan_gradients);
    for (auto& p : peaks_)
      FTotal += p.second.eval_grad(data.x, chan_gradients).all();

    double t3 = -2.0 * (data.y - FTotal) / square(data.weight_true);
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((data.y - FTotal) / data.weight_true);
  }
  //Chisq /= df

  return Chisq / degrees_of_freedom();
}

//Calculates the Chi-square over a region
double Region::chi_sq(const Eigen::VectorXd& fit) const
{
  double ChiSq = 0;
  for (const auto& data : data_.data)
  {
    double FTotal = background.eval_at(data.x, fit);
    for (auto& p : peaks_)
      FTotal += p.second.eval_at(data.x, fit).all();
    ChiSq += square((data.y - FTotal) / data.weight_phillips_marlow);
  }
  return ChiSq / degrees_of_freedom();
}

//Calculates the Chi-square and its gradient
double Region::operator ()(const Eigen::VectorXd& fit,
                           Eigen::VectorXd& gradients) const
{
  gradients.setConstant(fit.size(), 0.0);
  Eigen::VectorXd chan_gradients;
  chan_gradients.setConstant(fit.size(), 0.0);
  double Chisq{0.0};
  for (const auto& data : data_.data)
  {
    chan_gradients.setConstant(fit.size(), 0.0);

    double FTotal = background.eval_grad_at(data.x, fit, chan_gradients);
    for (auto& p : peaks_)
      FTotal += p.second.eval_grad_at(data.x, fit, chan_gradients).all();

    double t3 = -2.0 * (data.y - FTotal) / square(data.weight_phillips_marlow);
    for (size_t var = 0; var < fit_var_count(); ++var)
      gradients[var] += chan_gradients[var] * t3;
    Chisq += square((data.y - FTotal) / data.weight_phillips_marlow);
  }
  //Chisq /= df
  return Chisq / degrees_of_freedom();
}

double Region::eval(const fitter_vector& m) const
{
  Eigen::VectorXd v;
  v.setConstant(m.size(), 0.0);
  for (long i = 0; i < m.size(); ++i)
    v[i] = m(i);
  return chi_sq(v);
}

fitter_vector Region::derivative(const fitter_vector& m) const
{
  Eigen::VectorXd v;
  v.setConstant(m.size(), 0.0);
  for (long i = 0; i < m.size(); ++i)
    v[i] = m(i);
  Eigen::VectorXd g;
  g.setConstant(g.size(), 0.0);
  operator()(v, g);
  return dlib::mat(g);
}

fitter_matrix Region::hessian(const fitter_vector& m) const
{
  // \todo implement this
  (void) m;
  return fitter_matrix();
}

std::string Region::to_string(std::string prepend) const
{
  std::stringstream ss;
  ss << prepend
     << fmt::format("data_points={} on channels [{},{}] vars={} {}\n",
                               data_.data.size(), left(), right(),
                               variable_count_, (dirty_ ? " DIRTY" : ""));
  ss << prepend << "SUM4/LB: " << LB_.to_string() << "\n";
  ss << prepend << "SUM4/RB: " << RB_.to_string() << "\n";
  ss << prepend << "Background:\n";
  ss << background.to_string(prepend + " ");
  ss << prepend << "Default peak:\n";
  ss << default_peak_.to_string(prepend + " ");
  if (!peaks_.empty())
  {
    ss << prepend << "Peaks:\n";
    for (const auto& p : peaks_)
    {
      ss << prepend << " Peak at " << p.first << "\n";
      ss << p.second.to_string(prepend + "  ");
    }
  }
  return ss.str();
}

void to_json(nlohmann::json& j, const Region& s)
{
  j["background"] = s.background;
  j["default_peak"] = s.default_peak_;
  j["LB"] = s.LB_;
  j["RB"] = s.RB_;
  if (!s.peaks_.empty())
  {
    auto& peaks = j["peaks"];
    for (const auto& p : s.peaks_)
      peaks.push_back(p.second);
  }

  // \todo dirty, data, var_count
}

void from_json(const nlohmann::json& j, Region& s)
{
  s.background = j["background"];
  s.default_peak_ = j["default_peak"];
  s.LB_ = j["LB"];
  s.RB_ = j["RB"];
  if (j.count("peaks"))
  {
    const auto& peaks = j["peaks"];
    for (const auto& p : peaks)
    {
      Peak peak = p;
      s.peaks_[peak.id()] = peak;
    }
  }
}

}
