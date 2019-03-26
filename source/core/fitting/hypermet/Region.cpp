#include <core/fitting/hypermet/Region.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

Region::Region()
{
//  default_peak_ = default_peak_.gaussian_only();
}

Region::Region(const WeightedData& new_data, uint16_t background_samples)
    : Region()
{
  if (new_data.data.empty())
    throw std::runtime_error("Attempting to construct Region from empty sample");

  data = new_data;
  LB_ = SUM4Edge(data.left(background_samples));
  RB_ = SUM4Edge(data.right(background_samples));
  init_background();
}

void Region::replace_data(const WeightedData& new_data)
{
  auto lb = SUM4Edge(new_data.subset(LB_.left(), LB_.right()));
  auto rb = SUM4Edge(new_data.subset(RB_.left(), RB_.right()));
  replace_data(new_data, lb, rb);
}

void Region::replace_data(const WeightedData& new_data, uint16_t left_samples, uint16_t right_samples)
{
  replace_data(new_data, SUM4Edge(new_data.left(left_samples)), SUM4Edge(new_data.right(right_samples)));
}

void Region::replace_data(const WeightedData& new_data, const SUM4Edge& lb, const SUM4Edge& rb)
{
  if (new_data.empty())
    throw std::runtime_error("Attempting to construct Region from empty sample");

  if ((lb.left() < new_data.data.front().channel) || (new_data.data.back().channel < rb.right()))
    throw std::runtime_error("Sum4 edges outside of region");

  // \todo only if edge values changed
  LB_ = lb;
  RB_ = rb;
  data = new_data;

  init_background();
  cull_peaks();
}

void Region::adjust_LB(const SUM4Edge& lb)
{
  if ((lb.left() < left()) || (lb.right() > right()))
    return;

  // \todo only if edge values changed
  LB_ = lb;
  init_background();
}

void Region::adjust_RB(const SUM4Edge& rb)
{
  if ((rb.left() < left()) || (rb.right() > right()))
    return;

  // \todo only if edge values changed
  RB_ = rb;
  init_background();
}

double Region::left() const
{
  if (!data.empty())
    return data.data.front().channel;
  return std::numeric_limits<double>::quiet_NaN();
}

double Region::right() const
{
  if (!data.empty())
    return data.data.back().channel;
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
  double center = 0.5 * (l + r);
  double radius = (r - l) * 0.5;
  p.position.bound(center - 0.5 * radius, center + 0.5 * radius);
  p.position.val(center);

  data.subset(l, r);
  double max_val{0.0};
  double min_bkg{std::numeric_limits<double>::max()};
  for (const auto& v : data.data)
  {
    max_val = std::max(max_val, v.count);
    min_bkg = std::min(min_bkg, background.eval(v.channel));
  }
  double amp_max = max_val - min_bkg;

  //p.amplitude.bound(0, 1.1 * amp_max);
  if (amp_hint > 0.0)
    p.amplitude.val(amp_hint);
  else
    p.amplitude.val(amp_max);
  peaks_[p.id()] = p;
  dirty_ = true;
  return true;
}

bool Region::adjust_sum4(double peakID, double left, double right)
{
  if (!peaks_.count(peakID))
    return false;
  auto subregion = data.subset(left, right);
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
  {
    if (p.second.sanity_check(left(), right()))
      peaks[p.first] = p.second;
  }
  if (peaks.size() != peaks_.size())
    dirty_ = true;
  peaks_ = peaks;
}

void Region::init_background()
{
  background = PolyBackground(data, LB_, RB_);
  // \todo invalidate uncertanties

  dirty_ = true;
}

void Region::reindex_peaks()
{
  std::map<double, Peak> new_peaks;
  for (auto& p : peaks_)
    new_peaks[p.second.id()] = p.second;
  peaks_ = new_peaks;
}

void Region::update_indices()
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

  variable_count = 0;
  background.update_indices(variable_count);

  default_peak_.reset_indices();
  if (!peaks_.empty())
  {
    if (unique_widths < peaks_.size())
      default_peak_.width.update_index(variable_count);

    if (default_peak_.short_tail.enabled &&
        (unique_short_tails < peaks_.size()))
      default_peak_.short_tail.update_indices(variable_count);

    if (default_peak_.right_tail.enabled &&
        (unique_right_tails < peaks_.size()))
      default_peak_.right_tail.update_indices(variable_count);

    if (default_peak_.long_tail.enabled &&
        (unique_long_tails < peaks_.size()))
      default_peak_.long_tail.update_indices(variable_count);

    if (default_peak_.step.enabled &&
        (unique_steps < peaks_.size()))
      default_peak_.step.update_indices(variable_count);
  }

  for (auto& p : peaks_)
  {
    p.second.apply_defaults(default_peak_);
    p.second.update_indices(variable_count);
  }
}

Eigen::VectorXd Region::variables() const
{
  Eigen::VectorXd ret;
  ret.setConstant(variable_count, 0.0);

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
}

void Region::save_fit(const FitResult& result)
{
  save_fit(result.variables);
  //auto chi_sq_norm = chi_sq(result.variables);

  if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
    return;

  Eigen::VectorXd diagonals;
  diagonals.resize(result.variables.size());

  double df = degrees_of_freedom();
  for (size_t i = 0; i < static_cast<size_t>(result.variables.size()); ++i)
    diagonals[i] = result.inv_hessian.coeff(i, i) * df;

  // \todo should this be done here?
  // \todo check if optimizer does this 0.5 already?
  double chisq_norm = std::max(this->chi_sq(result.variables) / df, 1.0) * 0.5;

  background.get_uncerts(diagonals, chisq_norm);
  if (!peaks_.empty())
    default_peak_.get_uncerts(diagonals, chisq_norm);

  for (auto& p : peaks_)
    p.second.get_uncerts(diagonals, chisq_norm);

  dirty_ = !result.converged;
}

double Region::eval(double chan) const
{
  double ret = background.eval(chan);
  for (const auto& p : peaks_)
    ret += p.second.eval(chan).all();
  return ret;
}

double Region::eval_at(double chan, const Eigen::VectorXd& fit) const
{
  double ret = background.eval_at(chan, fit);
  for (const auto& p : peaks_)
    ret += p.second.eval_at(chan, fit).all();
  return ret;
}

double Region::eval_grad_at(double chan, const Eigen::VectorXd& fit, Eigen::VectorXd& grads) const
{
  double ret = background.eval_grad_at(chan, fit, grads);
  for (const auto& p : peaks_)
    ret += p.second.eval_grad_at(chan, fit, grads).all();
  return ret;
}

bool Region::perturb(std::mt19937& rng)
{
//  if (background.base.valid_index())
//    background.base.x(background.base.x() + x_dist(rng));
//  if (background.slope.valid_index())
//    background.slope.x(background.slope.x() + x_dist(rng));
//  if (background.curve.valid_index())
//    background.curve.x(background.curve.x() + x_dist(rng));

  for (auto& p : peaks_)
  {
    if (p.second.width.valid_index())
      p.second.width.x(x_dist(rng));
    if (p.second.position.valid_index())
      p.second.position.x(x_dist(rng));
    if (p.second.amplitude.valid_index())
      p.second.amplitude.x(p.second.amplitude.x() + x_dist(rng));

    if (p.second.short_tail.amplitude.valid_index())
      p.second.short_tail.amplitude.x(x_dist(rng));
    if (p.second.short_tail.slope.valid_index())
      p.second.short_tail.slope.x(x_dist(rng));

    if (p.second.right_tail.amplitude.valid_index())
      p.second.right_tail.amplitude.x(x_dist(rng));
    if (p.second.right_tail.slope.valid_index())
      p.second.right_tail.slope.x(x_dist(rng));

    if (p.second.long_tail.amplitude.valid_index())
      p.second.long_tail.amplitude.x(x_dist(rng));
    if (p.second.long_tail.slope.valid_index())
      p.second.long_tail.slope.x(x_dist(rng));

//    if (p.second.step.amplitude.valid_index())
//      p.second.step.amplitude.x(x_dist(rng));
  }

  return true;
}

bool Region::sane() const
{
  for (const auto& p : peaks_)
    if (!p.second.sane())
      return false;
  return true;
}

std::string Region::to_string(std::string prepend) const
{
  std::stringstream ss;
  ss << prepend
     << fmt::format("data_points={} on channels [{},{}] vars={} {}\n",
                               data.data.size(), left(), right(),
                               variable_count, (dirty_ ? " DIRTY" : ""));
  ss << prepend << "chi2: " << chi_sq() << "\n";
  ss << prepend << "chi2/dof: " << chi_sq() / degrees_of_freedom() << "\n";
  ss << prepend << "SUM4.LB: " << LB_.to_string() << "\n";
  ss << prepend << "SUM4.RB: " << RB_.to_string() << "\n";
  ss << prepend << "Background:\n";
  ss << background.to_string(prepend + " ");
  ss << prepend << "Common:\n";
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
