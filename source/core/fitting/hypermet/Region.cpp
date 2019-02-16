#include <core/fitting/hypermet/Region.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

#include <dlib/matrix/matrix_mat.h>

namespace DAQuiri
{

Region::Region()
{
  default_peak_.width.bound(0.8, 4.0);

  default_peak_.short_tail.amplitude.bound(0.02, 1.5);
  default_peak_.short_tail.slope.bound(0.2, 0.5);

  default_peak_.right_tail.amplitude.bound(0.01, 0.9);
  default_peak_.right_tail.slope.bound(0.3, 1.5);

  default_peak_.long_tail.amplitude.bound(0.0001, 0.15);
  default_peak_.long_tail.slope.bound(2.5, 50);

  default_peak_.step.amplitude.bound(0.000001, 0.05);
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
  auto lb = SUM4Edge(data.subset(LB_.left(), LB_.right()));
  auto rb = SUM4Edge(data.subset(RB_.left(), RB_.right()));
  replace_data(data, lb, rb);
}

void Region::replace_data(const WeightedData& data, uint16_t left_samples, uint16_t right_samples)
{
  replace_data(data, SUM4Edge(data.left(left_samples)), SUM4Edge(data.right(right_samples)));
}

void Region::replace_data(const WeightedData& data, const SUM4Edge& lb, const SUM4Edge& rb)
{
  if (data.data.empty())
    throw std::runtime_error("Attempting to construct Region from empty sample");

  if ((lb.left() < data.data.front().channel) || (data.data.back().channel < rb.right()))
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
    return data_.data.front().channel;
  return std::numeric_limits<double>::quiet_NaN();
}

double Region::right() const
{
  if (!data_.empty())
    return data_.data.back().channel;
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
    max_val = std::max(max_val, v.count);
    min_bkg = std::min(min_bkg, background.eval(v.channel));
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
  double global_min = data_.data.front().count;
  for (const auto& d : data_.data)
    global_min = std::min(global_min, d.count);
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
  background.x_offset = data_.data.front().channel;

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

  reindex_peaks();
}

void Region::save_fit(const FitResult& result)
{
  save_fit(result.variables);
  //auto chi_sq_norm = chi_sq(result.variables);

  Eigen::VectorXd diagonals;
  diagonals.resize(result.variables.size());

  double df = degrees_of_freedom();
  for (size_t i = 0; i < static_cast<size_t>(result.variables.size()); ++i)
    diagonals[i] = result.inv_hessian.coeff(i, i) * df;

  double chisq_norm = std::max(this->chi_sq(result.variables), 1.0) * 0.5;

  background.get_uncerts(diagonals, chisq_norm);
  if (!peaks_.empty())
    default_peak_.get_uncerts(diagonals, chisq_norm);

  for (auto& p : peaks_)
    p.second.get_uncerts(diagonals, chisq_norm);

  dirty_ = false;
}

double Region::eval(double chan) const
{
  double ret = background.eval(chan);
  for (auto& p : peaks_)
    ret += p.second.eval(chan).all();
  return ret;
}

double Region::eval_at(double chan, const Eigen::VectorXd& fit) const
{
  double ret = background.eval_at(chan, fit);
  for (auto& p : peaks_)
    ret += p.second.eval_at(chan, fit).all();
  return ret;
}

double Region::eval_grad_at(double chan, const Eigen::VectorXd& fit, Eigen::VectorXd& grads) const
{
  double ret = background.eval_grad_at(chan, fit, grads);
  for (auto& p : peaks_)
    ret += p.second.eval_grad_at(chan, fit, grads).all();
  return ret;
}

std::string Region::to_string(std::string prepend) const
{
  std::stringstream ss;
  ss << prepend
     << fmt::format("data_points={} on channels [{},{}] vars={} {}\n",
                               data_.data.size(), left(), right(),
                               variable_count, (dirty_ ? " DIRTY" : ""));
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
