#include <core/fitting/roi.h>
#include <core/util/custom_logger.h>
#include <core/util/timer.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri {

Fit::Fit(const SUM4Edge &lb, const SUM4Edge &rb,
         PolyBackground bkg,
         const std::map<double, Peak> &peaks,
         std::string descr)
         : LB_(lb)
         , RB_(rb)
         , background(bkg)
{
  for (const auto& p : peaks)
    peaks_.push_back(p.second);

  description.description = descr;
  description.peaknum = peaks_.size();
  if (!peaks_.empty())
  {
    //description.chi_sq_norm = chi_sq_normalized();
    UncertainDouble tot_gross {0.0, 0.0};
    UncertainDouble tot_back {0.0, 0.0};
    for (auto &p : peaks_)
    {
      tot_gross += p.sum4().gross_area();
      tot_back  += p.sum4().background_area();
    }
    UncertainDouble tot_net = tot_gross - tot_back;
    description.sum4aggregate = tot_net.error();
  }
}

ROI::ROI(const FitSettings& fs,
    const Finder &parentfinder, double min, double max)
    : settings_(fs)
{
  settings_ = fs;
  set_data(parentfinder, min, max);
}

double ROI::ID() const
{
  return left_bin();
}

double ROI::left_bin() const
{
  if (finder_.x_.empty())
    return -1;
  else
    return finder_.x_.front();
}

double ROI::right_bin() const
{
  if (finder_.x_.empty())
    return -1;
  else
    return finder_.x_.back();
}

double ROI::width() const
{
  if (finder_.x_.empty())
    return 0;
  else
    return right_bin() - left_bin() + 1;
}

void ROI::set_data(const Finder &parentfinder, double l, double r)
{
  if (!finder_.cloneRange(parentfinder, l, r))
  {
    finder_.clear();
    return;
  }

  init_edges();
  init_background();
  render();
}

bool ROI::refit(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  if (peaks_.empty())
    return auto_fit(optimizer, interruptor);

  if (!rebuild(optimizer, interruptor))
    return false;

  save_current_fit("Refit");
  return true;
}


bool ROI::auto_fit(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  peaks_.clear();
  finder_.y_resid_ = finder_.y_;
  finder_.find_peaks();  //assumes default params!!!

  if (finder_.filtered.empty())
    return false;

  if ((LB_.width() == 0) || (RB_.width() == 0))
  {
    init_edges();
    init_background();
  }

  if (!settings_.sum4_only)
  {
    std::vector<double> y_nobkg = remove_background();

    for (size_t i=0; i < finder_.filtered.size(); ++i)
    {
      // \todo finder.subset
      std::vector<double> x_pk = std::vector<double>(
          finder_.x_.begin() + finder_.lefts[i],
          finder_.x_.begin() + finder_.rights[i] + 1);
      std::vector<double> y_pk = std::vector<double>(
          y_nobkg.begin() + finder_.lefts[i],
          y_nobkg.begin() + finder_.rights[i] + 1);

      auto gaussian = Hypermet().gaussian_only();
      //optimizer->fit(gaussian, x_pk, y_pk);

      if (gaussian.sanity_check(finder_.x_[finder().lefts[i]], finder_.x_[finder_.rights[i]]))
      {
        gaussian.force_defaults(default_peak_);
        Peak fitted(gaussian, {}, settings_.calib);
        peaks_[fitted.center().value()] = fitted;
      }
    }
    if (peaks_.empty())
      settings_.sum4_only = true;
  }

  if (!rebuild(optimizer, interruptor))
    return false;

  save_current_fit("Autofit");

  if (settings_.resid_auto)
    iterative_fit(optimizer, interruptor);

  return true;
}

void ROI::iterative_fit(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  if (!settings_.calib.cali_fwhm_.valid() || peaks_.empty())
    return;

  //double prev_chi_sq = peaks_.begin()->second.hypermet().chi2();

  for (int i=0; i < settings_.resid_max_iterations; ++i)
  {
    ROI new_fit = *this;

    DBG("Attempting add from resid with {} peaks", peaks_.size());

    if (!new_fit.add_from_resid(optimizer, interruptor, -1)) {
      //      DBG << "    failed add from resid";
      break;
    }
//    double new_rsq = new_fit.peaks_.begin()->second.hypermet().chi2();
//    double improvement = (prev_chi_sq - new_rsq) / prev_chi_sq * 100;
//    DBG("X2 new={} previous={} improved by {}", new_rsq, prev_chi_sq, improvement);
//
//    if ((new_rsq > prev_chi_sq) || std::isnan(new_rsq)) {
//      DBG(" not improved. reject new refit");
//      break;
//    }

    new_fit.save_current_fit("Iterative " + std::to_string(new_fit.peaks().size()));
//    prev_chi_sq = new_rsq;
    *this = new_fit;

    if (interruptor.load())
      break;
  }
}

bool ROI::add_from_resid(BFGS& optimizer, std::atomic<bool>& interruptor, int32_t centroid_hint)
{
  if (finder_.filtered.empty())
    return false;

  int64_t target_peak = -1;
  if (centroid_hint == -1) {
    double biggest = 0;
    for (size_t j=0; j < finder_.filtered.size(); ++j)
    {
      std::vector<double> x_pk = std::vector<double>(finder_.x_.begin() + finder_.lefts[j],
                                                     finder_.x_.begin() + finder_.rights[j] + 1);
      std::vector<double> y_pk = std::vector<double>(finder_.y_resid_.begin() + finder_.lefts[j],
                                                     finder_.y_resid_.begin() + finder_.rights[j] + 1);
      auto gaussian = Hypermet().gaussian_only();
      //optimizer->fit(gaussian, x_pk, y_pk);

      bool too_close = false;

      double lateral_slack = settings_.resid_too_close * gaussian.width_.val() * 2;
      for (auto &p : peaks_)
      {
        if ((p.second.center().value() > (gaussian.position.val() - lateral_slack))
            && (p.second.center().value() < (gaussian.position.val() + lateral_slack)))
          too_close = true;
      }

      //      if (too_close)
      //        DBG << "Too close at " << settings_.cali_nrg_.transform(gaussian.center_, settings_.bits_);

      if (!too_close &&
      gaussian.sanity_check(finder_.x_[finder_.lefts[j]], finder_.x_[finder_.rights[j]])
      && (gaussian.amplitude.val() > settings_.resid_min_amplitude) &&
           (gaussian.area().value() > biggest)
           )
      {
        target_peak = j;
        biggest = gaussian.area().value();
      }
    }
    //    DBG << "    biggest potential add at " << finder_.x_[finder_.filtered[target_peak]] << " with area=" << biggest;
  } else {

    //THIS NEVER HAPPENS
    double diff = std::abs(finder_.x_[finder_.filtered[target_peak]] - centroid_hint);
    for (size_t j=0; j < finder_.filtered.size(); ++j)
      if (std::abs(finder_.x_[finder_.filtered[j]] - centroid_hint) < diff) {
        target_peak = j;
        diff = std::abs(finder_.x_[finder_.filtered[j]] - centroid_hint);
      }
  }

  if (target_peak == -1)
    return false;

  std::vector<double> x_pk = std::vector<double>(finder_.x_.begin() + finder_.lefts[target_peak],
                                                 finder_.x_.begin() + finder_.rights[target_peak] + 1);
  std::vector<double> y_pk = std::vector<double>(finder_.y_resid_.begin() + finder_.lefts[target_peak],
                                                 finder_.y_resid_.begin() + finder_.rights[target_peak] + 1);
  auto gaussian = Hypermet().gaussian_only();
  //optimizer->fit(gaussian, x_pk, y_pk);

  if (gaussian.sanity_check(finder_.x_[finder().lefts[target_peak]],
      finder_.x_[finder_.rights[target_peak]]))
  {
    gaussian.apply_defaults(default_peak_);
    Peak fitted(gaussian, {}, settings_.calib);
    peaks_[fitted.center().value()] = fitted;
    rebuild(optimizer, interruptor);
    return true;
  }
  else
    return false;
}


//Peak ROI::peak(double peakID) const
//{
//  if (contains(peakID))
//    return peaks_.at(peakID);
//  else
//    return Peak();
//}


bool ROI::overlaps(double bin) const
{
  if (!width())
    return false;
  return ((bin >= left_bin()) && (bin <= right_bin()));
}

bool ROI::overlaps(double Lbin, double Rbin) const
{
  if (finder_.x_.empty())
    return false;
  if (overlaps(Lbin) || overlaps(Rbin))
    return true;
  if ((Lbin <= left_bin()) && (Rbin >= right_bin()))
    return true;
  return false;
}

bool ROI::overlaps(const ROI& other) const
{
  if (!other.width())
    return false;
  return overlaps(other.left_bin(), other.right_bin());
}

const std::map<double, Peak> &ROI::peaks() const
{
  return peaks_;
}


bool ROI::adjust_sum4(double &peakID, double left, double right)
{
  if (!peaks_.count(peakID))
    return false;

  Peak pk = peaks_.at(peakID);
  SUM4 new_sum4(finder_.weighted_data.subset(left, right), LB_, RB_);
  pk = Peak(pk.hypermet(), new_sum4, settings_.calib);
  remove_peak(peakID);
  peakID = pk.center().value();
  peaks_[peakID] = pk;
  render();
  save_current_fit("SUM4 adjusted on " + std::to_string(pk.energy().value()));
  return true;
}


bool ROI::replace_hypermet(double &peakID, Hypermet hyp)
{
  if (!peaks_.count(peakID))
    return false;

  Peak pk = peaks_.at(peakID);
  pk = Peak(hyp, pk.sum4(), settings_.calib);
  remove_peak(peakID);
  peakID = pk.center().value();
  peaks_[peakID] = pk;
  //set chi_sq_norm to 0 for all peaks?

  render();
  save_current_fit("Hypermet adjusted on " + std::to_string(pk.energy().value()));
  return true;
}

bool ROI::override_energy(double peakID, double energy)
{
  if (!peaks_.count(peakID))
    return false;

   peaks_[peakID].override_energy(energy);

   render();
   save_current_fit("Peak energy override " + std::to_string(peaks_.at(peakID).center().value())
                    + "->" + std::to_string(peaks_.at(peakID).energy().value()));
   return true;
}


bool ROI::add_peak(const Finder &parentfinder,
                   double left, double right,
                   BFGS& optimizer,
                   std::atomic<bool>& interruptor)
{
  uint16_t center_prelim = (left+right) * 0.5; //assume down the middle

  if (overlaps(left) && overlaps(right)) {
    ROI new_fit = *this;

    if (!settings_.sum4_only && new_fit.add_from_resid(optimizer, interruptor, center_prelim))
    {
      *this = new_fit;
      save_current_fit("Added from residuals");
      return true;
    }
    else
    {
      Peak fitted(Hypermet(), SUM4(finder_.weighted_data.subset(left, right), LB_, RB_), settings_.calib);
      peaks_[fitted.center().value()] = fitted;
      render();
      save_current_fit("Manually added " + std::to_string(fitted.energy().value()));
      return true;
    }
  }
  else if (width()) //changing region bounds
  {
    if (!finder_.x_.empty())
    {
      left  = std::min(left, left_bin());
      right = std::max(right, right_bin());
    }
    if (!finder_.cloneRange(parentfinder, left, right))
      return false;

    init_edges();
    init_background();
    finder_.y_resid_ = remove_background();
    render();
    finder_.find_peaks();  //assumes default params!!!

    ROI new_fit = *this;
    if (settings_.sum4_only)
    {
      Peak fitted(Hypermet(), SUM4(finder_.weighted_data.subset(left, right), LB_, RB_), settings_.calib);
      peaks_[fitted.center().value()] = fitted;
      render();
      save_current_fit("Manually added " + std::to_string(fitted.energy().value()));
      return true;
    }
    else if (new_fit.add_from_resid(optimizer, interruptor, finder_.find_index(center_prelim)))
    {
      *this = new_fit;
      save_current_fit("Added from residuals");
      return true;
    }
    else
      return auto_fit(optimizer, interruptor);
  }

  DBG("<ROI> cannot add to empty ROI");
  return false;
}

bool ROI::remove_peaks(const std::set<double> &pks, BFGS& optimizer, std::atomic<bool>& interruptor)
{
  bool found = false;
  for (auto &q : pks)
    if (remove_peak(q))
      found = true;

  if (!found)
    return false;

  if (peaks_.size() && !rebuild(optimizer, interruptor))
    return false;

  render();
  save_current_fit("Peaks removed");
  return true;
}

bool ROI::remove_peak(double bin)
{
  if (peaks_.count(bin)) {
    peaks_.erase(bin);
    return true;
  }
  return false;
}

bool ROI::override_settings(const FitSettings &fs, std::atomic<bool>& interruptor)
{
  settings_ = fs;
  settings_.overriden = true; //do this in fitter if different?
  save_current_fit("Fit settings overriden");

  //propagate to peaks

  //render if calibs changed?
  return true;
}

void ROI::save_current_fit(std::string description)
{
  Fit thisfit(LB_, RB_, background_, peaks_, description);
  fits_.push_back(thisfit);
  current_fit_ = fits_.size() - 1;
}

bool ROI::rebuild(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  rendering_.clear();

  bool hypermet_fit = false;
  for (auto &q : peaks_)
    if (!q.second.hypermet().is_gaussian_only())
    {
      hypermet_fit = true;
      break;
    }

  bool success = hypermet_fit ?
                 rebuild_as_hypermet(optimizer, interruptor)
                              :
                 rebuild_as_gaussian(optimizer, interruptor);

  if (!success)
    return false;

  render();
  return true;
}

bool ROI::rebuild_as_hypermet(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  Timer timer(true);

  std::map<double, Peak> new_peaks;

  std::vector<Hypermet> old_hype;
  for (auto &q : peaks_) {
    if (q.second.hypermet().amplitude.val())
      old_hype.push_back(q.second.hypermet());
    else if (q.second.sum4().peak_width())
    {
      Peak s4only({},
                  SUM4(finder_.weighted_data.subset(q.second.sum4().left(), q.second.sum4().right()), LB_, RB_),
                  settings_.calib);
      new_peaks[s4only.center().value()] = s4only;
    }
  }

  if (old_hype.empty())
    return false;

//  std::vector<Hypermet> hype = optimizer->fit_multiplet(finder_.x_, finder_.y_,
//                                                        old_hype, background_,
//                                                        settings_);
//
//  for (size_t i=0; i < hype.size(); ++i) {
//    double edge =  hype[i].width().value() * sqrt(log(2)) * 3; //use const from settings
//    double left = hype[i].center().value() - edge;
//    double right = hype[i].center().value() + edge;
//    Peak one(hype[i], SUM4(left, right, finder_, LB_, RB_), settings_);
//    new_peaks[one.center()] = one;
//  }

  peaks_ = new_peaks;
  return true;
}

bool ROI::rebuild_as_gaussian(BFGS& optimizer, std::atomic<bool>& interruptor)
{
  Timer timer(true);

  std::map<double, Peak> new_peaks;

  std::vector<Hypermet> old_gauss;
  for (auto &q : peaks_)
  {
    if (q.second.hypermet().amplitude.val())
      old_gauss.push_back(q.second.hypermet().gaussian_only());
    else if (q.second.sum4().peak_width())
    {
      Peak s4only({},
                  SUM4(finder_.weighted_data.subset(q.second.sum4().left(), q.second.sum4().right()), LB_, RB_),
                  settings_.calib);
      //      q.second.construct(settings_);
      new_peaks[s4only.center().value()] = s4only;
    }
  }

  if (old_gauss.empty())
    return false;

//  std::vector<Gaussian> gauss = optimizer->fit_multiplet(finder_.x_, finder_.y_,
//                                                         old_gauss, background_,
//                                                         settings_);
//
//  for (size_t i=0; i < gauss.size(); ++i)
//  {
//    double edge =  gauss[i].hwhm().value() * 3; //use const from settings
//    double left = gauss[i].center().value() - edge;
//    double right = gauss[i].center().value() + edge;
//    Peak one(Hypermet(gauss[i], settings_),
//             SUM4(left, right, finder_, LB_, RB_),
//             settings_);
//    new_peaks[one.center()] = one;
//  }

  peaks_ = new_peaks;
  return true;
}

void PeakRendering::clear()
{
  peak.clear();
  full_fit.clear();
}

void PeakRendering::render(const Hypermet& h)
{

}

void RegionRendering::reserve(size_t count)
{
  channel.assign(count, 0.0);
  energy.assign(count, 0.0);
  background.assign(count, 0.0);
  back_steps.assign(count, 0.0);
  full_fit.assign(count, 0.0);
  sum4_background.assign(count, 0.0);
  for (auto& p : peaks)
  {
    p.second.peak.assign(count, 0.0);
    p.second.full_fit.assign(count, 0.0);
  }
}

void RegionRendering::clear()
{
  channel.clear();
  energy.clear();
  background.clear();
  back_steps.clear();
  full_fit.clear();
  sum4_background.clear();
  peaks.clear();
}

void RegionRendering::sum4only(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const DAQuiri::Calibration& energy_calib,
    const DAQuiri::Polynomial& sum4back,
    const std::map<double, DAQuiri::Peak>& pks)
{
  clear();
  channel = x;
  full_fit = y;
  energy = energy_calib.transform(channel);

  for (auto& p : pks)
    peaks[p.first].full_fit = peaks[p.first].peak = full_fit;

  sum4_background = sum4back.eval(channel);
}

void RegionRendering::with_hypermet(
    double start, double end,
    const Calibration& energy_calib,
    const Polynomial& sum4back,
    const PolyBackground& hyp_back,
    const std::map<double, Peak>& pks)
{
  clear();

  auto count = static_cast<size_t>(std::ceil((end - start) * subdivisions));
  double step = 1.0 / static_cast<double>(subdivisions);

  reserve(count);

  for (size_t i = 0; i < count; ++i)
  {
    double x = start + static_cast<double>(i) * step;
    channel[i] = x;
    energy[i] = energy_calib.transform(x);
    full_fit[i] = back_steps[i] = background[i] = hyp_back.eval(x);
    for (auto& p : pks)
    {
      auto vals = p.second.hypermet().eval(x);
      back_steps[i] += vals.step_tail();
      full_fit[i] += vals.all();
    }
    sum4_background[i] = sum4back(x);
  }

  for (auto& p : pks)
  {
    auto& peak = peaks[p.first];
    const auto& hyp = p.second.hypermet();
    peak.full_fit = back_steps;
    for (size_t i = 0; i < channel.size(); ++i)
    {
      auto vals = hyp.eval(channel[i]);
      peak.peak[i] += vals.peak_skews();
      peak.full_fit[i] += vals.peak_skews();
    }
  }
}

void ROI::render()
{
  Polynomial sum4back = SUM4Edge::sum4_background(LB_, RB_);
  std::vector<double> lowres_backsteps, lowres_fullfit;

  if (settings_.sum4_only)
  {
    rendering_.sum4only(finder_.x_, finder_.y_,
                        settings_.calib.cali_nrg_, sum4back, peaks_);
    lowres_backsteps = sum4back.eval(finder_.x_);
    lowres_fullfit = sum4back.eval(finder_.x_);
  }
  else
  {
    rendering_.with_hypermet(finder_.x_.front(), finder_.x_.back(),
                             settings_.calib.cali_nrg_, sum4back,
                             background_, peaks_);
    background_.eval_add(finder_.x_, lowres_backsteps);
    background_.eval_add(finder_.x_, lowres_fullfit);
    for (auto& p : peaks_)
    {
      for (size_t i = 0; i < finder_.x_.size(); ++i)
      {
        auto vals = p.second.hypermet().eval(finder_.x_[i]);
        lowres_backsteps[i] += vals.step_tail();
        lowres_fullfit[i] += vals.all();
      }
    }
  }

  finder_.reset();
  finder_.setFit(finder_.x_, lowres_fullfit, lowres_backsteps);
}

// \todo belongs in another class?
std::vector<double> ROI::remove_background()
{
  std::vector<double> y_nobkg(finder_.x_.size());
  for (size_t i = 0; i < finder_.y_.size(); ++i)
    y_nobkg[i] = finder_.y_[i] - background_.eval(finder_.x_[i]);
  return y_nobkg;
}

bool ROI::adjust_LB(const Finder &parentfinder, double left, double right,
                     BFGS& optimizer, std::atomic<bool>& interruptor)
{
  SUM4Edge edge(parentfinder.weighted_data.subset(left, right));
  if (!edge.width() || (edge.right() >= RB_.left()))
    return false;

  if ((edge.left() != left_bin()) && !finder_.cloneRange(parentfinder, left, right_bin()))
    return false;

  LB_ = edge;
  init_background();
  cull_peaks();
  render();
  rebuild(optimizer, interruptor);
  save_current_fit("Left baseline adjusted");
  return true;
}

bool ROI::adjust_RB(const Finder &parentfinder, double left, double right,
                    BFGS& optimizer, std::atomic<bool>& interruptor) {
  SUM4Edge edge(parentfinder.weighted_data.subset(left, right));
  if (!edge.width() || (edge.left() <= LB_.right()))
    return false;

  if ((edge.right() != right_bin()) && !finder_.cloneRange(parentfinder, left_bin(), right))
    return false;

  RB_ = edge;
  init_background();
  cull_peaks();
  render();
  rebuild(optimizer, interruptor);
  save_current_fit("Right baseline adjusted");
  return true;
}

void ROI::init_edges()
{
  init_LB();
  init_RB();
}

void ROI::init_LB()
{
  LB_ = SUM4Edge(finder_.weighted_data.left(settings_.background_edge_samples));
}

void ROI::init_RB()
{
  LB_ = SUM4Edge(finder_.weighted_data.right(settings_.background_edge_samples));
}

void ROI::init_background()
{
  if (finder_.x_.empty())
    return;

  background_ = PolyBackground();
  background_.x_offset = finder_.x_.front();

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

  double maxcurve = (square(run) - std::min(LB_.min(), RB_.min())) / std::max(LB_.max(), RB_.max());
  double slope = (RB_.average() - LB_.average()) / run;

//  background_.base. set_coeff(0, {ymin, ymax, yav});
//  background_.slope. set_coeff(1, {0.5 * minslope, 2 * maxslope, slope});
//  background_.curve. set_coeff(2, {-maxcurve, maxcurve, 0});

  background_.base.val(yav);
  background_.slope.val(slope);
  background_.curve.val(0);
}

size_t ROI::current_fit() const
{
  return current_fit_;
}

std::vector<FitDescription> ROI::history() const
{
  std::vector<FitDescription> ret;
  for (auto &f : fits_)
    ret.push_back(f.description);
  return ret;
}


bool ROI::rollback(const Finder &parent_finder, size_t i)
{
  if (i >= fits_.size())
    return false;

  //settings_ = fits_[i].settings_;
  set_data(parent_finder, fits_[i].LB_.left(), fits_[i].RB_.right());
  background_ = fits_[i].background;
  LB_ = fits_[i].LB_;
  RB_ = fits_[i].RB_;
  peaks_.clear();
  for (const auto& p : fits_[i].peaks_)
    peaks_[p.center().value()] = p;
  render();

  current_fit_ = i;

  return true;
}

void ROI::cull_peaks()
{
  std::map<double, Peak> peaks;
  for (auto &p : peaks_) {
    if ((p.first > LB_.right()) &&
        (p.first < RB_.left()))
      peaks[p.first] = p.second;
  }
  peaks_ = peaks;
}

nlohmann::json ROI::to_json(const Finder &parent_finder) const
{
  nlohmann::json j;

  if (fits_.empty())
    return j;

  j["current_fit"] = current_fit_;

  ROI temp(*this);

  for (size_t i=0; i < temp.fits_.size(); ++i)
  {
    nlohmann::json jj;

    jj["description"] = temp.fits_[i].description.description;
    temp.rollback(parent_finder, i);

    if (settings_.overriden)
      jj["settings"] = settings_;

    jj["background_left"] = temp.LB();
    jj["background_right"] = temp.RB();
    jj["background_poly"] = temp.background_;

    for (auto &p : temp.peaks_)
      jj["peaks"].push_back(p.second);

    j["fits"].push_back(jj);
  }

  return j;
}

ROI::ROI(const nlohmann::json& j, const Finder &finder, const FitSettings& fs)
{
  settings_ = fs;
  if (finder.x_.empty() || (finder.x_.size() != finder.y_.size()))
    return;

  if (j.count("fits"))
  {
    nlohmann::json o = j["fits"];
    for (nlohmann::json::iterator it = o.begin(); it != o.end(); ++it)
    {
      SUM4Edge LB = it.value()["background_left"];
      SUM4Edge RB = it.value()["background_right"];

      LB = SUM4Edge(finder.weighted_data.subset(LB.left(), LB.right()));
      RB = SUM4Edge(finder.weighted_data.subset(RB.left(), RB.right()));

      if (!LB.width() || !RB.width())
        return;

      if (it.value().count("settings"))
        settings_ = it.value()["settings"];
//      else
//        settings_ = finder.settings_;

      //validate background and edges?
      set_data(finder, LB.left(), RB.left());

      LB_ = LB;
      RB_ = RB;

      if (it.value().count("peaks"))
      {
        nlohmann::json p = it.value()["peaks"];
        for (nlohmann::json::iterator it2 = p.begin(); it2 != p.end(); ++it2)
        {
          Peak newpeak(it2.value(), settings_.calib, finder_, LB_, RB_);
          peaks_[newpeak.center().value()] = newpeak;
        }
      }

      background_ = it.value()["background_poly"];
      render();
      save_current_fit(it.value()["description"]);
      peaks_.clear();
    }
  }

  rollback(finder, j["current_fit"]);
}

}
