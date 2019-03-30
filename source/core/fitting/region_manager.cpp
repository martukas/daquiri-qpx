#include <core/fitting/region_manager.h>
#include <core/util/more_math.h>
#include <core/fitting/finders/finder_kon_naive.h>

#include <core/util/custom_logger.h>

namespace DAQuiri {

Fit::Fit(const Region& r, std::string descr)
    : region(r)
{
  description.description = descr;
  description.peaknum = region.peaks_.size();
  if (!region.peaks_.empty())
  {
    description.chi_sq_norm = region.chi_sq() / region.degrees_of_freedom();
    UncertainDouble tot_gross {0.0, 0.0};
    UncertainDouble tot_back {0.0, 0.0};
    for (const auto &p : region.peaks_)
    {
      tot_gross += p.second.sum4.gross_area();
      tot_back  += p.second.sum4.background_area();
    }
    UncertainDouble tot_net = tot_gross - tot_back;
    description.sum4aggregate = tot_net.error();
  }
}


void PeakRendering::clear()
{
  peak.clear();
  full_fit.clear();
}

void PeakRendering::render(const Peak& /*h*/)
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

void RegionRendering::render(const Region& r,
                             const Calibration& energy_calib)
{
  auto sum4back = SUM4Edge::sum4_background(r.LB_, r.RB_);

  clear();

  double start = r.left();
  double end = r.right();

  auto count = static_cast<size_t>(std::ceil((end - start) * subdivisions));
  double step = 1.0 / static_cast<double>(subdivisions);

  reserve(count);

  for (size_t i = 0; i < count; ++i)
  {
    double x = start + static_cast<double>(i) * step;
    channel[i] = x;
    energy[i] = energy_calib.transform(x);
    full_fit[i] = back_steps[i] = background[i] = r.background.eval(x);
    for (auto& p : r.peaks_)
    {
      auto vals = p.second.eval(x);
      back_steps[i] += vals.step_tail();
      full_fit[i] += vals.all();
    }
    sum4_background[i] = sum4back(x);
  }

  for (auto& p : r.peaks_)
  {
    auto& peak = peaks[p.first];
    const auto& hyp = p.second;
    peak.full_fit = back_steps;
    peak.peak.assign(count, 0.0);
    for (size_t i = 0; i < channel.size(); ++i)
    {
      auto vals = hyp.eval(channel[i]);
      peak.peak[i] += vals.peak_skews();
      peak.full_fit[i] += vals.peak_skews();
    }
  }
}


RegionManager::RegionManager(const FitSettings& fs,
    const FitEvaluation &parentfinder, double min, double max)
    : settings_(fs)
{
  settings_ = fs;
  region_.default_peak_ = settings_.default_peak;
  set_data(parentfinder, min, max);
  save_current_fit("Region created");
}

double RegionManager::id() const
{
  return left_bin();
}

bool RegionManager::dirty() const
{
  return region_.dirty();
}

double RegionManager::left_bin() const
{
  if (fit_eval_.x_.empty())
    return -1;
  else
    return fit_eval_.x_.front();
}

double RegionManager::right_bin() const
{
  if (fit_eval_.x_.empty())
    return -1;
  else
    return fit_eval_.x_.back();
}

double RegionManager::width() const
{
  if (fit_eval_.x_.empty())
    return 0;
  else
    return right_bin() - left_bin() + 1;
}

void RegionManager::set_data(const FitEvaluation &parentfinder, double l, double r)
{
  fit_eval_.cloneRange(parentfinder, l, r);
  region_ = Region(fit_eval_.weighted_data, settings_.background_edge_samples);
  render();
}

bool RegionManager::refit(AbstractOptimizer* optimizer)
{
  if (region_.peaks_.empty())
    return find_and_fit(optimizer);
  else
    return rebuild(optimizer);
}

bool RegionManager::find_and_fit(AbstractOptimizer* optimizer)
{
  fit_eval_.reset();
  NaiveKON kon(fit_eval_.x_, fit_eval_.y_,
      settings_.kon_settings.width,
      settings_.kon_settings.sigma_spectrum);
  region_ = Region(fit_eval_.weighted_data, settings_.background_edge_samples);

  if (kon.filtered.empty())
  {
    render();
    return false;
  }

  //std::vector<double> y_nobkg = remove_background();

  for (const auto& p : kon.filtered)
  {
    double height = kon.highest_residual(p.left, p.right);
    height -= region_.background.eval(0.5 * (p.left + p.right));
    region_.add_peak(p.left, p.right, height);
  }
  save_current_fit("Autofind");

  if (!rebuild(optimizer))
  {
    render();
    return false;
  }

  if (settings_.resid_auto)
    iterative_fit(optimizer);

  return true;
}

void RegionManager::iterative_fit(AbstractOptimizer* optimizer)
{
  if (!settings_.calib.cali_fwhm_.valid() || region_.peaks_.empty())
    return;

  //double prev_chi_sq = peaks_.begin()->second.hypermet().chi2();

  for (int i=0; i < settings_.resid_max_iterations; ++i)
  {
    DBG("Attempting add from resid with {} peaks", region_.peaks_.size());

    if (!add_from_resid(optimizer)) {
      //      DBG << "    failed add from resid";
      break;
    }

//    if (optimizer.cancel.load())
//      break;
  }
}

bool RegionManager::add_from_resid(AbstractOptimizer* optimizer)
{
  if (fit_eval_.empty())
    return false;

  NaiveKON kon(fit_eval_.x_, fit_eval_.y_resid_,
               settings_.kon_settings.width,
               settings_.kon_settings.sigma_resid);

  DetectedPeak target_peak = kon.tallest_detected();

  if (target_peak.highest_y == 0.0)
    return false;

  if (!region_.add_peak(target_peak.left, target_peak.right, target_peak.highest_y))
    return false;

  save_current_fit("Added peak from residuals");

  if (region_.dirty())
    rebuild(optimizer);
  return true;
}

//Peak ROI::peak(double peakID) const
//{
//  if (contains(peakID))
//    return peaks_.at(peakID);
//  else
//    return Peak();
//}

bool RegionManager::overlaps(double bin) const
{
  if (!width())
    return false;
  return ((left_bin() <= bin) && (bin <= right_bin()));
}

bool RegionManager::overlaps(double Lbin, double Rbin) const
{
  if (fit_eval_.x_.empty())
    return false;
  if (overlaps(Lbin) || overlaps(Rbin))
    return true;
  return ((Lbin <= left_bin()) && (right_bin() <= Rbin));
}

bool RegionManager::overlaps(const RegionManager& other) const
{
  if (other.width() <= 0)
    return false;
  return overlaps(other.left_bin(), other.right_bin());
}

size_t RegionManager::peak_count() const
{
  return region_.peaks_.size();
}

bool RegionManager::contains(double peakID) const
{
  return (region_.peaks_.count(peakID) > 0);
}

Peak RegionManager::peak(double peakID) const
{
  if (contains(peakID))
    return region_.peaks_.at(peakID);
  else
    return Peak();
}

const std::map<double, Peak> &RegionManager::peaks() const
{
  return region_.peaks_;
}

bool RegionManager::adjust_sum4(double peakID, double left, double right)
{
  if (region_.adjust_sum4(peakID, left, right))
  {
    render();
    save_current_fit("SUM4 adjusted on " + std::to_string(peakID));
    return true;
  }
  return false;
}

bool RegionManager::replace_hypermet(double &peakID, Peak hyp)
{
  if (region_.replace_hypermet(peakID, hyp))
  {
    peakID = hyp.id();
    render();
    save_current_fit("Hypermet adjusted on " + std::to_string(hyp.id()));
    return true;
  }
  return false;
}

//bool ROI::override_energy(double peakID, double energy)
//{
//  if (!peaks_.count(peakID))
//    return false;
//
//   peaks_[peakID].override_energy(energy);
//
//   render();
//   save_current_fit("Peak energy override " + std::to_string(peaks_.at(peakID).center().value())
//                    + "->" + std::to_string(peaks_.at(peakID).energy().value()));
//   return true;
//}

bool RegionManager::add_peak(const FitEvaluation &parentfinder,
                   double left, double right,
                   AbstractOptimizer* optimizer)
{
  bool added = false;
  if (overlaps(left) && overlaps(right))
  {
    NaiveKON kon(fit_eval_.x_, fit_eval_.y_resid_,
                 settings_.kon_settings.width,
                 settings_.kon_settings.sigma_resid);

    if (region_.add_peak(left, right, kon.highest_residual(left, right)))
    {
      save_current_fit("Manually added peak");
      render();
//      if (region_.dirty())
//        rebuild(optimizer);
      return true;
    }
  }
  else if (width()) //changing region bounds
  {
    double L  = std::min(left, left_bin());
    double R = std::max(right, right_bin());
    fit_eval_.cloneRange(parentfinder, L, R);
    render();
    save_current_fit("Implicitly expanded region");

    NaiveKON kon(fit_eval_.x_, fit_eval_.y_resid_,
                 settings_.kon_settings.width,
                 settings_.kon_settings.sigma_resid);

    if (L < region_.left())
    {
      SUM4Edge edge(fit_eval_.weighted_data.left(settings_.background_edge_samples));
      region_.replace_data(fit_eval_.weighted_data.subset(L, region_.right()), edge, RB());
    }

    if (R > region_.right())
    {
      SUM4Edge edge(fit_eval_.weighted_data.right(settings_.background_edge_samples));
      region_.replace_data(fit_eval_.weighted_data.subset(region_.left(), R), LB(), edge);
    }

    if (region_.add_peak(std::max(left, region_.left()),
                         std::min(right, region_.right()),
                         kon.highest_residual(left, right)))
    {
      save_current_fit("Manually added peak");
      render();
//      if (region_.dirty())
//        rebuild(optimizer);
      return true;
    }
  }

  return find_and_fit(optimizer);
}

bool RegionManager::remove_peaks(const std::set<double> &pks, AbstractOptimizer* optimizer)
{
  bool found = region_.remove_peaks(pks);

  if (!found)
    return false;

  // \todo save and refit if dirty
  if (region_.peaks_.size() && !rebuild(optimizer))
    return false;

  render();
  save_current_fit("Peaks removed");
  return true;
}

bool RegionManager::override_settings(const FitSettings &fs)
{
  settings_ = fs;
  settings_.overriden = true; //do this in fitter if different?

  for (auto& p : region_.peaks_)
    p.second.force_defaults(settings_.default_peak);

  save_current_fit("Region fit settings override");

  //propagate to peaks

  //render if calibs changed?
  return true;
}

void RegionManager::save_current_fit(std::string description)
{
  Fit thisfit(region_, description);
  fits_.push_back(thisfit);
  current_fit_ = fits_.size() - 1;
}

bool RegionManager::rebuild(AbstractOptimizer* optimizer)
{
  // \todo check for convergence before saving?
  if (region_.peaks_.empty())
    return false;

  region_.update_indices();
  INFO("Will rebuild\n{}", region_.to_string(" "));

  std::stringstream ss;
  auto vars = region_.variables();
  ss << vars.transpose();
  INFO("Initial vars: {}  chisq={}  chisq(at)={}", ss.str(),
      region_.chi_sq(), region_.chi_sq(vars));

  auto result = optimizer->minimize(&region_);
  //if (optimizer->verbose)
  INFO("Fitter result {}", result.to_string());
  region_.save_fit(result);
  region_.reindex_peaks();

  if (!region_.sane())
    INFO("!!!!! RESULT NOT SANE !!!!!");

  region_.auto_sum4();
  save_current_fit("Refit");
  INFO("Rebuilt as\n{}", region_.to_string(" "));

  render();
  return true;
}

void RegionManager::render()
{
  rendering_.render(region_, settings_.calib.cali_nrg_);

  std::vector<double> lowres_backsteps, lowres_fullfit;
  region_.background.eval_add(fit_eval_.x_, lowres_backsteps);
  region_.background.eval_add(fit_eval_.x_, lowres_fullfit);
  for (auto& p : region_.peaks_)
  {
    for (size_t i = 0; i < fit_eval_.x_.size(); ++i)
    {
      auto vals = p.second.eval(fit_eval_.x_[i]);
      lowres_backsteps[i] += vals.step_tail();
      lowres_fullfit[i] += vals.all();
    }
  }

  fit_eval_.update_fit(lowres_fullfit, lowres_backsteps);
}

// \todo belongs in another class?
std::vector<double> RegionManager::remove_background()
{
  std::vector<double> y_nobkg(fit_eval_.x_.size());
  for (size_t i = 0; i < fit_eval_.y_.size(); ++i)
    y_nobkg[i] = fit_eval_.y_[i] - region_.background.eval(fit_eval_.x_[i]);
  return y_nobkg;
}

bool RegionManager::adjust_LB(const FitEvaluation &parentfinder, double left, double right,
                     AbstractOptimizer* optimizer)
{
  SUM4Edge edge(parentfinder.weighted_data.subset(left, right));
  if (!edge.width() || (edge.right() >= region_.RB_.left()))
    return false;

  if (edge.left() == left_bin())
    region_.adjust_LB(edge);
  else
  {
    fit_eval_.cloneRange(parentfinder, left, right_bin());
    region_.replace_data(parentfinder.weighted_data.subset(left, right_bin()),
                         edge, RB());
  }

  save_current_fit("Left baseline adjusted");

//  if (region_.dirty())
//    rebuild(optimizer);

  render();
  return true;
}

bool RegionManager::adjust_RB(const FitEvaluation &parentfinder, double left, double right,
                    AbstractOptimizer* optimizer) {
  SUM4Edge edge(parentfinder.weighted_data.subset(left, right));
  if (!edge.width() || (edge.left() <= region_.LB_.right()))
    return false;

  if (edge.right() == right_bin())
    region_.adjust_RB(edge);
  else
  {
    fit_eval_.cloneRange(parentfinder, left_bin(), right);
    region_.replace_data(parentfinder.weighted_data.subset(left_bin(), right),
                         LB(), edge);
  }

  save_current_fit("Right baseline adjusted");

//  if (region_.dirty())
//    rebuild(optimizer);

  render();
  return true;
}

size_t RegionManager::current_fit() const
{
  return current_fit_;
}

std::vector<FitDescription> RegionManager::history() const
{
  std::vector<FitDescription> ret;
  for (auto &f : fits_)
    ret.push_back(f.description);
  return ret;
}


bool RegionManager::rollback(const FitEvaluation &/*parent_finder*/, size_t i)
{
  if (i >= fits_.size())
    return false;

  //settings_ = fits_[i].settings_;
  current_fit_ = i;
  region_ = fits_[current_fit_].region;
  // \todo reapply spectrum data
  render();

  return true;
}

nlohmann::json RegionManager::to_json(const FitEvaluation &parent_finder) const
{
  nlohmann::json j;

  if (fits_.empty())
    return j;

  j["current_fit"] = current_fit_;

  RegionManager temp(*this);

  for (size_t i=0; i < temp.fits_.size(); ++i)
  {
    nlohmann::json jj;

    jj["description"] = temp.fits_[i].description.description;
    temp.rollback(parent_finder, i);

    if (settings_.overriden)
      jj["settings"] = settings_;

    jj["background_left"] = temp.LB();
    jj["background_right"] = temp.RB();
    jj["background_poly"] = temp.region_.background;

    for (auto &p : temp.region_.peaks_)
      jj["peaks"].push_back(p.second);

    j["fits"].push_back(jj);
  }

  return j;
}

RegionManager::RegionManager(const nlohmann::json& j, const FitEvaluation &finder, const FitSettings& fs)
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

      region_.LB_ = LB;
      region_.RB_ = RB;

      if (it.value().count("peaks"))
      {
        nlohmann::json p = it.value()["peaks"];
        for (nlohmann::json::iterator it2 = p.begin(); it2 != p.end(); ++it2)
        {
          Peak newpeak = it2.value();
          newpeak.sum4 = SUM4(
              finder.weighted_data.subset(newpeak.sum4.left(), newpeak.sum4.right()),
              region_.LB_, region_.RB_);
          region_.peaks_[newpeak.id()] = newpeak;
        }
      }

      region_.background = it.value()["background_poly"];
      render();
      save_current_fit(it.value()["description"]);
      region_.peaks_.clear();
    }
  }

  rollback(finder, j["current_fit"]);
}

}
