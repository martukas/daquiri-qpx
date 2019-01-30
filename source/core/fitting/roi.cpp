#include <core/fitting/roi.h>
#include <core/util/custom_logger.h>
#include <core/util/timer.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri {

Fit::Fit(const Region& r, std::string descr)
    : region(r)
{
  description.description = descr;
  description.peaknum = region.peaks_.size();
  if (!region.peaks_.empty())
  {
    //description.chi_sq_norm = chi_sq_normalized();
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

void PeakRendering::render(const Peak& h)
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

void RegionRendering::render(const Region& r,
                             const Calibration& energy_calib)
{
  Polynomial sum4back = SUM4Edge::sum4_background(r.LB_, r.RB_);

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
    for (size_t i = 0; i < channel.size(); ++i)
    {
      auto vals = hyp.eval(channel[i]);
      peak.peak[i] += vals.peak_skews();
      peak.full_fit[i] += vals.peak_skews();
    }
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

  region_ = Region(finder_.weighted_data, settings_.background_edge_samples);
  render();
}

bool ROI::refit(BFGS& optimizer)
{
  if (region_.peaks_.empty())
    return find_and_fit(optimizer);
  else
    return rebuild(optimizer);
}

bool ROI::find_and_fit(BFGS& optimizer)
{
  finder_.y_resid_ = finder_.y_;
  finder_.find_peaks();  //assumes default params!!!
  region_ = Region(finder_.weighted_data, settings_.background_edge_samples);

  if (finder_.filtered.empty())
  {
    render();
    return false;
  }

  //std::vector<double> y_nobkg = remove_background();

  for (const auto& p : finder_.filtered)
  {
//                   subset should be from y_nobkg
//    auto subset = finder_.weighted_data.subset(finder_.x_[finder_.lefts[i]],
//                                               finder_.x_[finder_.rights[i]]);
//    auto gaussian = Peak().gaussian_only();
//    //optimizer->fit(gaussian, x_pk, y_pk);
//
//    if (gaussian.sanity_check(finder_.x_[finder().lefts[i]], finder_.x_[finder_.rights[i]]))
//    {
//      gaussian.force_defaults(region_.default_peak_);
//      region_.peaks_[gaussian.id()] = gaussian;
//    }
    region_.add_peak(p.left, p.right);
  }
  region_.reindex_peaks();
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

void ROI::iterative_fit(BFGS& optimizer)
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

    if (optimizer.cancel.load())
      break;
  }
}

bool ROI::add_from_resid(BFGS& optimizer)
{
  if (finder_.filtered.empty())
    return false;

  DetectedPeak target_peak = finder_.tallest_detected();

  if (target_peak.highest_y == 0.0)
    return false;

  if (!region_.add_peak(target_peak.left, target_peak.right, target_peak.highest_y))
    return false;

  save_current_fit("Added peak from residuals");

  if (region_.dirty)
    rebuild(optimizer);
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

size_t ROI::peak_count() const
{
  return region_.peaks_.size();
}

bool ROI::contains(double peakID) const
{
  return (region_.peaks_.count(peakID) > 0);
}

Peak ROI::peak(double peakID) const
{
  if (contains(peakID))
    return region_.peaks_.at(peakID);
  else
    return Peak();
}

const std::map<double, Peak> &ROI::peaks() const
{
  return region_.peaks_;
}

bool ROI::adjust_sum4(double peakID, double left, double right)
{
  if (region_.adjust_sum4(peakID, left, right))
  {
    render();
    save_current_fit("SUM4 adjusted on " + std::to_string(peakID));
    return true;
  }
  return false;
}

bool ROI::replace_hypermet(double &peakID, Peak hyp)
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

bool ROI::add_peak(const Finder &parentfinder,
                   double left, double right,
                   BFGS& optimizer)
{
  double center_prelim = (left+right) * 0.5; //assume down the middle

  if (overlaps(left) && overlaps(right))
  {
    if (region_.add_peak(left, right, finder_.y_resid_[finder_.find_index(center_prelim)]))
    {
      save_current_fit("Manually added peak");
      return true;
    }
  }
  else if (width()) //changing region bounds
  {
    left  = std::min(left, left_bin());
    right = std::max(right, right_bin());
    if (!finder_.cloneRange(parentfinder, left, right))
      return false;

    finder_.find_peaks();  //assumes default params!!!
    region_.replace_data(finder_.weighted_data);

    if (region_.add_peak(left, right, finder_.y_resid_[finder_.find_index(center_prelim)]))
    {
      save_current_fit("Manually added peak");
      return true;
    }
    else
      return find_and_fit(optimizer); // \todo maybe not?
  }

  DBG("<ROI> could not add peak");
  return false;
}

bool ROI::remove_peaks(const std::set<double> &pks, BFGS& optimizer)
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

bool ROI::override_settings(const FitSettings &fs)
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
  Fit thisfit(region_, description);
  fits_.push_back(thisfit);
  current_fit_ = fits_.size() - 1;
}

bool ROI::rebuild(BFGS& optimizer)
{
  if (region_.peaks_.empty())
    return false;

  region_.map_fit();
  auto result = optimizer.BFGSMin(&region_, 0.0001);
  // \todo check for convergence?
  region_.save_fit_uncerts(result);
  region_.auto_sum4();
  save_current_fit("Rebuild");

  render();
  return true;
}

void ROI::render()
{
  rendering_.render(region_, settings_.calib.cali_nrg_);

  std::vector<double> lowres_backsteps, lowres_fullfit;
  region_.background.eval_add(finder_.x_, lowres_backsteps);
  region_.background.eval_add(finder_.x_, lowres_fullfit);
  for (auto& p : region_.peaks_)
  {
    for (size_t i = 0; i < finder_.x_.size(); ++i)
    {
      auto vals = p.second.eval(finder_.x_[i]);
      lowres_backsteps[i] += vals.step_tail();
      lowres_fullfit[i] += vals.all();
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
    y_nobkg[i] = finder_.y_[i] - region_.background.eval(finder_.x_[i]);
  return y_nobkg;
}

bool ROI::adjust_LB(const Finder &parentfinder, double left, double right,
                     BFGS& optimizer)
{
  SUM4Edge edge(parentfinder.weighted_data.subset(left, right));
  if (!edge.width() || (edge.right() >= region_.RB_.left()))
    return false;

  if ((edge.left() != left_bin()) && !finder_.cloneRange(parentfinder, left, right_bin()))
    return false;

  region_.replace_data(parentfinder.weighted_data.subset(left, right_bin()),
      edge, RB());
  save_current_fit("Left baseline adjusted");

  if (region_.dirty)
    rebuild(optimizer);

  render();
  return true;
}

bool ROI::adjust_RB(const Finder &parentfinder, double left, double right,
                    BFGS& optimizer) {
  SUM4Edge edge(parentfinder.weighted_data.subset(left, right));
  if (!edge.width() || (edge.left() <= region_.LB_.right()))
    return false;

  if ((edge.right() != right_bin()) && !finder_.cloneRange(parentfinder, left_bin(), right))
    return false;

  region_.replace_data(parentfinder.weighted_data.subset(left_bin(), right),
                       LB(), edge);

  save_current_fit("Right baseline adjusted");

  if (region_.dirty)
    rebuild(optimizer);

  render();
  return true;
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
  current_fit_ = i;
  region_ = fits_[current_fit_].region;
  // \todo reapply spectrum data
  render();

  return true;
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
    jj["background_poly"] = temp.region_.background;

    for (auto &p : temp.region_.peaks_)
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
