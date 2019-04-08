#include <core/gamma/region_manager.h>
#include <core/util/more_math.h>
#include <core/gamma/finders/finder_kon_naive.h>

#include <core/util/logger.h>

namespace DAQuiri
{

Fit::Fit(const Region& r, std::string descr)
    : region(r)
{
  description.description = descr;
  description.peaknum = region.peaks_.size();
  if (!region.peaks_.empty())
  {
    description.chi_sq_norm = region.chi_sq() / region.degrees_of_freedom();
    UncertainDouble tot_gross{0.0, 0.0};
    UncertainDouble tot_back{0.0, 0.0};
    for (const auto& p : region.peaks_)
    {
      tot_gross += p.second.sum4.gross_area();
      tot_back += p.second.sum4.background_area();
    }
    UncertainDouble tot_net = tot_gross - tot_back;
    description.sum4aggregate = tot_net.error();
  }
}

RegionManager::RegionManager(const Region& initial_region)
{
  modify_region(initial_region, "Region created");
}

double RegionManager::id() const
{
  return region_.left();
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


const Region& RegionManager::region() const
{
  return region_;
}

void RegionManager::modify_region(const Region& new_region, std::string message)
{
  if (message.empty())
    message = "User modified region";
  region_ = new_region;
  save_current_fit(message);
}

size_t RegionManager::current_fit() const
{
  return current_fit_;
}

std::vector<FitDescription> RegionManager::history() const
{
  std::vector<FitDescription> ret;
  for (auto& f : fits_)
    ret.push_back(f.description);
  return ret;
}

bool RegionManager::rollback(size_t i)
{
  if (i >= fits_.size())
    return false;

  current_fit_ = i;
  region_ = fits_[current_fit_].region;
  // \todo reapply spectrum data
  return true;
}


void RegionManager::save_current_fit(std::string description)
{
  Fit thisfit(region_, description);
  fits_.push_back(thisfit);
  current_fit_ = fits_.size() - 1;
}

bool RegionManager::refit(AbstractOptimizer* optimizer)
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

  return true;
}

FitEvaluation RegionManager::eval() const
{
  FitEvaluation ret(region_.data);

  std::vector<double> lowres_backsteps, lowres_fullfit;
  region_.background.eval_add(region_.data.chan, lowres_backsteps);
  region_.background.eval_add(region_.data.chan, lowres_fullfit);
  for (auto& p : region_.peaks_)
  {
    for (size_t i = 0; i < region_.data.chan.size(); ++i)
    {
      auto vals = p.second.eval(region_.data.chan[i]);
      lowres_backsteps[i] += vals.step_tail();
      lowres_fullfit[i] += vals.all();
    }
  }

  ret.update_fit(lowres_fullfit, lowres_backsteps);
  return ret;
}


//bool RegionManager::add_from_resid(AbstractOptimizer* optimizer)
//{
//  if (fit_eval_.empty())
//    return false;
//
//  NaiveKON kon(fit_eval_.x_, fit_eval_.y_resid_,
//               settings_.kon_settings.width,
//               settings_.kon_settings.sigma_resid);
//
//  DetectedPeak target_peak = kon.tallest_detected();
//
//  if (target_peak.highest_y == 0.0)
//    return false;
//
//  if (!region_.add_peak(target_peak.left, target_peak.right, target_peak.highest_y))
//    return false;
//
//  save_current_fit("Added peak from residuals");
//
//  if (region_.dirty())
//    rebuild(optimizer);
//  return true;
//}
//
//
//void RegionManager::iterative_fit(AbstractOptimizer* optimizer)
//{
//  if (!settings_.calib.cali_fwhm_.valid() || region_.peaks_.empty())
//    return;
//
//  //double prev_chi_sq = peaks_.begin()->second.hypermet().chi2();
//
//  for (int i = 0; i < settings_.resid_max_iterations; ++i)
//  {
//    DBG("Attempting add from resid with {} peaks", region_.peaks_.size());
//
//    if (!add_from_resid(optimizer))
//    {
//      //      DBG << "    failed add from resid";
//      break;
//    }
//
//    if (optimizer.cancel.load())
//      break;
//  }
//}

nlohmann::json RegionManager::to_json(const FitEvaluation& parent_finder) const
{
  nlohmann::json j;

  if (fits_.empty())
    return j;

  j["current_fit"] = current_fit_;

  RegionManager temp(*this);

  for (size_t i = 0; i < temp.fits_.size(); ++i)
  {
    nlohmann::json jj;

    jj["description"] = temp.fits_[i].description.description;
    temp.rollback(i);

    jj["background_left"] = temp.region_.LB_;
    jj["background_right"] = temp.region_.RB_;
    jj["background_poly"] = temp.region_.background;

    for (auto& p : temp.region_.peaks_)
      jj["peaks"].push_back(p.second);

    j["fits"].push_back(jj);
  }

  return j;
}

RegionManager::RegionManager(const nlohmann::json& j, const WeightedData& super_region)
{
  if (!super_region.valid())
    return;

  if (j.count("fits"))
  {
    nlohmann::json o = j["fits"];
    for (nlohmann::json::iterator it = o.begin(); it != o.end(); ++it)
    {
      SUM4Edge LB = it.value()["background_left"];
      SUM4Edge RB = it.value()["background_right"];

      LB = SUM4Edge(super_region.subset(LB.left(), LB.right()));
      RB = SUM4Edge(super_region.subset(RB.left(), RB.right()));

      if (!LB.width() || !RB.width())
        return;

      //validate background and edges?
      region_ = Region(super_region.subset(LB.left(), RB.left()), 1);
      region_.LB_ = LB;
      region_.RB_ = RB;

      if (it.value().count("peaks"))
      {
        nlohmann::json p = it.value()["peaks"];
        for (nlohmann::json::iterator it2 = p.begin(); it2 != p.end(); ++it2)
        {
          Peak newpeak = it2.value();
          newpeak.sum4 = SUM4(
              super_region.subset(newpeak.sum4.left(), newpeak.sum4.right()),
              region_.LB_, region_.RB_);
          region_.peaks_[newpeak.id()] = newpeak;
        }
      }

      region_.background = it.value()["background_poly"];
      save_current_fit(it.value()["description"]);
      region_.peaks_.clear();
    }
  }

  rollback(j["current_fit"]);
}

}
