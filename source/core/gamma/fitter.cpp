#include <core/gamma/fitter.h>
#include <core/gamma/finders/finder_kon_naive.h>
#include <core/fitting/data_model/weight_strategies.h>
#include <core/util/time_extensions.h>
#include <algorithm>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void Fitter::setData(ConsumerPtr spectrum)
{
//  clear();
  if (spectrum)
  {
    ConsumerMetadata md = spectrum->metadata();
    auto data = spectrum->data();
    // \todo check downshift
    //Setting res = md.get_attribute("resolution");
    if (!data || (data->dimensions() != 1)
      //|| (res.get_int() <= 0)
        )
      return;

    metadata_ = md;

//    finder_.settings_.bits_ = res.get_int();

    if (!md.detectors.empty())
      detector_ = md.detectors[0];

    settings_.calib.cali_nrg_ = detector_.get_calibration({"energy", detector_.id()}, {"energy"});
    settings_.calib.cali_fwhm_ = detector_.get_calibration({"energy", detector_.id()}, {"fwhm"});
    settings_.live_time = md.get_attribute("live_time").duration();

    auto spectrum_dump = data->all_data();
    std::vector<double> x;
    std::vector<double> y;

    int i = 0, j = 0;
    int x_bound = 0;
    bool go = false;
    for (auto it : *spectrum_dump)
    {
      if (it.second > 0)
        go = true;
      if (go)
      {
        x.push_back(static_cast<double>(i));
        y.push_back(static_cast<double>(it.second));
        if (it.second > 0)
          x_bound = j + 1;
        j++;
      }
      i++;
    }

    x.resize(x_bound);
    y.resize(x_bound);

    WeightedData wd(x, y);
    // \todo allow variations
    wd.count_weight = weight_true(wd.count);

    fit_eval_ = FitEvaluation(wd);
    apply_settings(settings_);
  }
}

bool Fitter::empty() const
{
  return regions_.empty();
}

void Fitter::clear()
{
  detector_ = Detector();
  metadata_ = ConsumerMetadata();
  fit_eval_.clear();
  regions_.clear();
}

void Fitter::find_regions()
{
  regions_.clear();
//  DBG << "Fitter: looking for " << filtered.size()  << " peaks";

  NaiveKON kon(fit_eval_.x_, fit_eval_.y_,
               settings_.kon_settings.width,
               settings_.kon_settings.sigma_spectrum);

  if (kon.filtered.empty())
    return;

  std::vector<DetectedPeak> new_regions;

  DetectedPeak bounds = kon.filtered[0];
  for (const auto& p : kon.filtered)
  {
    double margin = settings_.kon_settings.width; // \todo use another param?
//    if (!finder_.fw_theoretical_bin.empty())
//      margin = settings_.ROI_extend_background * finder_.fw_theoretical_bin[R];

    if (p.left < (bounds.right + 2 * margin))
    {
//      DBG << "cat ROI " << L << " " << R << " " << finder_.lefts[i] << " " << finder_.rights[i];
      bounds.left = std::min(bounds.left, p.left);
      bounds.right = std::max(bounds.right, p.right);
//      DBG << "postcat ROI " << L << " " << R;
    }
    else
    {
//      DBG << "<Fitter> Creating ROI " << L << "-" << R;
      bounds.left -= margin; //if (L < 0) L = 0;
      bounds.right += margin;
      if (settings().calib.cali_nrg_.transform(bounds.right) > settings().finder_cutoff_kev)
      {
//        DBG << "<Fitter> region " << L << "-" << R;
        new_regions.push_back(bounds);
      }
      bounds = p;
    }
  }
  double margin = settings_.kon_settings.width; // \todo use another param?
//  if (!finder_.fw_theoretical_bin.empty())
//    margin = finder_.settings_.ROI_extend_background * finder_.fw_theoretical_bin[R];
  bounds.right += margin;
  new_regions.push_back(bounds);


  //extend limits of ROI to edges of neighbor ROIs (grab more background)
  if (new_regions.size() > 2)
  {
    for (size_t i = 0; (i + 1) < new_regions.size(); ++i)
    {
      if (new_regions[i].right < new_regions[i + 1].left)
      {
        double mid = (new_regions[i].right + new_regions[i + 1].left) / 2;
        new_regions[i].right = mid - 1;
        new_regions[i + 1].left = mid + 1;
      }

//      int32_t Rthis = Rs[i];
//      Rs[i] = Ls[i+1] - 1;
//      Ls[i+1] = Rthis + 1;
    }
  }

  for (const auto& r : new_regions)
  {
    RegionManager newROI(settings_, fit_eval_, r.left, r.right);
    if (newROI.width())
      regions_[newROI.id()] = newROI;
  }
//  DBG << "<Fitter> Created " << regions_.size() << " regions";

}

size_t Fitter::peak_count() const
{
  size_t tally = 0;
  for (auto& r : regions_)
    tally += r.second.peak_count();
  return tally;
}

bool Fitter::contains_peak(double bin) const
{
  for (auto& r : regions_)
    if (r.second.contains(bin))
      return true;
  return false;
}

Peak Fitter::peak(double peakID) const
{
  for (auto& r : regions_)
    if (r.second.contains(peakID))
      return r.second.peaks().at(peakID);
  return Peak();
}

size_t Fitter::region_count() const
{
  return regions_.size();
}

bool Fitter::contains_region(double bin) const
{
  return (regions_.count(bin) > 0);
}

RegionManager Fitter::region(double bin) const
{
  if (contains_region(bin))
    return regions_.at(bin);
  else
    return RegionManager();
}

const std::map<double, RegionManager>& Fitter::regions() const
{
  return regions_;
}

std::map<double, Peak> Fitter::peaks() const
{
  std::map<double, Peak> peaks;
  for (auto& q : regions_)
    if (q.second.peak_count())
    {
      std::map<double, Peak> roipeaks(q.second.peaks());
      peaks.insert(roipeaks.begin(), roipeaks.end());
    }
  return peaks;
}

std::set<double> Fitter::relevant_regions(double left, double right)
{
  auto L = std::min(left, right);
  auto R = std::max(left, right);
  std::set<double> ret;
  for (auto& r : regions_)
    if (r.second.overlaps (L, R))
      ret.insert(r.second.id());
  return ret;
}

bool Fitter::delete_ROI(double regionID)
{
  if (!contains_region(regionID))
    return false;

  auto it = regions_.find(regionID);
  regions_.erase(it);
  render_all();
  return true;
}

void Fitter::clear_all_ROIs()
{
  regions_.clear();
  render_all();
}

RegionManager Fitter::parent_region(double peakID) const
{
  for (auto& m : regions_)
    if (m.second.contains(peakID))
      return m.second;
  return RegionManager();
}

RegionManager* Fitter::parent_of(double peakID)
{
  RegionManager* parent = nullptr;
  for (auto& m : regions_)
    if (m.second.contains(peakID))
    {
      parent = &m.second;
      break;
    }
  return parent;
}

void Fitter::render_all()
{
  fit_eval_.reset();
  for (auto &r : regions_)
    fit_eval_.merge_fit(r.second.finder());
}

bool Fitter::find_and_fit(double regionID, AbstractOptimizer* optimizer)
{
  if (!regions_.count(regionID))
    return false;

  regions_[regionID].find_and_fit(optimizer);
  render_all();
  return true;
}

bool Fitter::refit_region(double regionID, AbstractOptimizer* optimizer)
{
  if (!contains_region(regionID))
    return false;

  regions_[regionID].refit(optimizer);
  render_all();
  return true;
}

double Fitter::adj_LB(double regionID, double left, double right,
                    AbstractOptimizer* optimizer)
{
  if (!contains_region(regionID))
    return -1;

  RegionManager newROI = regions_[regionID];
  if (!newROI.adjust_LB(fit_eval_, left, right, optimizer))
    return -1;
  regions_.erase(regionID);
  regions_[newROI.id()] = newROI;
  render_all();
  return newROI.id();
}

bool Fitter::adj_RB(double regionID, double left, double right,
                    AbstractOptimizer* optimizer)
{
  if (!contains_region(regionID))
    return false;

  RegionManager newROI = regions_[regionID];
  if (!newROI.adjust_RB(fit_eval_, left, right, optimizer))
    return false;
  regions_.erase(regionID);
  regions_[newROI.id()] = newROI;
  render_all();
  return true;
}

bool Fitter::override_ROI_settings(double regionID, const FitSettings& fs)
{
  if (!contains_region(regionID))
    return false;

  if (!regions_[regionID].override_settings(fs))
    return false;
  //refit?

  render_all();
  return true;

}

bool Fitter::merge_regions(double left, double right,
                           AbstractOptimizer* optimizer)
{
  std::set<double> rois = relevant_regions(left, right);
  double min = std::min(left, right);
  double max = std::max(left, right);

  for (auto& r : rois)
  {
    if (regions_.count(r) && (regions_.at(r).left_bin() < min))
      min = regions_.at(r).left_bin();
    if (regions_.count(r) && (regions_.at(r).right_bin() > max))
      max = regions_.at(r).right_bin();
    regions_.erase(r);
  }

  RegionManager newROI(settings_, fit_eval_, min, max);

  //add old peaks?
  newROI.find_and_fit(optimizer);
  if (!newROI.width())
    return false;

  regions_[newROI.id()] = newROI;
  render_all();
  return true;
}

bool Fitter::adjust_sum4(double& peak_center, double left, double right)
{
  RegionManager* parent = parent_of(peak_center);
  if (!parent)
    return false;

  return parent->adjust_sum4(peak_center, left, right);
}

bool Fitter::replace_hypermet(double& peak_center, Peak hyp)
{
  RegionManager* parent = parent_of(peak_center);
  if (!parent)
    return false;

  if (!parent->replace_hypermet(peak_center, hyp))
    return false;
  render_all();
  return true;
}

bool Fitter::rollback_ROI(double regionID, size_t point)
{
  if (!contains_region(regionID))
    return false;

  if (!regions_[regionID].rollback(fit_eval_, point))
    return false;

  render_all();
  return true;
}

double Fitter::add_peak(double left, double right, AbstractOptimizer* optimizer)
{
  if (fit_eval_.x_.empty())
    return -1;

  double id = -1;
  for (auto& q : regions_)
  {
    if (q.second.overlaps(left, right))
    {
      id = q.first;
      break;
    }
  }

  if (id != -1)
  {
    RegionManager newROI = regions_[id];
    if (newROI.add_peak(fit_eval_, left, right, optimizer))
    {
      regions_.erase(id);
      regions_[newROI.id()] = newROI;
      render_all();
      return newROI.id();
    }
    else
      return -1;
  }

//  DBG << "<Fitter> making new ROI to add peak manually " << left << " " << right;
  RegionManager newROI(settings_, fit_eval_, left, right);
//  newROI.add_peak(left, right);
  newROI.refit(optimizer);

  regions_[newROI.id()] = newROI;
  render_all();
  return newROI.id();
}

bool Fitter::remove_peaks(std::set<double> bins,
                          AbstractOptimizer* optimizer)
{
  bool changed = false;
  for (auto& m : regions_)
    if (m.second.remove_peaks(bins, optimizer))
      changed = true;
  if (changed)
    render_all();
  return changed;
}

void Fitter::apply_settings(FitSettings settings)
{
  settings_.clone(settings);
  //propagate to regions?
  // \todo reenable
//  if (regions_.empty())
//    finder_.find_peaks();
}

//bool Fitter::override_energy(double peakID, double energy)
//{
//  ROI *parent = parent_of(peakID);
//  if (!parent)
//    return false;
//
//  return parent->override_energy(peakID, energy);
//}

//void Fitter::apply_energy_calibration(Calibration cal) {
//  settings_.cali_nrg_ = cal;
//  apply_settings(settings_);
//}

//void Fitter::apply_fwhm_calibration(Calibration cal) {
//  settings_.cali_fwhm_ = cal;
//  apply_settings(settings_);
//}

std::set<double> Fitter::get_selected_peaks() const
{
  return selected_peaks_;
}

void Fitter::set_selected_peaks(std::set<double> selected_peaks)
{
  selected_peaks_ = selected_peaks;
  filter_selection();
}

void Fitter::filter_selection()
{
  std::set<double> sel;
  for (auto& p : selected_peaks_)
    if (contains_peak(p))
      sel.insert(p);
  selected_peaks_ = sel;
}

void Fitter::save_report(std::string filename)
{
  std::ofstream file(filename, std::ios::out | std::ios::app);
  file << "Spectrum \"" << metadata_.get_attribute("name").get_text() << "\"" << std::endl;
  file << "========================================================" << std::endl;
  //file << "Bits: " << finder_.settings_.bits_ << "    Resolution: " << pow(2,finder_.settings_.bits_) << std::endl;

//  file << "Match pattern:  ";
//  for (auto &q : metadata_.match_pattern)
//    file << q << " ";
//  file << std::endl;

//  file << "Add pattern:    ";
//  for (auto &q : metadata_.add_pattern)
//    file << q << " ";
//  file << std::endl;

  file << "Spectrum type: " << metadata_.type() << std::endl;

//  if (!metadata_.attributes.branches.empty()) {
//    file << "Attributes" << std::endl;
//    file << metadata_.attributes;
//  }

  if (!metadata_.detectors.empty())
  {
    file << "Detectors" << std::endl;
    for (auto& q : metadata_.detectors)
      file << "   " << q.id() << " (" << q.type() << ")" << std::endl;
  }

  file << "========================================================" << std::endl;
  file << std::endl;

  file << "Acquisition start time:  " << to_iso_extended(metadata_.get_attribute("start_time").time()) << std::endl;
  double lt = settings_.live_time.count() * 0.001;
//  double rt = finder_.settings_.real_time.total_milliseconds() * 0.001;
  file << "Live time(s):   " << lt << std::endl;
//  file << "Real time(s):   " << rt << std::endl;
//  if ((lt < rt) && (rt > 0))
//    file << "Dead time(%):   " << (rt-lt)/rt*100 << std::endl;
//  double tc = metadata_.total_count.convert_to<double>();
//  file << "Total count:    " << tc << std::endl;
//  if ((tc > 0) && (lt > 0))
//    file << "Count rate:     " << tc/lt << " cps(total/live)"<< std::endl;
//  if ((tc > 0) && (rt > 0))
//    file << "Count rate:     " << tc/rt << " cps(total/real)"<< std::endl;
//  file << std::endl;

  file.fill(' ');
  file << "========================================================" << std::endl;
  file << "================ Fitter analysis results ===============" << std::endl;
  file << "========================================================" << std::endl;

  file << std::endl;
  file.fill('-');
  file << std::setw(15) << "center(Hyp)" << "--|"
       << std::setw(15) << "energy(Hyp)" << "--|"
       << std::setw(15) << "FWHM(Hyp)" << "--|"
       << std::setw(25) << "area(Hyp)" << "--|"
       << std::setw(16) << "cps(Hyp)" << "-||"

       << std::setw(25) << "center(S4)" << "--|"
       << std::setw(15) << "cntr-err(S4)" << "--|"
       << std::setw(15) << "FWHM(S4)" << "--|"
       << std::setw(25) << "bckg-area(S4)" << "--|"
       << std::setw(15) << "bckg-err(S4)" << "--|"
       << std::setw(25) << "area(S4)" << "--|"
       << std::setw(15) << "area-err(S4)" << "--|"
       << std::setw(15) << "cps(S4)" << "--|"
       << std::setw(5) << "CQI" << "--|"
       << std::endl;
  file.fill(' ');
  for (auto& q : peaks())
  {
    file << std::setw(16) << std::setprecision(10) << std::to_string(q.second.peak_position().value()) << " | "
         << std::setw(15) << std::setprecision(10)
         << std::to_string(q.second.peak_energy(settings_.calib.cali_nrg_).value()) << " | "
         //         << std::setw( 15 ) << std::setprecision( 10 ) << q.second.fwhm_hyp << " | "
         << std::setw(26) << std::to_string(q.second.area().value()) << " | ";
//         << std::setw( 15 ) << std::setprecision( 10 ) << q.second.cps_hyp << " || "

//      file << std::setw( 26 ) << q.second.sum4_.centroid.val_uncert(10) << " | "
//           << std::setw( 15 ) << q.second.sum4_.centroid.err(10) << " | "

//           << std::setw( 15 ) << std::setprecision( 10 ) << q.second.fwhm_sum4 << " | "
//           << std::setw( 26 ) << q.second.sum4_.background_area.val_uncert(10) << " | "
//           << std::setw( 15 ) << q.second.sum4_.background_area.err(10) << " | "
//           << std::setw( 26 ) << q.second.sum4_.peak_area.val_uncert(10) << " | "
//           << std::setw( 15 ) << q.second.sum4_.peak_area.err(10) << " | "
//           << std::setw( 15 ) << std::setprecision( 10 ) << q.second.cps_sum4 << " | "
//           << std::setw( 5 ) << std::setprecision( 10 ) << q.second.good() << " |";

    file << std::endl;
  }

  file.close();
}

void to_json(json& j, const Fitter& s)
{
  j["settings"] = s.settings_;
  if (!s.selected_peaks_.empty())
    j["selected_peaks"] = s.selected_peaks_;
  for (auto& r : s.regions_)
    j["regions"].push_back(r.second.to_json(s.fit_eval_));
}

Fitter::Fitter(const json& j, ConsumerPtr spectrum)
{
  if (!spectrum)
    return;

  if (j.count("selected_peaks"))
  {
    auto o = j["selected_peaks"];
    for (json::iterator it = o.begin(); it != o.end(); ++it)
      selected_peaks_.insert(it.value().get<double>());
  }

  settings_ = j["settings"];

  setData(spectrum);

  if (j.count("regions"))
  {
    auto o = j["regions"];
    for (json::iterator it = o.begin(); it != o.end(); ++it)
    {
      RegionManager region(it.value(), fit_eval_, settings_);
      if (region.width())
        regions_[region.id()] = region;
    }
  }

  render_all();
}

}
