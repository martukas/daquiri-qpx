#include <core/fitting/fitter.h>
#include <algorithm>

#include <core/util/time_extensions.h>

#include <core/util/custom_logger.h>

namespace DAQuiri {

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

    finder_.settings_.cali_nrg_ = detector_.get_calibration({"energy", detector_.id()}, {"energy"});
    finder_.settings_.cali_fwhm_ = detector_.get_calibration({"energy", detector_.id()}, {"fwhm"});
    finder_.settings_.live_time = md.get_attribute("live_time").duration();

    auto spectrum_dump = data->all_data();
    std::vector<double> x;
    std::vector<double> y;

    int i = 0, j = 0;
    int x_bound = 0;
    bool go = false;
    for (auto it : *spectrum_dump) {
      if (it.second > 0)
        go = true;
      if (go) {
        x.push_back(static_cast<double>(i));
        y.push_back(static_cast<double>(it.second));
        if (it.second > 0)
          x_bound = j+1;
        j++;
      }
      i++;
    }

    x.resize(x_bound);
    y.resize(x_bound);

    finder_ = Finder(x, y, finder_.settings_);
    apply_settings(finder_.settings_);
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
  finder_.clear();
  regions_.clear();
}

void Fitter::find_regions() {
  regions_.clear();
//  DBG << "Fitter: looking for " << filtered.size()  << " peaks";

  finder_.find_peaks();

  if (finder_.filtered.empty())
    return;

  std::vector<size_t> Ls;
  std::vector<size_t> Rs;

  size_t L = finder_.lefts[0];
  size_t R = finder_.rights[0];
  for (size_t i=1; i < finder_.filtered.size(); ++i) {
    double margin = 0;
    if (!finder_.fw_theoretical_bin.empty())
      margin = finder_.settings_.ROI_extend_background * finder_.fw_theoretical_bin[R];
//    DBG << "Margin = " << margin;

    if (finder_.lefts[i] < (R + 2 * margin) ) {
//      DBG << "cat ROI " << L << " " << R << " " << finder_.lefts[i] << " " << finder_.rights[i];
      L = std::min(L, finder_.lefts[i]);
      R = std::max(R, finder_.rights[i]);
//      DBG << "postcat ROI " << L << " " << R;
    } else {
//      DBG << "<Fitter> Creating ROI " << L << "-" << R;
      L -= margin; //if (L < 0) L = 0;
      R += margin; if (R >= finder_.x_.size()) R = finder_.x_.size() - 1;
      if (settings().cali_nrg_.transform(R) > settings().finder_cutoff_kev) {
//        DBG << "<Fitter> region " << L << "-" << R;
        Ls.push_back(L);
        Rs.push_back(R);
      }
      L = finder_.lefts[i];
      R = finder_.rights[i];
    }
  }
  double margin = 0;
  if (!finder_.fw_theoretical_bin.empty())
    margin = finder_.settings_.ROI_extend_background * finder_.fw_theoretical_bin[R];
  R += margin; if (R >= finder_.x_.size()) R = finder_.x_.size() - 1;
  Ls.push_back(L);
  Rs.push_back(R);


  //extend limits of ROI to edges of neighbor ROIs (grab more background)
  if (Ls.size() > 2) {
    for (size_t i=0; (i+1) < Ls.size(); ++i) {
      if (Rs[i] < Ls[i+1]) {
        int mid = (Ls[i+1] + Rs[i]) / 2;
        Rs[i] = mid - 1;
        Ls[i+1] = mid + 1;
      }

//      int32_t Rthis = Rs[i];
//      Rs[i] = Ls[i+1] - 1;
//      Ls[i+1] = Rthis + 1;
    }
  }

  for (size_t i=0; i < Ls.size(); ++i) {
    ROI newROI(finder_, finder_.x_[Ls[i]], finder_.x_[Rs[i]]);
    if (newROI.width())
      regions_[newROI.ID()] = newROI;
  }
//  DBG << "<Fitter> Created " << regions_.size() << " regions";

}

size_t Fitter::peak_count() const
{
  size_t tally = 0;
  for (auto &r : regions_)
    tally += r.second.peak_count();
  return tally;
}

bool Fitter::contains_peak(double bin) const
{
  for (auto &r : regions_)
    if (r.second.contains(bin))
      return true;
  return false;
}

Peak Fitter::peak(double peakID) const
{
  for (auto &r : regions_)
    if (r.second.contains(peakID))
      return r.second.peaks().at(peakID);
  return Peak();
}

size_t Fitter::region_count() const
{
  return regions_.size();
}

bool Fitter::contains_region(double bin) const {
  return (regions_.count(bin) > 0);
}

ROI Fitter::region(double bin) const
{
  if (contains_region(bin))
    return regions_.at(bin);
  else
    return ROI();
}


const std::map<double, ROI> &Fitter::regions() const
{
  return regions_;
}

std::map<double, Peak> Fitter::peaks() const
{
  std::map<double, Peak> peaks;
  for (auto &q : regions_)
    if (q.second.peak_count()) {
      std::map<double, Peak> roipeaks(q.second.peaks());
      peaks.insert(roipeaks.begin(), roipeaks.end());
    }
  return peaks;
}

std::set<double> Fitter::relevant_regions(double left, double right)
{
  std::set<double> ret;
  for (auto & r : regions_)
  {
    if (
        ((left <= r.second.left_bin()) && (r.second.left_bin() <= right))
        ||
        ((left <= r.second.right_bin()) && (r.second.right_bin() <= right))
        ||
        r.second.overlaps(left)
        ||
        r.second.overlaps(right)
       )
      ret.insert(r.second.ID());
  }
  return ret;
}


bool Fitter::delete_ROI(double regionID) {
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

ROI Fitter::parent_region(double peakID) const
{
  for (auto &m : regions_)
    if (m.second.contains(peakID))
      return m.second;
  return ROI();
}

ROI *Fitter::parent_of(double peakID) {
  ROI *parent = nullptr;
  for (auto &m : regions_)
    if (m.second.contains(peakID)) {
      parent = &m.second;
      break;
    }
  return parent;
}

void Fitter::render_all()
{
  finder_.reset();
  for (auto &r : regions_)
    finder_.setFit(r.second.finder().x_,
                   r.second.finder().y_fit_,
                   r.second.finder().y_background_);
}

bool Fitter::auto_fit(double regionID, OptimizerPtr optimizer, std::atomic<bool>& interruptor)
{
  if (!regions_.count(regionID))
    return false;

  regions_[regionID].auto_fit(optimizer, interruptor);
  render_all();
  return true;
}

bool Fitter::refit_region(double regionID, OptimizerPtr optimizer, std::atomic<bool>& interruptor)
{
  if (!contains_region(regionID))
    return false;

  regions_[regionID].refit(optimizer, interruptor);
  render_all();
  return true;
}


bool Fitter::adj_LB(double regionID, double left, double right,
                    OptimizerPtr optimizer, std::atomic<bool>& interruptor)
{
  if (!contains_region(regionID))
    return false;

  ROI newROI = regions_[regionID];
  if (!newROI.adjust_LB(finder_, left, right, optimizer, interruptor))
    return false;
  regions_.erase(regionID);
  regions_[newROI.ID()] = newROI;
  render_all();
  return true;
}

bool Fitter::adj_RB(double regionID, double left, double right,
                    OptimizerPtr optimizer, std::atomic<bool>& interruptor)
{
  if (!contains_region(regionID))
    return false;

  ROI newROI = regions_[regionID];
  if (!newROI.adjust_RB(finder_, left, right, optimizer, interruptor))
    return false;
  regions_.erase(regionID);
  regions_[newROI.ID()] = newROI;
  render_all();
  return true;
}

bool Fitter::override_ROI_settings(double regionID, const FitSettings &fs, std::atomic<bool>& interruptor)
{
  if (!contains_region(regionID))
    return false;

  if (!regions_[regionID].override_settings(fs, interruptor))
    return false;
  //refit?

  render_all();
  return true;

}

bool Fitter::merge_regions(double left, double right,
                           OptimizerPtr optimizer, std::atomic<bool>& interruptor)
{
  std::set<double> rois = relevant_regions(left, right);
  double min = std::min(left, right);
  double max = std::max(left, right);

  for (auto & r : rois)
  {
    if (regions_.count(r) && (regions_.at(r).left_bin() < min))
      min = regions_.at(r).left_bin();
    if (regions_.count(r) && (regions_.at(r).right_bin() > max))
      max = regions_.at(r).right_bin();
    regions_.erase(r);
  }

  ROI newROI(finder_, min, max);

  //add old peaks?
  newROI.auto_fit(optimizer, interruptor);
  if (!newROI.width())
    return false;

  regions_[newROI.ID()] = newROI;
  render_all();
  return true;
}


bool Fitter::adjust_sum4(double &peak_center, double left, double right)
{
  ROI *parent = parent_of(peak_center);
  if (!parent)
    return false;

  return parent->adjust_sum4(peak_center, left, right);
}

bool Fitter::replace_hypermet(double &peak_center, Hypermet hyp)
{
  ROI *parent = parent_of(peak_center);
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

  if (!regions_[regionID].rollback(finder_, point))
    return false;

  render_all();
  return true;
}

bool Fitter::add_peak(double left, double right,
                      OptimizerPtr optimizer, std::atomic<bool>& interruptor)
{
  if (finder_.x_.empty())
    return false;

  for (auto &q : regions_) {
    if (q.second.overlaps(left, right)) {
      q.second.add_peak(finder_, left, right, optimizer, interruptor);
      render_all();
      return true;
    }
  }

//  DBG << "<Fitter> making new ROI to add peak manually " << left << " " << right;
  ROI newROI(finder_, left, right);
//  newROI.add_peak(finder_.x_, finder_.y_, left, right, interruptor);
  newROI.auto_fit(optimizer, interruptor);
  if (!newROI.width())
    return false;

  regions_[newROI.ID()] = newROI;
  render_all();
  return true;
}

bool Fitter::remove_peaks(std::set<double> bins,
                          OptimizerPtr optimizer, std::atomic<bool>& interruptor)
{
  bool changed = false;
  for (auto &m : regions_)
    if (m.second.remove_peaks(bins, optimizer, interruptor))
      changed = true;
  if (changed)
    render_all();
  return changed;
}

void Fitter::apply_settings(FitSettings settings) {
  finder_.settings_.clone(settings);
  //propagate to regions?
  if (regions_.empty())
    finder_.find_peaks();
}

bool Fitter::override_energy(double peakID, double energy)
{
  ROI *parent = parent_of(peakID);
  if (!parent)
    return false;

  return parent->override_energy(peakID, energy);
}

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
  for (auto &p : selected_peaks_)
    if (contains_peak(p))
      sel.insert(p);
  selected_peaks_ = sel;
}

void Fitter::save_report(std::string filename) {
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
    for (auto &q : metadata_.detectors)
      file << "   " << q.id() << " (" << q.type() << ")" << std::endl;
  }
  
  file << "========================================================" << std::endl;
  file << std::endl;

  file << "Acquisition start time:  " << to_iso_extended(metadata_.get_attribute("start_time").time()) << std::endl;
  double lt = finder_.settings_.live_time.count() * 0.001;
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
  file << "===========QPX Fitter analysis results===========" << std::endl;
  file << "========================================================" << std::endl;

  file << std::endl;
  file.fill('-');
  file << std::setw( 15 ) << "center(Hyp)" << "--|"
       << std::setw( 15 ) << "energy(Hyp)" << "--|"
       << std::setw( 15 ) << "FWHM(Hyp)" << "--|"
       << std::setw( 25 ) << "area(Hyp)" << "--|"
       << std::setw( 16 ) << "cps(Hyp)"  << "-||"
      
       << std::setw( 25 ) << "center(S4)" << "--|"
       << std::setw( 15 ) << "cntr-err(S4)" << "--|"
       << std::setw( 15 ) << "FWHM(S4)" << "--|"
       << std::setw( 25 ) << "bckg-area(S4)" << "--|"
       << std::setw( 15 ) << "bckg-err(S4)" << "--|"
       << std::setw( 25 ) << "area(S4)" << "--|"
       << std::setw( 15 ) << "area-err(S4)" << "--|"
       << std::setw( 15 ) << "cps(S4)"  << "--|"
       << std::setw( 5 ) << "CQI"  << "--|"
       << std::endl;
  file.fill(' ');
  for (auto &q : peaks()) {
    file << std::setw( 16 ) << std::setprecision( 10 ) << std::to_string(q.second.center()) << " | "
         << std::setw( 15 ) << std::setprecision( 10 ) << std::to_string(q.second.energy()) << " | "
//         << std::setw( 15 ) << std::setprecision( 10 ) << q.second.fwhm_hyp << " | "
         << std::setw( 26 ) << std::to_string(q.second.area_hyp()) << " | "
//         << std::setw( 15 ) << std::setprecision( 10 ) << q.second.cps_hyp << " || "

//      file << std::setw( 26 ) << q.second.sum4_.centroid.val_uncert(10) << " | "
//           << std::setw( 15 ) << q.second.sum4_.centroid.err(10) << " | "

//           << std::setw( 15 ) << std::setprecision( 10 ) << q.second.fwhm_sum4 << " | "
//           << std::setw( 26 ) << q.second.sum4_.background_area.val_uncert(10) << " | "
//           << std::setw( 15 ) << q.second.sum4_.background_area.err(10) << " | "
//           << std::setw( 26 ) << q.second.sum4_.peak_area.val_uncert(10) << " | "
//           << std::setw( 15 ) << q.second.sum4_.peak_area.err(10) << " | "
//           << std::setw( 15 ) << std::setprecision( 10 ) << q.second.cps_sum4 << " | "
           << std::setw( 5 ) << std::setprecision( 10 ) << q.second.good() << " |";

    file << std::endl;
  }
  
  file.close();
}

void to_json(json& j, const Fitter &s)
{
  j["settings"] = s.finder_.settings_;
  if (!s.selected_peaks_.empty())
    j["selected_peaks"] = s.selected_peaks_;
  for (auto &r : s.regions_)
    j["regions"].push_back(r.second.to_json(s.finder_));
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

  finder_.settings_ = j["settings"];

  setData(spectrum);

  if (j.count("regions"))
  {
    auto o = j["regions"];
    for (json::iterator it = o.begin(); it != o.end(); ++it)
    {
      ROI region(it.value(), finder_);
      if (region.width())
        regions_[region.ID()] = region;
    }
  }

  render_all();
}


}
