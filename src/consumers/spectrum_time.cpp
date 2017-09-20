#include "spectrum_time.h"
#include "sparse_map2d.h"

#include "custom_logger.h"

#define kDimensions 2

TimeSpectrum::TimeSpectrum()
{
  data_ = std::make_shared<SparseMap2D>(); //use dense 2d

  Setting base_options = metadata_.attributes();
  metadata_ = ConsumerMetadata(my_type(), "Spectra in time series");

  SettingMeta ds("downsample", SettingType::integer);
  ds.set_val("units", "bits");
  ds.set_flag("preset");
  ds.set_val("min", 0);
  ds.set_val("max", 31);
  base_options.branches.add(ds);

  SettingMeta cutoff_bin("cutoff", SettingType::integer);
  cutoff_bin.set_val("description", "Hits rejected below minimum value");
  cutoff_bin.set_val("min", 0);
  cutoff_bin.set_flag("preset");
  base_options.branches.add(cutoff_bin);

  SettingMeta val_name("value_name", SettingType::text);
  val_name.set_flag("preset");
  val_name.set_val("description", "Name of event value to bin");
  base_options.branches.add(val_name);

  SettingMeta add_channels("add_channels", SettingType::pattern, "Channels to bin");
  add_channels.set_flag("preset");
  add_channels.set_val("chans", 1);
  base_options.branches.add(add_channels);

  metadata_.overwrite_all_attributes(base_options);
}

bool TimeSpectrum::_initialize()
{
  Spectrum::_initialize();

  downsample_ = metadata_.get_attribute("downsample").get_number();
  cutoff_ = metadata_.get_attribute("cutoff").get_number();
  val_name_ = metadata_.get_attribute("value_name").get_text();
  channels_ = metadata_.get_attribute("add_channels").pattern();

  int adds = 1;//0;
  //  std::vector<bool> gts = add_channels_.gates();
  //  for (size_t i=0; i < gts.size(); ++i)
  //    if (gts[i])
  //      adds++;

  if (adds != 1)
  {
    WARN << "<TimeSpectrum> Cannot initialize. Add pattern must have 1 selected channel.";
    return false;
  }

  this->_recalc_axes();
  return true;
}

void TimeSpectrum::_init_from_file()
{
  metadata_.set_attribute(Setting::integer("downsample", downsample_));

  channels_.resize(1);
  channels_.set_gates(std::vector<bool>({true}));

  metadata_.set_attribute(Setting("add_channels", channels_));

  Spectrum::_init_from_file();
}

void TimeSpectrum::_set_detectors(const std::vector<Detector>& dets)
{
  metadata_.detectors.resize(kDimensions, Detector());

  if (dets.size() == kDimensions)
    metadata_.detectors = dets;

  if (dets.size() >= kDimensions)
  {
    for (size_t i=0; i < dets.size(); ++i)
    {
      if (metadata_.chan_relevant(i))
      {
        metadata_.detectors[0] = dets[i];
        break;
      }
    }
  }

  this->_recalc_axes();
}

void TimeSpectrum::_recalc_axes()
{
  data_->set_axis(0, DataAxis(Calibration(), 0));
  data_->set_axis(1, DataAxis(Calibration(), 0));

  if (data_->dimensions() != metadata_.detectors.size())
    return;

  for (size_t i=0; i < metadata_.detectors.size(); ++i)
  {
    auto det = metadata_.detectors[i];
    std::string valname = (i == 0) ? "v1" : "v2";
    CalibID from(det.id(), valname, "", 0);
    CalibID to("", valname, "", 0);
    auto calib = det.get_preferred_calibration(from, to);
    data_->set_axis(i, DataAxis(calib, 0));
  }

  data_->recalc_axes(0);
}

bool TimeSpectrum::channel_relevant(int16_t channel) const
{
  return ((channel >= 0) && channels_.relevant(channel));
}

void TimeSpectrum::_push_event(const Event& e)
{
  const auto& c = e.channel();

  if (!this->channel_relevant(c) ||
      !value_relevant(c, value_idx_))
    return;

  if (downsample_)
  {
    coords_[0] = (e.value(value_idx_[c]) >> downsample_);
    coords_[1] = (e.value(value_idx_[c]) >> downsample_);
  }
  else
  {
    coords_[0] = e.value(value_idx_[c]);
    coords_[1] = e.value(value_idx_[c]);
  }

  data_->add_one(coords_);
  total_count_++;
  recent_count_++;
}


void TimeSpectrum::_push_stats(const Status& newBlock)
{
  if (!this->channel_relevant(newBlock.channel()))
    return;

  if (newBlock.channel() >= static_cast<int16_t>(value_idx_.size()))
    value_idx_.resize(newBlock.channel() + 1, -1);
  if (newBlock.event_model().name_to_val.count(val_name_))
    value_idx_[newBlock.channel()] = newBlock.event_model().name_to_val.at(val_name_);

  PreciseFloat real = 0;
  PreciseFloat live = 0;
  PreciseFloat percent_dead = 0;
  PreciseFloat tot_time = 0;

//  spectra_.push_back(std::vector<PreciseFloat>(pow(2, bits_)));

  if (!updates_.empty())
  {
    if ((newBlock.type() == StatusType::stop) &&
        (updates_.back().type() == StatusType::running))
    {
      updates_.pop_back();
      seconds_.pop_back();
    }

    boost::posix_time::time_duration rt ,lt;
    lt = rt = newBlock.time() - updates_.back().time();

    Status diff;// = newBlock - updates_.back();
    PreciseFloat scale_factor = 1;
    if (diff.stats().count("native_time") && (diff.stats().at("native_time") > 0))
      scale_factor = rt.total_microseconds() / diff.stats()["native_time"];

    if (diff.stats().count("live_time"))
    {
      PreciseFloat scaled_live = diff.stats().at("live_time") * scale_factor;
      lt = boost::posix_time::microseconds(static_cast<long>(to_double(scaled_live)));
    }

    real     = rt.total_milliseconds()  * 0.001;
    tot_time = (newBlock.time() - updates_.front().time()).total_milliseconds() * 0.001;

    live = lt.total_milliseconds() * 0.001;

    percent_dead = (real - live) / real * 100;
  }


  if (seconds_.empty() || (tot_time != 0))
  {
    seconds_.push_back(to_double(tot_time));
    updates_.push_back(newBlock);

    //      spectrum_.push_back(count);

//    axes_[0].clear();
//    for (auto &q : seconds_)
//      axes_[0].push_back(to_double(q));
  }

  Spectrum::_push_stats(newBlock);
}


