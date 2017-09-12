#include "consumer.h"
#include "h5json.h"
#include "ascii_tree.h"

#include "custom_logger.h"
#include "custom_timer.h"

namespace DAQuiri {

Consumer::Consumer()
{
  Setting attributes = metadata_.attributes();

  SettingMeta name("name", SettingType::text);
  name.set_val("description", "Short label");
  attributes.branches.add(name);

  SettingMeta vis("visible", SettingType::boolean);
  vis.set_val("description", "Plot visible");
  attributes.branches.add(vis);

  SettingMeta start_time("start_time", SettingType::time);
  start_time.set_val("description", "Start time");
  start_time.set_flag("readonly");
  attributes.branches.add(start_time);

  metadata_.overwrite_all_attributes(attributes);
}

bool Consumer::_initialize()
{
  metadata_.disable_presets();
  return false; //abstract sink indicates failure to init
}

void Consumer::_init_from_file()
{
  this->_initialize();
  this->_recalc_axes();
  this->_flush();
}

bool Consumer::from_prototype(const ConsumerMetadata& newtemplate)
{
  UNIQUE_LOCK_EVENTUALLY_ST

  if (metadata_.type() != newtemplate.type())
    return false;

//  metadata_.overwrite_all_attributes(newtemplate.attributes());
  metadata_.set_attributes(newtemplate.attributes_flat());

  metadata_.detectors.clear(); // really?

  return (this->_initialize());
//  DBG << "<Consumer::from_prototype>" << metadata_.get_attribute("name").value_text << " made with dims=" << metadata_.dimensions();
//  DBG << "from prototype " << metadata_.debug();
//  mutex_.unlock();
}

void Consumer::push_spill(const Spill& spill)
{
  UNIQUE_LOCK_EVENTUALLY_ST
  this->_push_spill(spill);
}

void Consumer::_push_spill(const Spill& spill)
{
  CustomTimer addspill_timer(true);

  if (!spill.detectors.empty())
    this->_set_detectors(spill.detectors);

  for (auto &q : spill.events)
    this->_push_event(q);

  for (auto &q : spill.stats)
    this->_push_stats(q.second);

//  DBG << "<" << metadata_.get_attribute("name").get_text() << "> added "
//      << spill.events.size() << " events in "
//      << addspill_timer.ms() << " ms at "
//      << addspill_timer.us() / double(spill.events.size()) << " us/hit";
}

void Consumer::flush()
{
  UNIQUE_LOCK_EVENTUALLY_ST
  this->_flush();
}

bool Consumer::changed() const
{
  SHARED_LOCK_ST
  return changed_;
}

void Consumer::set_detectors(const std::vector<Detector>& dets)
{
  UNIQUE_LOCK_EVENTUALLY_ST
  this->_set_detectors(dets);
  changed_ = true;
}

void Consumer::reset_changed()
{
  UNIQUE_LOCK_EVENTUALLY_ST
  changed_ = false;
}

//accessors for various properties
ConsumerMetadata Consumer::metadata() const
{
  SHARED_LOCK_ST
  return metadata_;
}

DataspacePtr Consumer::data() const
{
  SHARED_LOCK_ST
  return DataspacePtr(data_->clone());
}

std::string Consumer::type() const
{
  SHARED_LOCK_ST
  return my_type();
}

uint16_t Consumer::dimensions() const
{
  SHARED_LOCK_ST
  if (data_)
    return data_->dimensions();
  return 0;
}

std::string Consumer::debug(std::string prepend, bool verbose) const
{
  SHARED_LOCK_ST
  std::stringstream ss;
  ss << prepend << my_type();
  if (changed_)
    ss << "(changed)";
  ss << "\n";
  ss << prepend << k_branch_mid_B
     << metadata_.debug(prepend + k_branch_pre_B, verbose);
  if (data_)
    ss << data_->debug(prepend + k_branch_end_B);
  else
    ss << prepend << k_branch_end_B << "NODATA";
  return ss.str();
}

//change stuff

void Consumer::set_attribute(const Setting &setting, bool greedy)
{
  UNIQUE_LOCK_EVENTUALLY_ST
  metadata_.set_attribute(setting, greedy);
  changed_ = true;
}

void Consumer::set_attributes(const Setting &settings)
{
  UNIQUE_LOCK_EVENTUALLY_ST
  metadata_.set_attributes(settings.branches.data(), true);
  changed_ = true;
}




/////////////////////
/// Save and load ///
/////////////////////

bool Consumer::load(H5CC::Group& g, bool withdata)
{
  UNIQUE_LOCK_EVENTUALLY_ST
  if (!g.has_group("metadata"))
    return false;

  from_json(json(g.open_group("metadata")), metadata_);
//  metadata_.from_json(g.open_group("metadata"));

  bool ret = this->_initialize();

  if (ret && withdata && data_)
    data_->load(g);

  if (ret)
    this->_recalc_axes();

  return ret;
}

void Consumer::save(H5CC::Group& g) const
{
  SHARED_LOCK_ST
  g.write_attribute("type", this->my_type());

//  json j = metadata_.to_json();
  auto mdg = g.require_group("metadata");

  H5CC::from_json(json(metadata_), mdg);

  if (data_)
    data_->save(g);
}

}
