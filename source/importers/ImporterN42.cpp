#include "ImporterN42.h"
#include <boost/algorithm/string.hpp>
#include <core/consumer_factory.h>
#include <importers/string_to_chans.h>
#include <pugixml.hpp>

#include <core/util/custom_logger.h>

bool ImporterN42::validate(const boost::filesystem::path& path) const
{
  return true;
}

void ImporterN42::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  pugi::xml_document doc;

  if (!doc.load_file(path.string().c_str()))
    throw std::runtime_error("ImporterN42 could not load XML document");

  pugi::xml_node root = doc.first_child();
  if (!root || (std::string(root.name()) != "N42InstrumentData"))
    throw std::runtime_error("ImporterN42 could not find N42InstrumentData node");

  pugi::xml_node meas_node = root.child("Measurement");
  if (!meas_node)
    throw std::runtime_error("ImporterN42 could not find Measurement node");

  pugi::xml_node node = meas_node.child("Spectrum");
  if (!node)
    throw std::runtime_error("ImporterN42 could not find Spectrum node");

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterN42 could not get a valid Histogram 1D from factory");

  std::string data;

  data = node.child_value("RealTime");
  if (!data.empty()) {
    boost::algorithm::trim(data);
    if (data.size() > 3)
      data = data.substr(2, data.size()-3); //to trim PTnnnS to nnn
    auto rt_ms = static_cast<int64_t>(std::stod(data) * 1000.0);
    hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));
  }

  data = node.child_value("LiveTime");
  if (!data.empty()) {
    boost::algorithm::trim(data);
    if (data.size() > 3)
      data = data.substr(2, data.size()-3); //to trim PTnnnS to nnn
    auto lt_ms = static_cast<int64_t>(std::stod(data) * 1000.0);
    hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));
  }

  data = node.child_value("StartTime");
  if (!data.empty()) {
    hist->set_attribute(DAQuiri::Setting("start_time", from_iso_extended(data)));
  }

  data = node.child_value("ChannelData");
  boost::algorithm::trim(data);
  std::stringstream channel_data(data);

  entry_list = string_to_chans_zero_suppressed(channel_data);

//  Detector newdet;
//  if (node.attribute("Detector"))
//    newdet.set_name(node.attribute("Detector").value());
//  if (node.child_value("DetectorType"))
//    newdet.set_type(node.child_value("DetectorType"));
//
//  //have more claibrations? fwhm,etc..?
//  if (node.child("Calibration"))
//  {
//    Calibration newcalib(bits_);
//    newcalib.from_xml(node.child("Calibration"));
//    newdet.set_energy_calibration(newcalib);
//  }
//
//  metadata_.detectors.resize(1);
//  metadata_.detectors[0] = newdet;

  hist->set_attribute(DAQuiri::Setting::text("name", path.stem().string()));
  hist->set_attribute(DAQuiri::Setting::boolean("visible", true));

  hist->import(*this);

  project->add_consumer(hist);
}
