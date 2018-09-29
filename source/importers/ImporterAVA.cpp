#include "ImporterAVA.h"
#include <boost/algorithm/string.hpp>
#include <core/consumer_factory.h>
#include <importers/string_to_chans.h>
#include <pugixml.hpp>
#include <date/date.h>

#include <core/util/custom_logger.h>

bool ImporterAVA::validate(const boost::filesystem::path& path) const
{
  return true;
}

void ImporterAVA::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  pugi::xml_document doc;
  if (!doc.load_file(path.string().c_str()))
    throw std::runtime_error("ImporterAVA could not load XML document");

  pugi::xml_node root = doc.first_child();
  if (!root || (std::string(root.name()) != "mcadata"))
    throw std::runtime_error("ImporterAVA could not find mcadata node");

  pugi::xml_node node = root.child("spectrum");
  if (!node)
    throw std::runtime_error("ImporterAVA could not find Spectrum node");

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterAVA could not get a valid Histogram 1D from factory");

  std::string data;

  data = node.attribute("start").value();
  if (!data.empty()) {
    std::istringstream stream{data};
    date::sys_time<std::chrono::seconds> t;
    stream >> date::parse("%a %b %d %T GMT%z %Y", t);
    if (stream.fail())
      throw std::runtime_error("failed to parse " + data);
    hist->set_attribute(DAQuiri::Setting("start_time", t));
  }

  // \todo extract ms part separately

  data = node.attribute("elapsed_real").value();
  if (!data.empty()) {
    auto rt_ms = static_cast<int64_t>(std::stod(boost::algorithm::trim_copy(data)) * 1000.0);
    //DBG("parsing real time = '{}' -> {}", data, rt_ms);
    hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));
  }

  data = node.attribute("elapsed_live").value();
  if (!data.empty()) {
    auto lt_ms = static_cast<int64_t>(std::stod(boost::algorithm::trim_copy(data)) * 1000.0);
    //DBG("parsing live time = '{}' -> {}", data, lt_ms);
    hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));
  }

  data = root.child_value("spectrum");
  std::replace( data.begin(), data.end(), ',', ' ');
  boost::algorithm::trim(data);

  std::stringstream channel_data(data);
  entry_list = string_to_chans(channel_data);

//  Detector newdet;
//  for (auto &q : root.children("calibration_details"))
//  {
//    if ((std::string(q.attribute("type").value()) != "energy") ||
//        !q.child("model"))
//      continue;
//
//    Calibration newcalib(bits_);
//
//    std::string ctype;
//    if (std::string(q.child("model").attribute("type").value()) == "polynomial")
//      ctype = "Polynomial";
//    std::vector<double> encalib;
//    for (auto &p : q.child("model").children("coefficient"))
//    {
//      std::string coefvalstr = boost::algorithm::trim_copy(std::string(p.attribute("value").value()));
//      encalib.push_back(boost::lexical_cast<double>(coefvalstr));
//    }
//    newcalib.set_units("keV");
//    newcalib.set_function(ctype, encalib);
//    newdet.set_energy_calibration(newcalib);
//  }
//
//  metadata_.detectors.resize(1);
//  metadata_.detectors[0] = newdet;

  hist->import(*this);

  project->add_consumer(hist);
}
