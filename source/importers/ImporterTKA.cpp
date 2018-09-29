#include "ImporterTKA.h"
#include <boost/algorithm/string.hpp>
#include <core/consumer_factory.h>
#include <importers/string_to_chans.h>

#include <core/util/custom_logger.h>

bool ImporterTKA::validate(const boost::filesystem::path& path) const
{
  return true;
}

void ImporterTKA::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  std::ifstream myfile(path.string(), std::ios::in);
  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterTKA could not get a valid Histogram 1D from factory");

  using namespace boost::algorithm;
  std::string data;

  // \todo extract ms part separately

  myfile >> data;
  auto lt_ms = static_cast<int64_t>(std::stod(trim_copy(data)) * 1000.0);
  //DBG("parsing live time = '{}' -> {}", data, lt_ms);
  hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));

  myfile >> data;
  auto rt_ms = static_cast<int64_t>(std::stod(trim_copy(data)) * 1000.0);
  //DBG("parsing real time = '{}' -> {}", data, rt_ms);
  hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));

  entry_list = string_to_chans(myfile);

//  metadata_.detectors.resize(1);
//  init_from_file(name);

  hist->set_attribute(DAQuiri::Setting::text("name", path.stem().string()));

  hist->import(*this);

  project->add_consumer(hist);
}
