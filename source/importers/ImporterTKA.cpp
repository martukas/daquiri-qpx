#include "ImporterTKA.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <core/consumer_factory.h>
#include <importers/string_to_chans.h>

bool ImporterTKA::validate(const boost::filesystem::path& path) const
{
  return true;
}

void ImporterTKA::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  using namespace boost::algorithm;
  std::ifstream myfile(path.string(), std::ios::in);
  std::string data;
  double timed;

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterTKA could not get a valid Histogram 1D from factory");

  myfile >> data;
  DAQuiri::Setting live_time(DAQuiri::SettingMeta("live_time", DAQuiri::SettingType::duration));
  timed = boost::lexical_cast<double>(trim_copy(data)) * 1000.0;
  live_time.set_duration(std::chrono::milliseconds(static_cast<int64_t>(timed)));
  hist->set_attribute(live_time);

  myfile >> data;
  DAQuiri::Setting real_time(DAQuiri::SettingMeta("real_time", DAQuiri::SettingType::duration));
  timed = boost::lexical_cast<double>(trim_copy(data)) * 1000.0;
  real_time.set_duration(std::chrono::milliseconds(static_cast<int64_t>(timed)));
  hist->set_attribute(real_time);

  entry_list = string_to_chans(myfile);

//  metadata_.detectors.resize(1);
//  init_from_file(name);

  hist->import(*this);

  project->add_consumer(hist);
}
