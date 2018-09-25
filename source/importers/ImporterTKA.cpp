#include "ImporterTKA.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <core/consumer_factory.h>

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

  channels_from_string(*hist, myfile, false);

//  metadata_.detectors.resize(1);
//  init_from_file(name);

  hist->import(*this);

  project->add_consumer(hist);
}

bool ImporterTKA::channels_from_string(DAQuiri::Consumer& consumer,
    std::istream& data_stream, bool compression)
{
  int i = 0;

  std::string numero, numero_z;
  if (compression)
  {
    while (data_stream.rdbuf()->in_avail())
    {
      data_stream >> numero;
      if (numero == "0")
      {
        data_stream >> numero_z;
        i += boost::lexical_cast<uint16_t>(numero_z);
      }
      else
      {
        DAQuiri::Entry new_entry;
        new_entry.first.resize(1);
        new_entry.first[0] = i;
        PreciseFloat nr{0};
        try { nr = std::stold(numero); }
        catch (...) {}
        new_entry.second = nr;
        entry_list.push_back(new_entry);
        i++;
      }
    }
  }
  else
  {
    while (data_stream.rdbuf()->in_avail())
    {
      data_stream >> numero;
      DAQuiri::Entry new_entry;
      new_entry.first.resize(1);
      new_entry.first[0] = i;
      PreciseFloat nr{0};
      try { nr = std::stold(numero); }
      catch (...) {}
      new_entry.second = nr;
      entry_list.push_back(new_entry);
      i++;
    }
  }

  if (i == 0)
    return false;

//  bits_ = log2(i);
//  if (pow(2, bits_) < i)
//    bits_++;
//  maxchan_ = i;
//
//  Setting res = metadata_.get_attribute("resolution");
//  res.value_int = bits_;
//  metadata_.set_attribute(res);

  return true;
}
