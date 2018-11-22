#include "ImporterMCA.h"
#include <core/consumer_factory.h>
#include <date/date.h>
#include <iostream>

#include <core/util/custom_logger.h>

struct __attribute__ ((packed)) TimeB
{
  int32_t time;
  int16_t milli_time;
  int16_t time_zone;
  int16_t DST_flag;
};

struct __attribute__ ((packed)) Elapsed
{
  int32_t live_time;
  int32_t true_time;
  int32_t sweeps;
  double comp;
};

struct __attribute__ ((packed)) XCal
{
  float ecal2;
  float ecal1;
  float ecal0;
  char unit_name[5];
  char unit_type;
  char cformat;
  char corder;
};

struct __attribute__ ((packed)) McaHeader16k
{
  int16_t read_type;
  int16_t mca_number;
  int16_t read_region;
  int32_t tagno;
  char spectrum_name[26];
  int16_t acquisition_mode;
  TimeB acquisition_start;
  Elapsed elapsed;
  XCal xcal;
  int16_t filler[19]; // \todo why is this higher than in specs?
  int16_t nchans;
};

bool ImporterMCA::validate(const boost::filesystem::path& path) const
{
  return true;
}

void ImporterMCA::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  std::ifstream file (path.string(), std::ios::binary);

  std::vector<char> buf(sizeof(McaHeader16k));

  INFO("Read array size: {}", buf.size());
  file.read(buf.data(), buf.size());

  McaHeader16k *header = reinterpret_cast<McaHeader16k *>(buf.data());

  INFO("Name: {}", std::string(header->spectrum_name));
  INFO("Mode: {}", header->acquisition_mode);

  INFO("XCal.units: {}", std::string(header->xcal.unit_name));
  INFO("XCal.ecal0: {}", header->xcal.ecal0);
  INFO("XCal.ecal1: {}", header->xcal.ecal1);
  INFO("XCal.ecal2: {}", header->xcal.ecal2);

  INFO("time.time: {}", header->acquisition_start.time);
  INFO("time.ms_time: {}", header->acquisition_start.milli_time);
  INFO("time.time_zone: {}", header->acquisition_start.time_zone);

  INFO("elapsed.live: {}", header->elapsed.live_time);
  INFO("elapsed.true: {}", header->elapsed.true_time);


  INFO("Nchans: {}", header->nchans);

  // \todo why is this signed?
  // \todo should value of nchans be used?
  std::vector<int32_t> spectrum(16384);
  file.read(reinterpret_cast<char*>(spectrum.data()), spectrum.size()*4);

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterMCA could not get a valid Histogram 1D from factory");

  hist->set_attribute(DAQuiri::Setting::text("description", std::string(header->xcal.unit_name)));

  auto lt_ms = static_cast<int64_t>(header->elapsed.live_time * 1000.0);
  hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));

  auto rt_ms = static_cast<int64_t>(header->elapsed.true_time * 1000.0);
  hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));

  for (size_t i=0; i < spectrum.size(); ++i)
  {
    DAQuiri::Entry new_entry;
    new_entry.first.resize(1);
    new_entry.first[0] = i;
    new_entry.second = PreciseFloat(spectrum[i]);
    entry_list.push_back(new_entry);
  }

  hist->set_attribute(DAQuiri::Setting::text("name", path.stem().string()));
  hist->set_attribute(DAQuiri::Setting::boolean("visible", true));

  hist->import(*this);

  project->add_consumer(hist);
}
