#include "ImporterCNF.h"
#include <core/consumer_factory.h>
#include "xylib/xylib.h"
#include <date/date.h>
#include <core/util/string_extensions.h>
#include <core/calibration/polynomial.h>

#include <core/util/custom_logger.h>

bool ImporterCNF::validate(const boost::filesystem::path& path) const
{
  (void) path;
  return true;
}

void ImporterCNF::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  std::unique_ptr<xylib::DataSet> newdata(xylib::load_file(path.string(), "canberra_cnf"));

  if (!newdata)
    throw std::runtime_error("ImporterCNF could not load file with xylib");

  std::string descr;
  date::sys_time<std::chrono::seconds> start_time;
  std::vector<double> calibration;
  int64_t rt_ms {0}, lt_ms {0};

//  DBG << "xylib.blocks =  " << newdata->get_block_count();
  for (int i = 0; i < newdata->get_block_count(); i++)
  {
    calibration.clear();

    for (uint32_t j = 0; j < newdata->get_block(i)->meta.size(); j++)
    {
      std::string key = newdata->get_block(i)->meta.get_key(j);
      std::string value = newdata->get_block(i)->meta.get(key);
//      DBG << "xylib.meta " << key << " = " << value;
      if (key.substr(0, 12) == "energy calib")
      {
//        DBG("calib = {}", value);
        calibration.push_back(std::stod(value));
      }
      else if (key == "description")
      {
//        DBG("descr = {}", value);
        descr = value;
      }
      else if (key == "date and time")
      {
        std::istringstream stream{value};
        stream >> date::parse("%a, %Y-%m-%d %H:%M:%S", start_time);
        if (stream.fail())
          throw std::runtime_error("failed to parse " + value);
//        DBG("parsing start time = '{}' -> {}", value, date::format("%FT%TZ", t));
      }
      else if (key == "real time (s)")
      {
        rt_ms = static_cast<int64_t>(std::stod(trim_copy(value)) * 1000.0);
//        DBG("parsing real time = '{}' -> {}", value, rt_ms);
      }
      else if (key == "live time (s)")
      {
        lt_ms = static_cast<int64_t>(std::stod(trim_copy(value)) * 1000.0);
//        DBG("parsing live time = '{}' -> {}", value, lt_ms);
      }
    }

    int column = newdata->get_block(i)->get_column_count();

    if (column == 2)
    {
//      DBG << "xylib.points = " << newdata->get_block(i)->get_point_count();
      for (int k = 0; k < newdata->get_block(i)->get_point_count(); k++)
      {
        double data = newdata->get_block(i)->get_column(column).get_value(k);
        DAQuiri::Entry new_entry;
        new_entry.first.resize(1);
        new_entry.first[0] = k;
        new_entry.second = PreciseFloat(data);

        entry_list.push_back(new_entry);
      }
    }
  }

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterCNF could not get a valid Histogram 1D from factory");

  hist->set_attribute(DAQuiri::Setting::text("description", descr));
  hist->set_attribute(DAQuiri::Setting("start_time", start_time));
  hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));
  hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));

  DAQuiri::CalibID from("energy","unknown","");
  DAQuiri::CalibID to("energy","unknown","keV");
  DAQuiri::Calibration new_calib(from, to);
  new_calib.function(std::make_shared<DAQuiri::Polynomial>(calibration));
//  DBG("calib = {}", new_calib.debug());
  DAQuiri::Detector det;
  det.set_calibration(new_calib);
  hist->set_detectors({det});

  hist->set_attribute(DAQuiri::Setting::text("value_latch", "energy"));
  hist->set_attribute(DAQuiri::Setting::text("name", path.stem().string()));
  hist->set_attribute(DAQuiri::Setting::boolean("visible", true));

  hist->import(*this);

  project->add_consumer(hist);
}
