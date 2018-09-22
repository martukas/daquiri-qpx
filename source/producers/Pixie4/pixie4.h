/*******************************************************************************
 *
 * This software was developed at the National Institute of Standards and
 * Technology (NIST) by employees of the Federal Government in the course
 * of their official duties. Pursuant to title 17 Section 105 of the
 * United States Code, this software is not subject to copyright protection
 * and is in the public domain. NIST assumes no responsibility whatsoever for
 * its use by other parties, and makes no guarantees, expressed or implied,
 * about its quality, reliability, or any other characteristic.
 *
 * This software can be redistributed and/or modified freely provided that
 * any derivative works bear some notice that they are derived from it, and
 * any modified versions bear some notice that they have been modified.
 *
 * Author(s):
 *      Martin Shetty (NIST)
 *
 ******************************************************************************/

#pragma once

#include <atomic>
#include <thread>
#include <core/producer.h>
#include "pixie4_api_wrapper.h"

using namespace DAQuiri;

class Pixie4 : public Producer
{
  
public:
  Pixie4();
  ~Pixie4();

  std::string plugin_name() const override {return "Pixie4";}

  void settings(const Setting &set) override;
  Setting settings() const override;

  void get_all_settings() override;
  void boot() override;
  void die() override;

  OscilData oscilloscope() override;
  StreamManifest stream_manifest() const override;

  bool daq_start(SpillQueue out_queue) override;
  bool daq_stop() override;
  bool daq_running() override;

protected:
  Pixie4Wrapper PixieAPI;
  std::string XIA_file_directory_;
  std::vector<std::string> boot_files_;

  struct RunSetup
  {
    std::vector<std::vector<int32_t>> indices;
    uint16_t type {0x103};
    int  run_poll_interval_ms {100};

    void set_num_modules(uint16_t nmod);
  }
  run_setup;

  // Threads
  std::atomic<bool> running_{false};
  std::atomic<bool> terminating_{false};
  std::thread runner_;
  std::thread parser_;
  SpillQueue raw_queue_ {nullptr};
  void worker_parse(SpillQueue in_queue, SpillQueue out_queue);
  void worker_run_dbl(SpillQueue spill_queue);

  // Helpers for daq
//  void fill_stats(std::map<int16_t, StatsUpdate>&, uint8_t module);
  static EventModel model_hit(uint16_t runtype);

  // Helpers for settings
  bool execute_command(const std::string& id);
  bool change_XIA_path(const std::string& xia_path);

  void read_run_settings(Setting &set) const;
  void read_files(Setting &set) const;
  void read_system(Setting &set) const;
  void read_module(Setting &set) const;
  void read_channel(Setting &set, uint16_t modnum, int filterrange) const;

  void write_run_settings(Setting &set);
  void write_files(Setting &set);
  void write_system(Setting &set);
  void write_module(Setting &set);
  void write_channel(Setting &set, uint16_t modnum);

  void rebuild_structure(Setting &set);
  void reindex_modules(Setting &set);
  Setting default_module() const;

  static double get_value(const Setting& s);
  static void set_value(Setting& s, double val);

  static SettingMeta px_setting(uint32_t address, std::string id_prefix, std::string name,
                                SettingType t, bool readonly = false);


private:
  //no copying
  void operator=(Pixie4 const&);
  Pixie4(const Pixie4&);

};
