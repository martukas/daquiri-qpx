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

#include "pixie4.h"
#include <core/producer_factory.h>
#include <boost/filesystem.hpp>
#include <core/util/custom_logger.h>
#include <core/util/timer.h>
#include <bitset>

#define SLOT_WAVE_OFFSET      7
#define NUMBER_OF_CHANNELS    4

void Pixie4::RunSetup::set_num_modules(uint16_t nmod)
{
  indices.resize(nmod);
  for (auto& i : indices)
    if (i.size() != NUMBER_OF_CHANNELS)
      i.resize(NUMBER_OF_CHANNELS, -1);
}

Pixie4::Pixie4()
{
  status_ = ProducerStatus::loaded | ProducerStatus::can_boot;
  boot_files_.resize(7);
  running_.store(false);
  terminating_.store(false);
}

Pixie4::~Pixie4()
{
  daq_stop();
  if (runner_.joinable())
    runner_.join();
  if (parser_.joinable())
    parser_.join();
  if (raw_queue_ != nullptr)
  {
    raw_queue_->stop();
    delete raw_queue_;
  }
}

StreamManifest Pixie4::stream_manifest() const
{
  StreamManifest ret;
  // \todo make this work
//  for (auto& s : streams_)
//  {
//    if (!s.parser)
//      continue;
//    for (auto m : s.parser->stream_manifest())
//      ret[m.first] = m.second;
//  }
  return ret;
}

bool Pixie4::daq_start(SpillQueue out_queue)
{
  if (running_.load() || parser_.joinable() || runner_.joinable())
    return false;

  terminating_.store(false);

  for (size_t i = 0; i < run_setup.indices.size(); i++)
  {
    PixieAPI.clear_EM(i);
    // double-buffered:
    PixieAPI.set_mod(i, "DBLBUFCSR", 1);
    PixieAPI.set_mod(i, "MODULE_CSRA", 0);
  }

  for (size_t i = 0; i < run_setup.indices.size(); ++i)
  {
    PixieAPI.set_mod(i, "RUN_TYPE", run_setup.type);
    PixieAPI.set_mod(i, "MAX_EVENTS", 0);
  }

  PixieAPI.reset_counters_next_run(); //assume new run

  raw_queue_ = new SpillMultiqueue(false, 100); // \todo params?

  parser_ = std::thread(&Pixie4::worker_parse, this, raw_queue_, out_queue);

  runner_ = std::thread(&Pixie4::worker_run_dbl, this, raw_queue_);

  return true;
}

bool Pixie4::daq_stop()
{
  if (!running_.load())
    return false;

  terminating_.store(true);

  if (runner_.joinable())
    runner_.join();

  Timer::wait_ms(500);
  while (raw_queue_->size() > 0)
    Timer::wait_ms(1000);
  Timer::wait_ms(500);
  raw_queue_->stop();
  Timer::wait_ms(500);

  if (parser_.joinable())
    parser_.join();

  delete raw_queue_;
  raw_queue_ = nullptr;

  terminating_.store(false);
  return true;
}

bool Pixie4::daq_running()
{
  return (running_.load());
}

EventModel Pixie4::model_hit(uint16_t runtype)
{
  EventModel h;
  h.timebase = TimeBase(1000, 75);
  h.add_value("energy", 16);
  h.add_value("front", 1);

  if (runtype < 259)
  {
    h.add_value("XIA_PSA", 16);
    h.add_value("user_PSA", 16);
  }

  if (runtype == 256)
    h.add_trace("waveform", {{1024}});

  return h;
}

//void Pixie4::fill_stats(std::map<int16_t, StatsUpdate> &all_stats, uint8_t module)
//{
//  StatsUpdate stats;
//  stats.items["native_time"] = PixieAPI.get_mod(module, "TOTAL_TIME");
//  stats.model_hit = model_hit(run_setup.type);
//  //tracelength?!?!?!
//  for (uint16_t i=0; i < run_setup.indices[module].size(); ++i)
//  {
//    stats.source_channel         = run_setup.indices[module][i];
//    stats.items["trigger_count"] = PixieAPI.get_chan(module, i, "FAST_PEAKS");
//    double live_time  = PixieAPI.get_chan(module, i, "LIVE_TIME");
//    stats.items["live_time"]    = live_time -
//        PixieAPI.get_chan(module, i, "SFDT");
//    stats.items["live_trigger"] = live_time -
//        PixieAPI.get_chan(module, i, "FTDT");
//    all_stats[stats.source_channel] = stats;
//  }
//}


Setting Pixie4::settings() const
{
  Setting set;

  if (set.id() != plugin_name())
    return set;

  for (auto& q : set.branches)
  {
    if (set.metadata().type() == SettingType::command)
    {
      if (status_ & ProducerStatus::booted)
        set.metadata().remove_flag("readonly");
      else
        set.metadata().set_flag("readonly");
    }

    if (q.id() == "Pixie4/Run settings")
      read_run_settings(q);
    else if (q.id() == "Pixie4/Files")
      read_files(q);
    else if (q.id() == "Pixie4/System")
      read_system(q);
  }
}

void Pixie4::read_run_settings(Setting& set) const
{
  for (auto& k : set.branches)
  {
    if (k.id() == "Pixie4/Run type")
      k.set_int(run_setup.type);
    if (k.id() == "Pixie4/Poll interval")
      k.set_int(run_setup.run_poll_interval_ms);
  }
}

void Pixie4::read_files(Setting& set) const
{
  for (auto& k : set.branches)
  {
    if (status_ & ProducerStatus::booted)
      k.metadata().set_flag("readonly");
    else
      k.metadata().remove_flag("readonly");
    if (k.id() == "Pixie4/Files/XIA_path")
      k.set_text(XIA_file_directory_);
    else if ((k.metadata().type() == SettingType::text) &&
        (k.metadata().get_num("address", 0) > 0) && (k.metadata().get_num("address", 9) < 8))
      k.set_text(boot_files_[k.metadata().get_num("address", 0) - 1]);
  }
}

void Pixie4::read_system(Setting& set) const
{
  for (auto& k : set.branches)
  {
    if (!(status_ & ProducerStatus::booted) &&
        setting_definitions_.count(k.id()) &&
        !setting_definitions_.at(k.id()).has_flag("readonly"))
      k.metadata().remove_flag("readonly");
    else
      k.metadata().set_flag("readonly");

    if (k.metadata().type() == SettingType::stem)
      read_module(k);
    else
      set_value(k, PixieAPI.get_sys(k.metadata().get_num("address", -1)));
  }
}

void Pixie4::read_module(Setting& set) const
{
  int16_t modnum = set.metadata().get_num("address", -1);
  if (!PixieAPI.module_valid(modnum))
  {
    WARN("<Pixie4> module address out of bounds, ignoring branch {}", modnum);
    return;
  }
  int filterrange = PixieAPI.get_mod(modnum, "FILTER_RANGE");
  for (auto& p : set.branches)
  {
    if (p.metadata().type() == SettingType::stem)
      read_channel(p, modnum, filterrange);
    else
      set_value(p, PixieAPI.get_mod(modnum, p.metadata().get_num("address", -1)));
  }
}

void Pixie4::read_channel(Setting& set, uint16_t modnum, int filterrange) const
{
  int16_t channum = set.metadata().get_num("address", -1);
  if (!PixieAPI.channel_valid(channum))
  {
    WARN("<Pixie4> channel address out of bounds, ignoring branch {}", channum);
    return;
  }
  for (auto& o : set.branches)
  {
    set_value(o, PixieAPI.get_chan(modnum, channum, o.metadata().get_num("address", -1)));
    if (o.metadata().get_string("name", "") == "ENERGY_RISETIME")
    {
      o.metadata().set_val("step", static_cast<double>(pow(2, filterrange)) / 75.0);
      o.metadata().set_val("minimum", 2 * o.metadata().step<double>());
      o.metadata().set_val("maximum", 124 * o.metadata().step<double>());
    }
    else if (o.metadata().get_string("name", "") == "ENERGY_FLATTOP")
    {
      o.metadata().set_val("step", static_cast<double>(pow(2, filterrange)) / 75.0);
      o.metadata().set_val("minimum", 3 * o.metadata().step<double>());
      o.metadata().set_val("maximum", 125 * o.metadata().step<double>());
    }
  }
}

void Pixie4::rebuild_structure(Setting& set)
{
  auto maxmod = set.find(Setting("Pixie4/System/MAX_NUMBER_MODULES"), Match::id);
  maxmod.set_int(std::min(std::max(int(maxmod.get_int()), 1),
                          int(Pixie4Wrapper::hard_max_modules())));
  set.set(maxmod, Match::id);

  auto oslots = set.find_all(Setting("Pixie4/System/SLOT_WAVE"), Match::id);
  std::vector<Setting> all_slots(oslots.begin(), oslots.end());

  all_slots.resize(maxmod.get_int(), get_rich_setting("Pixie4/System/SLOT_WAVE"));

  for (size_t i = 0; i < all_slots.size(); ++i)
    all_slots[i].metadata().set_val("address", SLOT_WAVE_OFFSET + i);

  set.erase(Setting("Pixie4/System/SLOT_WAVE"), Match::id);
  for (auto& q : all_slots)
    set.branches.add_a(q);

  auto totmod = set.find(Setting("Pixie4/System/NUMBER_MODULES"), Match::id);
  totmod.metadata().set_val("step", 1);
  totmod.metadata().set_val("maximum", Pixie4Wrapper::hard_max_modules());
  totmod.set_number(0);
  for (auto s : all_slots)
    if (s.get_int() > 0)
      totmod++;
  set.set(totmod, Match::id);

  auto omods = set.find_all(Setting("Pixie4/System/module"), Match::id);
  std::vector<Setting> old_modules(omods.begin(), omods.end());

  old_modules.resize(totmod.get_int(), default_module());

  set.erase(Setting("Pixie4/System/module"), Match::id);
  for (auto& q : old_modules)
    set.branches.add_a(q);

  PixieAPI.set_num_modules(totmod.get_int());
  run_setup.set_num_modules(totmod.get_int());
}

Setting Pixie4::default_module() const
{
  auto chan = get_rich_setting("Pixie4/System/module/channel");
  auto mod = get_rich_setting("Pixie4/System/module");
  for (int j = 0; j < NUMBER_OF_CHANNELS; ++j)
  {
    chan.metadata().set_val("address", j);
    mod.branches.add_a(chan);
  }
  return mod;
}

void Pixie4::reindex_modules(Setting& set)
{
  int module_address = 0;
  for (auto& m : set.branches)
  {
    if (m.id() != "Pixie4/System/module")
      continue;

    std::set<int32_t> new_set;
    int channel_address = 0;
    for (auto& c : m.branches)
    {
      if (c.id() != "Pixie4/System/module/channel")
        continue;

      if (c.indices().size() > 1)
      {
        int32_t i = *c.indices().begin();
        c.set_indices({i});
      }
      if (c.indices().size() > 0)
        new_set.insert(*c.indices().begin());
      c.metadata().set_val("address", channel_address++);
    }

    m.set_indices(new_set);
    m.metadata().set_val("address", module_address++);
  }
}

bool Pixie4::execute_command(const std::string& id)
{
  if (id == "Pixie4/Measure baselines")
    return PixieAPI.control_measure_baselines();
  else if (id == "Pixie4/Adjust offsets")
    return PixieAPI.control_adjust_offsets();
  else if (id == "Pixie4/Compute Tau")
    return PixieAPI.control_find_tau();
  else if (id == "Pixie4/Compute BLCUT")
    return PixieAPI.control_compute_BLcut();
  return false;
}

bool Pixie4::change_XIA_path(const std::string& xia_path)
{
  if (XIA_file_directory_ == xia_path)
    return false;

  boost::filesystem::path path(xia_path);
  for (auto& f : boot_files_)
  {
    boost::filesystem::path full_path = path / boost::filesystem::path(f).filename();
    f = full_path.string();
  }
  XIA_file_directory_ = xia_path;
  return true;
}

void Pixie4::settings(const Setting& setting)
{
  Setting set = setting;

  if (set.id() != plugin_name())
    return;

  set.enrich(setting_definitions_);

  for (auto& q : set.branches)
  {
    if ((q.metadata().type() == SettingType::command) && (q.get_int() == 1))
    {
      q.set_int(0);
      execute_command(q.id());
    }
    else if (q.id() == "Pixie4/Run settings")
      write_run_settings(q);
    else if ((q.id() == "Pixie4/Files") && !(status_ & ProducerStatus::booted))
      write_files(q);
    else if (q.id() == "Pixie4/System")
      write_system(q);
  }
}

void Pixie4::write_run_settings(Setting& set)
{
  for (auto& k : set.branches)
  {
    if (k.id() == "Pixie4/Run settings/Run type")
      run_setup.type = k.get_int();
    else if (k.id() == "Pixie4/Run settings/Poll interval")
      run_setup.run_poll_interval_ms = k.get_int();
  }
}

void Pixie4::write_files(Setting& set)
{
  for (auto& k : set.branches)
  {
    if ((k.id() == "Pixie4/Files/XIA_path") && change_XIA_path(k.get_text()))
      break;
    else if ((k.metadata().type() == SettingType::text) &&
        (k.metadata().get_num("address", 0) > 0) && (k.metadata().get_num("address", 9) < 8))
      boot_files_[k.metadata().get_num("address", 0) - 1] = k.get_text();
  }
}

void Pixie4::write_system(Setting& set)
{
  if (!(status_ & ProducerStatus::booted))
    rebuild_structure(set);

  reindex_modules(set);

  for (auto& k : set.branches)
  {
    if (k.metadata().type() == SettingType::stem)
      write_module(k);
    else if (!k.metadata().has_flag("readonly") &&
        (PixieAPI.get_sys(k.metadata().get_num("address", -1)) != get_value(k)))
      PixieAPI.set_sys(k.metadata().get_string("name", ""), get_value(k));
  }
}

void Pixie4::write_module(Setting& set)
{
  int16_t modnum = set.metadata().get_num("address", -1);
  if (!PixieAPI.module_valid(modnum))
    return;

  for (auto& p : set.branches)
  {
    if (p.metadata().type() != SettingType::stem)
      p.set_indices(set.indices());

    if (p.metadata().type() == SettingType::stem)
      write_channel(p, modnum);
    else if (!p.metadata().has_flag("readonly") &&
        (PixieAPI.get_mod(modnum, p.metadata().get_num("address", -1)) != get_value(p)))
      PixieAPI.set_mod(modnum, p.metadata().get_string("name", ""), get_value(p));
  }
}

void Pixie4::write_channel(Setting& set, uint16_t modnum)
{
  int16_t channum = set.metadata().get_num("address", -1);
  if (!Pixie4Wrapper::channel_valid(channum))
    return;

  int det = -1;
  if (!set.indices().empty())
    det = *set.indices().begin();
  run_setup.indices[modnum][channum] = det;

  for (auto& o : set.branches)
  {
    o.set_indices({det});

    if (!o.metadata().has_flag("readonly") &&
        (PixieAPI.get_chan(modnum, channum, o.metadata().get_num("address", -1)) != get_value(o)))
      PixieAPI.set_chan(modnum, channum, o.metadata().get_string("name", ""), get_value(o));
  }
}

void Pixie4::boot()
{
  if (!(status_ & ProducerStatus::can_boot))
  {
    WARN("<Pixie4> Cannot boot Pixie-4. Failed flag check (can_boot == 0)");
    return;
  }

  status_ = ProducerStatus::loaded | ProducerStatus::can_boot;

  bool valid_files = true;
  for (int i = 0; i < 7; i++)
    if (!boost::filesystem::exists(boot_files_[i]))
    {
      ERR("<Pixie4> Boot file {} not found", boot_files_[i]);
      valid_files = false;
    }

  if (!valid_files)
  {
    ERR("<Pixie4> Problem with boot files. Boot aborting.");
    return;
  }

  if (!PixieAPI.boot(boot_files_))
    return;

  status_ = ProducerStatus::loaded | ProducerStatus::booted
      | ProducerStatus::can_run | ProducerStatus::can_oscil;
}

void Pixie4::die()
{
  daq_stop();
  status_ = ProducerStatus::loaded | ProducerStatus::can_boot;
}

OscilData Pixie4::oscilloscope()
{
  OscilData result;

  uint32_t* oscil_data;

  for (size_t m = 0; m < run_setup.indices.size(); ++m)
  {
    if ((oscil_data = PixieAPI.control_collect_ADC(m)) != nullptr)
    {
      for (size_t i = 0; i < run_setup.indices[m].size(); i++)
      {
        std::vector<uint16_t> trace = std::vector<uint16_t>(oscil_data + (i * Pixie4Wrapper::max_buf_len),
                                                            oscil_data + ((i + 1) * Pixie4Wrapper::max_buf_len));
        if ((i < run_setup.indices[m].size()) &&
            (run_setup.indices[m][i] >= 0))
        {
          EventModel hm;
          hm.timebase = TimeBase(PixieAPI.get_chan(m, i, "XDT") * 1000, 1); //us to ns
          hm.add_trace("waveform", {{Pixie4Wrapper::max_buf_len}});
          // \todo use index run_setup.indices[m][i]
          Event tr(hm);
          // tr.trace(0) = trace;
          result["hack_id"] = tr;
        }
      }
    }

  }

  delete[] oscil_data;
  return result;
}

void Pixie4::get_all_settings()
{
  if (status_ & ProducerStatus::booted)
    PixieAPI.get_all_settings();
}

double Pixie4::get_value(const Setting& s)
{
  if (s.metadata().type() == SettingType::floating)
    return s.get_number();
  else if ((s.metadata().type() == SettingType::integer)
      || (s.metadata().type() == SettingType::boolean)
      || (s.metadata().type() == SettingType::menu)
      || (s.metadata().type() == SettingType::binary))
    return s.get_int();
}

void Pixie4::set_value(Setting& s, double val)
{
  if (s.metadata().type() == SettingType::floating)
    s.set_number(val);
  else if ((s.metadata().type() == SettingType::integer)
      || (s.metadata().type() == SettingType::boolean)
      || (s.metadata().type() == SettingType::menu)
      || (s.metadata().type() == SettingType::binary))
    s.set_int(val);
}

//////////////////////////////
//////////////////////////////
/////  T H R E A D S   ///////
//////////////////////////////
//////////////////////////////

void Pixie4::worker_run_dbl(SpillQueue spill_queue)
{
  std::string stream_id_ = "dummy_id";

  auto pixie = PixieAPI;
  auto setup = run_setup;
  SpillPtr fetched_spill = std::make_shared<Spill>(stream_id_, Spill::Type::start);

  //Start run;
  running_.store(true);
  if (!pixie.start_run(setup.type))
  {
    running_.store(false);
    return;
  }

  //Push spill with initial channel stats
  pixie.get_all_stats();
  // \todo bring this back
//  for (size_t i=0; i < setup.indices.size(); i++)
//    callback->fill_stats(fetched_spill.stats, i);
//  for (auto &q : fetched_spill.stats)
//  {
//    q.second.lab_time = fetched_spill.time;
//    q.second.stats_type = Spill::Type::start;
//  }
  spill_queue->enqueue(fetched_spill);

  //Main data acquisition loop
  bool timeout = false;
  std::set<uint16_t> triggered_modules;
  while (!timeout || !triggered_modules.empty())
  {
    //keep polling, pause if no trigger
    triggered_modules.clear();
    while (!timeout && triggered_modules.empty())
    {
      triggered_modules = pixie.get_triggered_modules();
      if (triggered_modules.empty())
        Timer::wait_ms(setup.run_poll_interval_ms);
      timeout = terminating_.load();
    };

    //get time ASAP after trigger
    fetched_spill->time = std::chrono::system_clock::now();

    //if timeout, stop run, wait, poll again
    if (timeout)
    {
      pixie.stop_run(setup.type);
      Timer::wait_ms(setup.run_poll_interval_ms);
      triggered_modules = pixie.get_triggered_modules();
    }

    //get stats ASAP
    pixie.get_all_settings(triggered_modules);

    //Read events, push spills
    bool success = false;
    for (auto& q : triggered_modules)
    {
      fetched_spill = std::make_shared<Spill>(stream_id_, Spill::Type::running);
      fetched_spill->raw.resize(Pixie4Wrapper::list_mem_len_bytes, 0); // buffer resolution multiply
      if (pixie.read_EM_double_buffered(
          reinterpret_cast<uint32_t*>(fetched_spill->raw.data()), q))
        success = true;

      //callback->fill_stats(fetched_spill.stats, q);
//      for (auto &p : fetched_spill.stats)
//      {
//        p.second.lab_time = fetched_spill.time;
//        if (timeout)
//          p.second.stats_type = StatsUpdate::Type::stop;
//      }
      spill_queue->enqueue(fetched_spill);
    }

    if (!success)
      break;
  }

  //Push spill with end of acquisition stats, indicate stop
  fetched_spill = std::make_shared<Spill>(stream_id_, Spill::Type::stop);
  pixie.get_all_stats();
//  for (size_t i=0; i < setup.indices.size(); i++)
//    callback->fill_stats(fetched_spill.stats, i);
//  for (auto &q : fetched_spill.stats)
//  {
//    q.second.lab_time = fetched_spill.time;
//    q.second.stats_type = StatsUpdate::Type::stop;
//  }
  spill_queue->enqueue(fetched_spill);
  running_.store(false);
}

void Pixie4::worker_parse(SpillQueue in_queue, SpillQueue out_queue)
{
  std::string stream_id_ = "dummy_id";

  SpillPtr spill = std::make_shared<Spill>(stream_id_, Spill::Type::running);

  uint64_t all_hits = 0, cycles = 0;
  double total_us{0};

  while ((spill = in_queue->dequeue()) != NULL)
  {
    Timer parse_timer(true);

    if (spill->raw.size() > 0)
    {
      cycles++;
      uint16_t* buff16 = reinterpret_cast<uint16_t*>(spill->raw.data());
      uint32_t idx = 0, spill_hits = 0;

      while (true)
      {
        uint16_t buf_ndata = buff16[idx++];
        uint32_t buf_end = idx + buf_ndata - 1;

        if ((buf_ndata == 0)
            || (buf_ndata > Pixie4Wrapper::max_buf_len)
            || (buf_end > Pixie4Wrapper::list_mem_len16))
          break;

        uint16_t buf_module = buff16[idx++];
        uint16_t buf_format = buff16[idx++];
        uint16_t buf_timehi = buff16[idx++];
        uint16_t buf_timemi = buff16[idx++];
        idx++; //uint16_t buf_timelo = buff16[idx++]; unused
        uint16_t task_a = (buf_format & 0x0F00);
        uint16_t task_b = (buf_format & 0x000F);

        while ((task_a == 0x0100) && (idx < buf_end))
        {
          std::bitset<16> pattern(buff16[idx++]);

          uint16_t evt_time_hi = buff16[idx++];
          uint16_t evt_time_lo = buff16[idx++];

          std::multimap<uint64_t, Event> ordered;

          for (size_t i = 0; i < NUMBER_OF_CHANNELS; i++)
          {
            if (!pattern[i])
              continue;

            int16_t sourcechan = -1;
            if ((buf_module < run_setup.indices.size()) &&
                (i < run_setup.indices[buf_module].size()) &&
                (run_setup.indices[buf_module][i] >= 0))
              sourcechan = run_setup.indices[buf_module][i];

            // \\todo use sourcechan
            Event one_hit(model_hit(run_setup.type));

            uint64_t hi = buf_timehi;
            uint64_t mi = evt_time_hi;
            uint64_t lo = evt_time_lo;
            uint16_t chan_trig_time = lo;
            uint16_t chan_time_hi = hi;

            one_hit.set_value(1, pattern[4]); //Front panel input value

            if (task_b == 0x0000)
            {
              uint16_t trace_len = buff16[idx++] - 9;
              one_hit.set_value(0, buff16[idx++]); //energy
              one_hit.set_value(2, buff16[idx++]); //XIA_PSA
              one_hit.set_value(3, buff16[idx++]); //user_PSA
              idx += 3;
              hi = buff16[idx++]; //not always?
              //one_hit.trace(0) = std::vector<uint16_t>(buff16 + idx, buff16 + idx + trace_len);
              idx += trace_len;
            }
            else if (task_b == 0x0001)
            {
              idx++;
              chan_trig_time = buff16[idx++];
              one_hit.set_value(0, buff16[idx++]); //energy
              one_hit.set_value(2, buff16[idx++]); //XIA_PSA
              one_hit.set_value(3, buff16[idx++]); //user_PSA
              idx += 3;
              hi = buff16[idx++];
            }
            else if (task_b == 0x0002)
            {
              chan_trig_time = buff16[idx++];
              one_hit.set_value(0, buff16[idx++]); //energy
              one_hit.set_value(2, buff16[idx++]); //XIA_PSA
              one_hit.set_value(3, buff16[idx++]); //user_PSA
            }
            else if (task_b == 0x0003)
            {
              chan_trig_time = buff16[idx++];
              one_hit.set_value(0, buff16[idx++]); //energy
            }
            else
              ERR("<Pixie4::parser> Parsed event type invalid");

            if (!pattern[i + 8])
              one_hit.set_value(0, 0); //energy invalid or approximate

            //Corrections for overflow, page 30 in Pixie-4 user manual
            if (chan_trig_time > evt_time_lo)
              mi--;
            if (evt_time_hi < buf_timemi)
              hi++;
            if ((task_b == 0x0000) || (task_b == 0x0001))
              hi = chan_time_hi;
            lo = chan_trig_time;
            uint64_t time = (hi << 32) + (mi << 16) + lo;

            one_hit.set_time(time);

            if (sourcechan >= 0)
              ordered.emplace(one_hit.timestamp(), one_hit);
          }

          for (auto& q : ordered)
          {
            spill->events.last() = q.second;
            spill->events++;
            spill_hits++;
          }

        };
      }
      all_hits += spill_hits;
    }
    spill->raw.clear();
    spill->events.finalize();
    out_queue->enqueue(spill);
    total_us += parse_timer.us();
  }

  if (cycles == 0)
    DBG("<Pixie4::parser> Buffer queue closed without events");
  else
    DBG("<Pixie4::parser> Parsed {} hits, with avg time/spill: {} us",
        all_hits, total_us / cycles);
}
