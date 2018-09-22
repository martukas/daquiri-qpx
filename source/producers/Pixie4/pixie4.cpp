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

SettingMeta Pixie4::px_setting(uint32_t address, std::string id_prefix, std::string name,
                              SettingType t, bool readonly)
{
  SettingMeta ret(id_prefix + name, t, name);
  ret.set_val("address", address);
  if (readonly)
    ret.set_flag("readonly");
  return ret;
}

Pixie4::Pixie4()
{
  std::string r{plugin_name() + "/"};

  SettingMeta root(plugin_name(), SettingType::stem);
  root.set_flag("saveworthy");
  root.set_enums(0, {r + "Measure baselines",
                     r + "Adjust offsets",
                     r + "Compute Tau",
                     r + "Compute BLCUT",
                     r + "Run settings",
                     r + "Files",
                     r + "System"});
  add_definition(root);
  add_definition(px_setting(0, r, "Measure baselines", SettingType::command));
  add_definition(px_setting(0, r, "Adjust offsets", SettingType::command));
  add_definition(px_setting(0, r, "Compute BLCUT", SettingType::command));
  add_definition(px_setting(0, r, "Compute Tau", SettingType::command));


  SettingMeta run_settings(r + "Run settings", SettingType::stem, "Run settings");
  run_settings.set_flag("saveworthy");
  run_settings.set_enums(0, {
    r + "Run settings/Run type",
    r + "Run settings/Poll interval"
  });
  add_definition(run_settings);

  SettingMeta run_type(r + "Run settings/Run type", SettingType::menu, "Run type");
  run_settings.set_enum(256, "Traces");
  run_settings.set_enum(257, "Full");
  run_settings.set_enum(258, "PSA only");
  run_settings.set_enum(259, "Compressed");
  add_definition(run_type);

  SettingMeta poll_interval(r + "Run settings/Poll interval", SettingType::integer, "Poll interval");
  poll_interval.set_bounds(5, 50, 5000);
  poll_interval.set_val("units", "ms");
  add_definition(poll_interval);


  SettingMeta files(r + "Files", SettingType::stem, "Files");
  files.set_flag("saveworthy");
  files.set_enums(0, {
      r + "Files/XIA_path",
      r + "Files/Fippi",
      r + "Files/SysPixie",
      r + "Files/PixieBin",
      r + "Files/PxiBin",
      r + "Files/SetFile",
      r + "Files/PxiVar",
      r + "Files/PxiLst",
  });
  add_definition(files);

  SettingMeta xia_path(r + "Files/XIA_path", SettingType::text, "XIA_path");
  xia_path.set_flag("directory");
  add_definition(xia_path);

  SettingMeta Fippi(r + "Files/Fippi", SettingType::text, "Fippi");
  Fippi.set_flag("file");
  Fippi.set_val("wildcards", "Fippi binary (FippiP500.bin)");
  add_definition(Fippi);

  SettingMeta PixieBin(r + "Files/PixieBin", SettingType::text, "PixieBin");
  PixieBin.set_flag("file");
  PixieBin.set_val("wildcards", "Pixie binary (pixie.bin)");
  add_definition(PixieBin);

  SettingMeta PxiBin(r + "Files/PxiBin", SettingType::text, "PxiBin");
  PxiBin.set_flag("file");
  PxiBin.set_val("wildcards", "PXI code binary (PXIcode.bin)");
  add_definition(PixieBin);

  SettingMeta PxiLst(r + "Files/PxiLst", SettingType::text, "PxiLst");
  PxiLst.set_flag("file");
  PxiLst.set_val("wildcards", "Pxi list file (PXIcode.lst)");
  add_definition(PxiLst);

  SettingMeta PxiVar(r + "Files/PxiVar", SettingType::text, "PxiVar");
  PxiVar.set_flag("file");
  PxiVar.set_val("wildcards", "Pxi variable file (PXIcode.var)");
  add_definition(PxiVar);

  SettingMeta SetFile(r + "Files/PxiVar", SettingType::text, "SetFile");
  SetFile.set_flag("file");
  SetFile.set_val("wildcards", "Pixie settings file (*.set)");
  add_definition(SetFile);

  SettingMeta SysPixie(r + "Files/SysPixie", SettingType::text, "SysPixie");
  SysPixie.set_flag("file");
  SysPixie.set_val("wildcards", "Syspixie binary (syspixie*.bin)");
  add_definition(SysPixie);

  std::string rs{r + "System/"};

  SettingMeta system(plugin_name() + "/System", SettingType::stem, "System");
  system.set_flag("saveworthy");
  system.set_enums(0, {rs + "NUMBER_MODULES",
                       rs + "OFFLINE_ANALYSIS",
                       rs + "AUTO_PROCESSLMDATA",
                       rs + "MAX_NUMBER_MODULES",
                       rs + "C_LIBRARY_RELEASE",
                       rs + "C_LIBRARY_BUILD",
                       rs + "KEEP_CW",
                       rs + "SLOT_WAVE",
                       rs + "module"});
  add_definition(system);

  std::string rm{rs + "module/"};
  std::string rc{rm + "channel/"};

  auto BASELINE_PERCENT = px_setting(18, rc, "BASELINE_PERCENT", SettingType::floating);
  BASELINE_PERCENT.set_bounds(0.0, 100.0);
  add_definition(BASELINE_PERCENT);

  auto BINFACTOR = px_setting(14, rc, "BINFACTOR", SettingType::integer);
  BINFACTOR.set_bounds(0, 65535);
  BINFACTOR.set_val("description", "downshift before binning (Pixie4 MCA)");
  add_definition(BINFACTOR);

  auto BLAVG = px_setting(25, rc, "BLAVG", SettingType::floating, true);
  add_definition(BLAVG);

  auto BLCUT = px_setting(16, rc, "BLCUT", SettingType::integer);
  BLCUT.set_bounds(0, 32767);
  add_definition(BLCUT);

  auto CFD_THRESHOLD = px_setting(19, rc, "CFD_THRESHOLD", SettingType::floating);
  CFD_THRESHOLD.set_bounds(0, 65535);
  add_definition(CFD_THRESHOLD);

  auto CHANNEL_CSRA = px_setting(0, rc, "CHANNEL_CSRA", SettingType::binary);
  CHANNEL_CSRA.set_val("bits", 16);
  CHANNEL_CSRA.set_enum(0, "Respond to group triggers only");
  CHANNEL_CSRA.set_enums(2, {"Good channel (contribute to list mode)",
                             "Read always (contribute to list mode even without hit)",
                             "Enable trigger (event trigger)",
                             "Trigger on positive slope, else on negative",
                             "GFLT (Veto) required",
                             "Histogram energies to on-board MCA"});
  CHANNEL_CSRA.set_enums(9, {"Allow negative number for pulse height",
                             "Compute constant fraction timing (PSA)"});
  CHANNEL_CSRA.set_enums(12, {"GATE required (Pixie4 only)",
                              "Local trigger to latch time stamp, else distributed group trigger",
                              "Estimate energy when not hit (with read always)"});
  add_definition(CHANNEL_CSRA);

  auto CHANNEL_CSRB = px_setting(1, rc, "CHANNEL_CSRB", SettingType::binary);
  CHANNEL_CSRB.set_val("bits", 16);
  CHANNEL_CSRB.set_enum(0, "Call user written DSP code");
  CHANNEL_CSRB.set_enum(1, "Overwrite channel header (except Ndata, TrigTime, Energy) with URETVAL");
  add_definition(CHANNEL_CSRB);

  auto CHANNEL_CSRC = px_setting(21, rc, "CHANNEL_CSRC", SettingType::binary);
  CHANNEL_CSRC.set_val("bits", 16);
  CHANNEL_CSRC.set_enums(0, {"GFLT acceptance polarity",
                             "GATE acceptance polarity",
                             "Use GFLT for GATE",
                             "Disable pileup inspection",
                             "Disable out-of-range rejection",
                             "Invert piuleup inspection",
                             "Pause pileup inspection (426 ns)",
                             "Gate edge polarity inverted",
                             "Gate statistics (live only if GATE or VETO)",
                             "Invert GDT for GATE or VETO presence",
                             "No gate pulse (pattern only reports status of gate)",
                             "4x traces - each waveform sample is avg of 4 ADC samples",
                             "Veto out-of-range (to backplane)"});
  add_definition(CHANNEL_CSRC);

  auto COINC_DELAY = px_setting(24, rc, "COINC_DELAY", SettingType::floating);
  COINC_DELAY.set_bounds(0.0, 0.013333, 65535.0);
  COINC_DELAY.set_val("units", "μs");
  add_definition(COINC_DELAY);

  auto CURRENT_ICR = px_setting(36, rc, "CURRENT_ICR", SettingType::floating, true);
  CURRENT_ICR.set_val("units", "count/s");
  add_definition(CURRENT_ICR);

  auto CURRENT_OORF = px_setting(37, rc, "CURRENT_OORF", SettingType::floating, true);
  CURRENT_OORF.set_val("units", "%");
  add_definition(CURRENT_OORF);

  auto EMIN = px_setting(13, rc, "EMIN", SettingType::integer);
  EMIN.set_bounds(0, 65535);
  EMIN.set_val("description", "subtracted from computed energy in list mode");
  add_definition(EMIN);

  auto ENERGY_FLATTOP = px_setting(3, rc, "ENERGY_FLATTOP", SettingType::floating);
  ENERGY_FLATTOP.set_bounds(0.0, 0.05, 500.0);
  ENERGY_FLATTOP.set_val("units", "μs");
  ENERGY_FLATTOP.set_flag("optimize");
  add_definition(ENERGY_FLATTOP);

  auto ENERGY_RISETIME = px_setting(2, rc, "ENERGY_RISETIME", SettingType::floating);
  ENERGY_RISETIME.set_bounds(0.0, 0.05, 500.0);
  ENERGY_RISETIME.set_val("units", "μs");
  ENERGY_RISETIME.set_flag("optimize");
  add_definition(ENERGY_RISETIME);

  auto FAST_PEAKS = px_setting(28, rc, "FAST_PEAKS", SettingType::floating, true);
  add_definition(FAST_PEAKS);

  auto FTDT = px_setting(33, rc, "FTDT", SettingType::floating, true);
  FTDT.set_val("units", "s");
  add_definition(FTDT);

  auto GATE_COUNTS = px_setting(32, rc, "GATE_COUNTS", SettingType::floating, true);
  add_definition(GATE_COUNTS);

  auto GATE_DELAY = px_setting(23, rc, "GATE_DELAY", SettingType::floating);
  GATE_DELAY.set_bounds(0.013333, 0.013333, 3.4);
  GATE_DELAY.set_val("units", "μs");
  add_definition(GATE_DELAY);

  auto GATE_RATE = px_setting(31, rc, "GATE_RATE", SettingType::floating, true);
  GATE_RATE.set_val("units", "counts/s");
  add_definition(GATE_RATE);

  auto GATE_WINDOW = px_setting(22, rc, "GATE_WINDOW", SettingType::floating);
  GATE_WINDOW.set_bounds(0.013333, 0.013333, 3.4);
  GATE_WINDOW.set_val("units", "μs");
  add_definition(GATE_WINDOW);

  auto GDT = px_setting(35, rc, "GDT", SettingType::floating, true);
  GDT.set_val("units", "s");
  add_definition(GDT);

  auto INPUT_COUNT_RATE = px_setting(27, rc, "INPUT_COUNT_RATE", SettingType::floating, true);
  INPUT_COUNT_RATE.set_val("units", "counts/s");
  add_definition(INPUT_COUNT_RATE);

  auto INTEGRATOR = px_setting(20, rc, "INTEGRATOR", SettingType::menu);
  INTEGRATOR.set_enums(0, {"trapezoidal", "gap sum", "step", "2x gap", "4x gap", "8x gap"});
  add_definition(INTEGRATOR);

  auto LIVE_TIME = px_setting(26, rc, "LIVE_TIME", SettingType::floating, true);
  LIVE_TIME.set_val("units", "s");
  add_definition(LIVE_TIME);

  auto NOUT = px_setting(30, rc, "NOUT", SettingType::floating, true);
  add_definition(NOUT);

  auto OUTPUT_COUNT_RATE = px_setting(29, rc, "OUTPUT_COUNT_RATE", SettingType::floating, true);
  OUTPUT_COUNT_RATE.set_val("units", "counts/s");
  add_definition(OUTPUT_COUNT_RATE);

  auto PSA_END = px_setting(12, rc, "PSA_END", SettingType::floating);
  PSA_END.set_bounds(0.0, 0.013333, 13.1936);
  PSA_END.set_val("units", "μs");
  add_definition(PSA_END);

  auto PSA_START = px_setting(11, rc, "PSA_START", SettingType::floating);
  PSA_START.set_bounds(0.0, 0.013333, 13.1936);
  PSA_START.set_val("units", "μs");
  add_definition(PSA_START);

  auto PSM_GAIN_AVG = px_setting(38, rc, "PSM_GAIN_AVG", SettingType::floating, true);
  add_definition(PSM_GAIN_AVG);

  auto PSM_GAIN_AVG_LEN = px_setting(39, rc, "PSM_GAIN_AVG_LEN", SettingType::floating, true);
  add_definition(PSM_GAIN_AVG_LEN);

  // is this correct?
  auto PSM_GAIN_CORR = px_setting(42, rc, "PSM_GAIN_CORR", SettingType::floating);
  PSM_GAIN_CORR.set_bounds(0.0, 0.1, 1.0);
  add_definition(PSM_GAIN_CORR);

  auto PSM_TEMP_AVG = px_setting(40, rc, "PSM_TEMP_AVG", SettingType::floating, true);
  PSM_TEMP_AVG.set_val("units", "°C");
  add_definition(PSM_TEMP_AVG);

  auto PSM_TEMP_AVG_LEN = px_setting(41, rc, "PSM_TEMP_AVG_LEN", SettingType::floating);
  PSM_TEMP_AVG_LEN.set_bounds(0.0, 0.1, 1.0);
  add_definition(PSM_TEMP_AVG_LEN);

  auto SFDT = px_setting(34, rc, "SFDT", SettingType::floating, true);
  SFDT.set_val("units", "s");
  add_definition(SFDT);

  auto TAU = px_setting(15, rc, "TAU", SettingType::floating);
  TAU.set_bounds(0.0, 1.0, 65535.0); // min=1.5e-05
  TAU.set_val("units", "μs");
  add_definition(TAU);

  auto TRACE_DELAY = px_setting(10, rc, "TRACE_DELAY", SettingType::floating);
  TRACE_DELAY.set_bounds(0.0, 0.013333, 13.1936);
  TRACE_DELAY.set_val("units", "μs");
  add_definition(TRACE_DELAY);

  auto TRACE_LENGTH = px_setting(9, rc, "TRACE_LENGTH", SettingType::floating);
  TRACE_LENGTH.set_bounds(0.0, 0.013333, 13.1936);
  TRACE_LENGTH.set_val("units", "μs");
  add_definition(TRACE_LENGTH);

  auto TRIGGER_FLATTOP = px_setting(5, rc, "TRIGGER_FLATTOP", SettingType::floating);
  TRIGGER_FLATTOP.set_bounds(0.0, 0.05, 500.0);
  TRIGGER_FLATTOP.set_val("units", "μs");
  add_definition(TRIGGER_FLATTOP);

  auto TRIGGER_RISETIME = px_setting(4, rc, "TRIGGER_RISETIME", SettingType::floating);
  TRIGGER_RISETIME.set_bounds(0.0, 0.05, 500.0);
  TRIGGER_RISETIME.set_val("units", "μs");
  add_definition(TRIGGER_RISETIME);

  auto TRIGGER_THRESHOLD = px_setting(6, rc, "TRIGGER_THRESHOLD", SettingType::floating);
  TRIGGER_THRESHOLD.set_bounds(0, 1, 4095);
  add_definition(TRIGGER_THRESHOLD);

  auto VGAIN = px_setting(7, rc, "VGAIN", SettingType::floating);
  VGAIN.set_bounds(0.0, 0.000076, 65535.0);
  VGAIN.set_val("units", "V/V");
  VGAIN.set_flag("gain");
  add_definition(VGAIN);

  auto VOFFSET = px_setting(8, rc, "VOFFSET", SettingType::floating);
  VOFFSET.set_bounds(-2.5, 0.000076, 2.5);
  VOFFSET.set_val("units", "V");
  add_definition(VOFFSET);

  auto XDT = px_setting(17, rc, "XDT", SettingType::floating);
  XDT.set_bounds(0.053333, 0.013333, 873.77333);
  XDT.set_val("units", "μs");
  add_definition(XDT);

  SettingMeta MODULE_CHANNEL(rm + "channel", SettingType::stem);
  MODULE_CHANNEL.set_flag("saveworthy");
  MODULE_CHANNEL.set_enums(0, {rc + "CHANNEL_CSRA",
                               rc + "CHANNEL_CSRB",
                               rc + "ENERGY_RISETIME",
                               rc + "ENERGY_FLATTOP",
                               rc + "TRIGGER_RISETIME",
                               rc + "TRIGGER_FLATTOP",
                               rc + "TRIGGER_THRESHOLD",
                               rc + "VGAIN",
                               rc + "VOFFSET",
                               rc + "TRACE_LENGTH",
                               rc + "TRACE_DELAY",
                               rc + "PSA_START",
                               rc + "PSA_END",
                               rc + "EMIN",
                               rc + "BINFACTOR",
                               rc + "TAU",
                               rc + "BLCUT",
                               rc + "XDT",
                               rc + "BASELINE_PERCENT",
                               rc + "CFD_THRESHOLD",
                               rc + "INTEGRATOR",
                               rc + "CHANNEL_CSRC",
                               rc + "GATE_WINDOW",
                               rc + "GATE_DELAY",
                               rc + "COINC_DELAY",
                               rc + "BLAVG",
                               rc + "LIVE_TIME",
                               rc + "INPUT_COUNT_RATE",
                               rc + "FAST_PEAKS",
                               rc + "OUTPUT_COUNT_RATE",
                               rc + "NOUT",
                               rc + "GATE_RATE",
                               rc + "GATE_COUNTS",
                               rc + "FTDT",
                               rc + "SFDT",
                               rc + "GDT",
                               rc + "CURRENT_ICR",
                               rc + "CURRENT_OORF",
                               rc + "PSM_GAIN_AVG",
                               rc + "PSM_GAIN_AVG_LEN",
                               rc + "PSM_TEMP_AVG",
                               rc + "PSM_TEMP_AVG_LEN",
                               rc + "PSM_GAIN_CORR"});
  add_definition(MODULE_CHANNEL);

  auto ACTUAL_COINCIDENCE_WAIT = px_setting(6, rm, "ACTUAL_COINCIDENCE_WAIT", SettingType::floating);
  ACTUAL_COINCIDENCE_WAIT.set_bounds(13.333333, 13.333333, 87374663.6448);
  ACTUAL_COINCIDENCE_WAIT.set_val("units", "ns");
  add_definition(ACTUAL_COINCIDENCE_WAIT);

  auto BOARD_VERSION = px_setting(24, rm, "BOARD_VERSION", SettingType::floating, true);
  add_definition(BOARD_VERSION);

  auto BUFFER_HEAD_LENGTH = px_setting(16, rm, "BUFFER_HEAD_LENGTH", SettingType::integer, true);
  add_definition(BUFFER_HEAD_LENGTH);

  auto CHANNEL_HEAD_LENGTH = px_setting(18, rm, "CHANNEL_HEAD_LENGTH", SettingType::integer, true);
  add_definition(CHANNEL_HEAD_LENGTH);

  auto COINCIDENCE_PATTERN = px_setting(5, rm, "COINCIDENCE_PATTERN", SettingType::binary);
  COINCIDENCE_PATTERN.set_val("bits", 16);
  COINCIDENCE_PATTERN.set_enums(0, {"oooo", "ooo+", "oo+o", "oo++",
                                    "o+oo", "o+o+", "o++o", "o+++",
                                    "+ooo", "+oo+", "+o+o", "+o++",
                                    "++oo", "++o+", "+++o", "++++"});
  add_definition(COINCIDENCE_PATTERN);

  auto DBLBUFCSR = px_setting(14, rm, "DBLBUFCSR", SettingType::binary);
  DBLBUFCSR.set_val("bits", 16);
  DBLBUFCSR.set_enum(0, "Enable double buffer mode");
  DBLBUFCSR.set_enum(1, "Host has read buffer");
  DBLBUFCSR.set_enum(3, "Host should read first block, else should read second block");
  add_definition(DBLBUFCSR);

  auto DSP_BUILD = px_setting(27, rm, "DSP_BUILD", SettingType::floating, true);
  add_definition(DSP_BUILD);

  auto DSP_RELEASE = px_setting(26, rm, "DSP_RELEASE", SettingType::floating, true);
  add_definition(DSP_RELEASE);

  auto EVENT_HEAD_LENGTH = px_setting(17, rm, "EVENT_HEAD_LENGTH", SettingType::integer, true);
  add_definition(EVENT_HEAD_LENGTH);

  auto EVENT_RATE = px_setting(22, rm, "EVENT_RATE", SettingType::floating, true);
  EVENT_RATE.set_val("units", "counts/s");
  add_definition(EVENT_RATE);

  auto FILTER_RANGE = px_setting(11, rm, "FILTER_RANGE", SettingType::menu);
  FILTER_RANGE.set_enums(1, {"2 bits", "4 bits", "8 bits", "16 bits", "32 bits", "64 bits"});
  add_definition(FILTER_RANGE);

  auto FIPPI_ID = px_setting(28, rm, "FIPPI_ID", SettingType::floating, true);
  add_definition(FIPPI_ID);

  auto IN_SYNCH = px_setting(9, rm, "IN_SYNCH", SettingType::boolean);
  IN_SYNCH.set_val("description", "modules synchronized; clear to have system reset clocks at start of daq");
  add_definition(IN_SYNCH);

  auto MAX_EVENTS = px_setting(4, rm, "MAX_EVENTS", SettingType::integer);
  MAX_EVENTS.set_bounds(0, 100);
  add_definition(MAX_EVENTS);

  auto MIN_COINCIDENCE_WAIT = px_setting(7, rm, "MIN_COINCIDENCE_WAIT", SettingType::integer, true);
  MIN_COINCIDENCE_WAIT.set_val("units", "ticks");
  add_definition(MIN_COINCIDENCE_WAIT);

  auto MODULEPATTERN = px_setting(12, rm, "MODULEPATTERN", SettingType::binary);
  MODULEPATTERN.set_val("bits", 16);
  MODULEPATTERN.set_enums(4, {"Gate event acceptance on FRONT panel input",
                              "Gate event acceptance on LOCAL coincidence test",
                              "Gate event acceptance on backplane STATUS line",
                              "Gate event acceptance on GLOBAL coincidence test"});
  add_definition(MODULEPATTERN);

  auto MODULE_CSRA = px_setting(1, rm, "MODULE_CSRA", SettingType::binary);
  MODULE_CSRA.set_val("bits", 16);
  MODULE_CSRA.set_enums(1, {"acquire 32 buffers and write to EM, else only one at a time",
                            "backpane trigger distribution: see manual...",
                            "Pixie native MCA: bin sums to addback spectrum",
                            "Pixie native MCA: individual spectra only singles",
                            "front panel DSP-OUT distributed as veto signal to backplane",
                            "chan3 hit status contributes to backplane STATUS",
                            "polarity of front panel pulse counter"});
  MODULE_CSRA.set_enums(9, {"write NNSHAREPATTERN to left neighbor (should be PDM) during ControlTask5",
                            "Pixie-500 only: time stamps as 2ns, else 13.3ns"});
  MODULE_CSRA.set_enums(12, {"drive low TOKEN backplane if local coincidence fails",
                             "send hit patten to slot2 using PXI STAR trigger for each event",
                             "front panel DSP-OUT as input to STATUS on backplane (wire-OR)",
                             "backpane trigger distribution: see manual..."});
  add_definition(MODULE_CSRA);

  auto MODULE_CSRB = px_setting(2, rm, "MODULE_CSRB", SettingType::binary);
  MODULE_CSRB.set_val("bits", 16);
  MODULE_CSRB.set_enums(0, {"Execute user code routines programmed by user dsp"});
  add_definition(MODULE_CSRB);

  auto MODULE_CSRC = px_setting(15, rm, "MODULE_CSRC", SettingType::binary);
  MODULE_CSRC.set_val("bits", 16);
  add_definition(MODULE_CSRC);

  auto MODULE_FORMAT = px_setting(3, rm, "MODULE_FORMAT", SettingType::binary, true);
  MODULE_FORMAT.set_val("bits", 16);
  MODULE_FORMAT.set_val("description", "not used");
  add_definition(MODULE_FORMAT);

  auto MODULE_NUMBER = px_setting(0, rm, "MODULE_NUMBER", SettingType::integer, true);
  add_definition(MODULE_NUMBER);

  auto NNSHAREPATTERN = px_setting(13, rm, "NNSHAREPATTERN", SettingType::integer);
  NNSHAREPATTERN.set_bounds(0, 65535);
  NNSHAREPATTERN.set_val("description", "User-defined control word for PXI-PDM");
  add_definition(NNSHAREPATTERN);

  auto NUMBER_EVENTS = px_setting(20, rm, "NUMBER_EVENTS", SettingType::integer, true);
  add_definition(NUMBER_EVENTS);

  auto OUTPUT_BUFFER_LENGTH = px_setting(19, rm, "OUTPUT_BUFFER_LENGTH", SettingType::integer, true);
  add_definition(OUTPUT_BUFFER_LENGTH);

  auto PDM_MASKA = px_setting(31, rm, "PDM_MASKA", SettingType::integer);
  PDM_MASKA.set_bounds(0, 65535);
  add_definition(PDM_MASKA);

  auto PDM_MASKB = px_setting(32, rm, "PDM_MASKB", SettingType::integer);
  PDM_MASKB.set_bounds(0, 65535);
  add_definition(PDM_MASKB);

  auto PDM_MASKC = px_setting(33, rm, "PDM_MASKC", SettingType::integer);
  PDM_MASKB.set_bounds(0, 65535);
  add_definition(PDM_MASKC);

  auto RUN_TIME = px_setting(21, rm, "RUN_TIME", SettingType::floating, true);
  RUN_TIME.set_val("units", "s");
  add_definition(RUN_TIME);

  auto RUN_TYPE = px_setting(10, rm, "RUN_TYPE", SettingType::menu);
  RUN_TYPE.set_enum(0, "Slow control run");
  RUN_TYPE.set_enum(256, "Traces");
  RUN_TYPE.set_enum(257, "Full");
  RUN_TYPE.set_enum(258, "PSA only");
  RUN_TYPE.set_enum(259, "Compressed");
  RUN_TYPE.set_enum(769, "Pixie MCA");
  add_definition(RUN_TYPE);

  auto SERIAL_NUMBER = px_setting(25, rm, "SERIAL_NUMBER", SettingType::floating, true);
  add_definition(SERIAL_NUMBER);

  auto SYNCH_WAIT = px_setting(8, rm, "SYNCH_WAIT", SettingType::boolean);
  SYNCH_WAIT.set_val("description", "wait for all modules ready before starting daq");
  add_definition(SYNCH_WAIT);

  auto SYSTEM_ID = px_setting(29, rm, "SYSTEM_ID", SettingType::floating, true);
  add_definition(SYSTEM_ID);

  auto TOTAL_TIME = px_setting(23, rm, "TOTAL_TIME", SettingType::floating, true);
  TOTAL_TIME.set_val("units", "s");
  add_definition(TOTAL_TIME);

  auto XET_DELAY = px_setting(30, rm, "XET_DELAY", SettingType::integer);
  XET_DELAY.set_bounds(0, 65535);
  XET_DELAY.set_val("description", "delay for generated event trigger from front panel to backplane");
  add_definition(XET_DELAY);

  SettingMeta MODULE(rs + "module", SettingType::stem);
  MODULE.set_flag("saveworthy");
  MODULE.set_enums(0, {rm + "MODULE_NUMBER",
                       rm + "MODULE_CSRA",
                       rm + "MODULE_CSRB",
                       rm + "MODULE_FORMAT",
                       rm + "MAX_EVENTS",
                       rm + "COINCIDENCE_PATTERN",
                       rm + "ACTUAL_COINCIDENCE_WAIT",
                       rm + "MIN_COINCIDENCE_WAIT",
                       rm + "SYNCH_WAIT",
                       rm + "IN_SYNCH",
                       rm + "RUN_TYPE",
                       rm + "FILTER_RANGE",
                       rm + "MODULEPATTERN",
                       rm + "NNSHAREPATTERN",
                       rm + "DBLBUFCSR",
                       rm + "MODULE_CSRC",
                       rm + "BUFFER_HEAD_LENGTH",
                       rm + "EVENT_HEAD_LENGTH",
                       rm + "CHANNEL_HEAD_LENGTH",
                       rm + "OUTPUT_BUFFER_LENGTH",
                       rm + "NUMBER_EVENTS",
                       rm + "RUN_TIME",
                       rm + "EVENT_RATE",
                       rm + "TOTAL_TIME",
                       rm + "BOARD_VERSION",
                       rm + "SERIAL_NUMBER",
                       rm + "DSP_RELEASE",
                       rm + "DSP_BUILD",
                       rm + "FIPPI_ID",
                       rm + "SYSTEM_ID",
                       rm + "XET_DELAY",
                       rm + "PDM_MASKA",
                       rm + "PDM_MASKB",
                       rm + "PDM_MASKC",
                       rm + "channel"});
  add_definition(MODULE);

  auto AUTO_PROCESSLMDATA = px_setting(2, rs, "AUTO_PROCESSLMDATA", SettingType::boolean, true);
  add_definition(AUTO_PROCESSLMDATA);

  auto C_LIBRARY_BUILD = px_setting(5, rs, "C_LIBRARY_BUILD", SettingType::floating, true);
  add_definition(C_LIBRARY_BUILD);

  auto C_LIBRARY_RELEASE = px_setting(4, rs, "C_LIBRARY_RELEASE", SettingType::floating, true);
  add_definition(C_LIBRARY_RELEASE);

  auto KEEP_CW = px_setting(6, rs, "KEEP_CW", SettingType::boolean);
  add_definition(KEEP_CW);

  auto MAX_NUMBER_MODULES = px_setting(3, rs, "MAX_NUMBER_MODULES", SettingType::integer);
  MAX_NUMBER_MODULES.set_bounds(1, 56);
  add_definition(MAX_NUMBER_MODULES);

  auto NUMBER_MODULES = px_setting(0, rs, "NUMBER_MODULES", SettingType::integer, true);
  add_definition(NUMBER_MODULES);

  auto OFFLINE_ANALYSIS = px_setting(1, rs, "OFFLINE_ANALYSIS", SettingType::boolean, true);
  add_definition(OFFLINE_ANALYSIS);

  auto SLOT_WAVE = px_setting(0, rs, "SLOT_WAVE", SettingType::integer);
  SLOT_WAVE.set_bounds(0, 7);
  add_definition(SLOT_WAVE);


//  root.set_flag("producer");

  boot_files_.resize(7);
  status_ = ProducerStatus::loaded | ProducerStatus::can_boot;
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
      o.metadata().set_val("min", 2 * o.metadata().step<double>());
      o.metadata().set_val("max", 124 * o.metadata().step<double>());
    }
    else if (o.metadata().get_string("name", "") == "ENERGY_FLATTOP")
    {
      o.metadata().set_val("step", static_cast<double>(pow(2, filterrange)) / 75.0);
      o.metadata().set_val("min", 3 * o.metadata().step<double>());
      o.metadata().set_val("max", 125 * o.metadata().step<double>());
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
  totmod.metadata().set_val("max", Pixie4Wrapper::hard_max_modules());
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
