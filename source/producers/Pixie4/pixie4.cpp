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
  std::string r{plugin_name() + "/"};
  std::string rm{r + "module/"};
  std::string rc{rm + "channel/"};

  SettingMeta BASELINE_PERCENT(rc + "BASELINE_PERCENT", SettingType::floating, "BASELINE_PERCENT");
  BASELINE_PERCENT.set_val("address", 18);
  BASELINE_PERCENT.set_val("min", 0);
  BASELINE_PERCENT.set_val("max", 100);
  BASELINE_PERCENT.set_val("step", 1);
  add_definition(BASELINE_PERCENT);

  SettingMeta BINFACTOR(rc + "BINFACTOR", SettingType::integer, "BINFACTOR");
  BINFACTOR.set_val("address", 14);
  BINFACTOR.set_val("min", 0);
  BINFACTOR.set_val("max", 65535);
  BINFACTOR.set_val("step", 1);
  BINFACTOR.set_val("description", "downshift before binning (Pixie4 MCA)");
  add_definition(BINFACTOR);

  SettingMeta BLAVG(rc + "BLAVG", SettingType::floating, "BLAVG");
  BLAVG.set_val("address", 25);
  BLAVG.set_flag("readonly");
  add_definition(BLAVG);

  SettingMeta BLCUT(rc + "BLCUT", SettingType::integer, "BLCUT");
  BLCUT.set_val("address", 16);
  BLCUT.set_val("min", 0);
  BLCUT.set_val("max", 32767);
  BLCUT.set_val("step", 1);
  add_definition(BLCUT);

  SettingMeta CFD_THRESHOLD(rc + "CFD_THRESHOLD", SettingType::floating, "CFD_THRESHOLD");
  CFD_THRESHOLD.set_val("address", 19);
  CFD_THRESHOLD.set_val("min", 0);
  CFD_THRESHOLD.set_val("max", 65535);
  CFD_THRESHOLD.set_val("step", 1);
  add_definition(CFD_THRESHOLD);

  SettingMeta CHANNEL_CSRA(rc + "CHANNEL_CSRA", SettingType::binary, "CHANNEL_CSRA");
  CHANNEL_CSRA.set_val("address", 0);
  CHANNEL_CSRA.set_val("bits", 16);
  CHANNEL_CSRA.set_enum(0, "Respond to group triggers only");
  CHANNEL_CSRA.set_enum(2, "Good channel (contribute to list mode)");
  CHANNEL_CSRA.set_enum(3, "Read always (contribute to list mode even without hit)");
  CHANNEL_CSRA.set_enum(4, "Enable trigger (event trigger)");
  CHANNEL_CSRA.set_enum(5, "Trigger on positive slope, else on negative");
  CHANNEL_CSRA.set_enum(6, "GFLT (Veto) required");
  CHANNEL_CSRA.set_enum(7, "Histogram energies to on-board MCA");
  CHANNEL_CSRA.set_enum(9, "Allow negative number for pulse height");
  CHANNEL_CSRA.set_enum(10, "Compute constant fraction timing (PSA)");
  CHANNEL_CSRA.set_enum(12, "GATE required (Pixie4 only)");
  CHANNEL_CSRA.set_enum(13, "Local trigger to latch time stamp, else distributed group trigger");
  CHANNEL_CSRA.set_enum(14, "Estimate energy when not hit (with read always)");
  add_definition(CHANNEL_CSRA);

  SettingMeta CHANNEL_CSRB(rc + "CHANNEL_CSRB", SettingType::binary, "CHANNEL_CSRB");
  CHANNEL_CSRB.set_val("address", 1);
  CHANNEL_CSRB.set_val("bits", 16);
  CHANNEL_CSRB.set_enum(0, "Call user written DSP code");
  CHANNEL_CSRB.set_enum(1, "Overwrite channel header (except Ndata, TrigTime, Energy) with URETVAL");
  add_definition(CHANNEL_CSRB);

  SettingMeta CHANNEL_CSRC(rc + "CHANNEL_CSRC", SettingType::binary, "CHANNEL_CSRC");
  CHANNEL_CSRC.set_val("address", 21);
  CHANNEL_CSRC.set_val("bits", 16);
  CHANNEL_CSRC.set_enum(0, "GFLT acceptance polarity");
  CHANNEL_CSRC.set_enum(1, "GATE acceptance polarity");
  CHANNEL_CSRC.set_enum(2, "Use GFLT for GATE");
  CHANNEL_CSRC.set_enum(3, "Disable pileup inspection");
  CHANNEL_CSRC.set_enum(4, "Disable out-of-range rejection");
  CHANNEL_CSRC.set_enum(5, "Invert piuleup inspection");
  CHANNEL_CSRC.set_enum(6, "Pause pileup inspection (426 ns)");
  CHANNEL_CSRC.set_enum(7, "Gate edge polarity inverted");
  CHANNEL_CSRC.set_enum(8, "Gate statistics (live only if GATE or VETO)");
  CHANNEL_CSRC.set_enum(9, "Invert GDT for GATE or VETO presence");
  CHANNEL_CSRC.set_enum(10, "No gate pulse (pattern only reports status of gate)");
  CHANNEL_CSRC.set_enum(11, "4x traces - each waveform sample is avg of 4 ADC samples");
  CHANNEL_CSRC.set_enum(12, "Veto out-of-range (to backplane)");
  add_definition(CHANNEL_CSRC);

  SettingMeta COINC_DELAY(rc + "COINC_DELAY", SettingType::floating, "COINC_DELAY");
  COINC_DELAY.set_val("address", 24);
  COINC_DELAY.set_val("min", 0);
  COINC_DELAY.set_val("max", 65535);
  COINC_DELAY.set_val("step", 0.013333);
  COINC_DELAY.set_val("units", "μs");
  add_definition(COINC_DELAY);

  SettingMeta CURRENT_ICR(rc + "CURRENT_ICR", SettingType::floating, "CURRENT_ICR");
  CURRENT_ICR.set_val("address", 36);
  CURRENT_ICR.set_flag("readonly");
  CURRENT_ICR.set_val("units", "count/s");
  add_definition(CURRENT_ICR);

  SettingMeta CURRENT_OORF(rc + "CURRENT_OORF", SettingType::floating, "CURRENT_OORF");
  CURRENT_OORF.set_val("address", 37);
  CURRENT_OORF.set_flag("readonly");
  CURRENT_OORF.set_val("units", "%");
  add_definition(CURRENT_OORF);

  SettingMeta EMIN(rc + "EMIN", SettingType::integer, "EMIN");
  EMIN.set_val("address", 13);
  EMIN.set_val("min", 0);
  EMIN.set_val("max", 65535);
  EMIN.set_val("step", 1);
  EMIN.set_val("description", "subtracted from computed energy in list mode");
  add_definition(EMIN);

  SettingMeta ENERGY_FLATTOP(rc + "ENERGY_FLATTOP", SettingType::floating, "ENERGY_FLATTOP");
  ENERGY_FLATTOP.set_val("address", 3);
  ENERGY_FLATTOP.set_val("min", 0);
  ENERGY_FLATTOP.set_val("max", 500);
  ENERGY_FLATTOP.set_val("step", 0.05);
  ENERGY_FLATTOP.set_val("units", "μs");
  ENERGY_FLATTOP.set_flag("optimize");
  add_definition(ENERGY_FLATTOP);

  SettingMeta ENERGY_RISETIME(rc + "ENERGY_RISETIME", SettingType::floating, "ENERGY_RISETIME");
  ENERGY_RISETIME.set_val("address", 2);
  ENERGY_RISETIME.set_val("min", 0);
  ENERGY_RISETIME.set_val("max", 500);
  ENERGY_RISETIME.set_val("step", 0.05);
  ENERGY_RISETIME.set_val("units", "μs");
  ENERGY_RISETIME.set_flag("optimize");
  add_definition(ENERGY_RISETIME);

  SettingMeta FAST_PEAKS(rc + "FAST_PEAKS", SettingType::floating, "FAST_PEAKS");
  FAST_PEAKS.set_val("address", 28);
  FAST_PEAKS.set_flag("readonly");
  add_definition(FAST_PEAKS);

  SettingMeta FTDT(rc + "FTDT", SettingType::floating, "FTDT");
  FTDT.set_val("address", 33);
  FTDT.set_val("units", "s");
  FTDT.set_flag("readonly");
  add_definition(FTDT);

  SettingMeta GATE_COUNTS(rc + "GATE_COUNTS", SettingType::floating, "GATE_COUNTS");
  GATE_COUNTS.set_val("address", 32);
  GATE_COUNTS.set_flag("readonly");
  add_definition(GATE_COUNTS);

  SettingMeta GATE_DELAY(rc + "GATE_DELAY", SettingType::floating, "GATE_DELAY");
  GATE_DELAY.set_val("address", 23);
  GATE_DELAY.set_val("min", 0.013333);
  GATE_DELAY.set_val("max", 3.4);
  GATE_DELAY.set_val("step", 0.013333);
  GATE_DELAY.set_val("units", "μs");
  add_definition(GATE_DELAY);

  SettingMeta GATE_RATE(rc + "GATE_RATE", SettingType::floating, "GATE_RATE");
  GATE_RATE.set_val("address", 31);
  GATE_RATE.set_flag("readonly");
  GATE_RATE.set_val("units", "counts/s");
  add_definition(GATE_RATE);

  SettingMeta GATE_WINDOW(rc + "GATE_WINDOW", SettingType::floating, "GATE_WINDOW");
  GATE_WINDOW.set_val("address", 22);
  GATE_WINDOW.set_val("min", 0.013333);
  GATE_WINDOW.set_val("max", 3.4);
  GATE_WINDOW.set_val("step", 0.013333);
  GATE_WINDOW.set_val("units", "μs");
  add_definition(GATE_WINDOW);

  SettingMeta GDT(rc + "GDT", SettingType::floating, "GDT");
  GDT.set_val("address", 35);
  GDT.set_flag("readonly");
  GDT.set_val("units", "s");
  add_definition(GDT);

  SettingMeta INPUT_COUNT_RATE(rc + "INPUT_COUNT_RATE", SettingType::floating, "INPUT_COUNT_RATE");
  INPUT_COUNT_RATE.set_val("address", 27);
  INPUT_COUNT_RATE.set_flag("readonly");
  INPUT_COUNT_RATE.set_val("units", "counts/s");
  add_definition(INPUT_COUNT_RATE);

  SettingMeta INTEGRATOR(rc + "INTEGRATOR", SettingType::menu, "INTEGRATOR");
  INTEGRATOR.set_enum(0, "trapezoidal");
  INTEGRATOR.set_enum(1, "gap sum");
  INTEGRATOR.set_enum(2, "step");
  INTEGRATOR.set_enum(3, "2x gap");
  INTEGRATOR.set_enum(4, "4x gap");
  INTEGRATOR.set_enum(5, "8x gap");
  add_definition(INTEGRATOR);

  SettingMeta LIVE_TIME(rc + "LIVE_TIME", SettingType::floating, "LIVE_TIME");
  LIVE_TIME.set_val("address", 26);
  LIVE_TIME.set_flag("readonly");
  LIVE_TIME.set_val("units", "s");
  add_definition(LIVE_TIME);

  SettingMeta NOUT(rc + "NOUT", SettingType::floating, "NOUT");
  NOUT.set_val("address", 30);
  NOUT.set_flag("readonly");
  add_definition(NOUT);

  SettingMeta OUTPUT_COUNT_RATE(rc + "OUTPUT_COUNT_RATE", SettingType::floating, "OUTPUT_COUNT_RATE");
  OUTPUT_COUNT_RATE.set_val("address", 29);
  OUTPUT_COUNT_RATE.set_flag("readonly");
  OUTPUT_COUNT_RATE.set_val("units", "counts/s");
  add_definition(OUTPUT_COUNT_RATE);

  SettingMeta PSA_END(rc + "PSA_END", SettingType::floating, "PSA_END");
  PSA_END.set_val("address", 12);
  PSA_END.set_val("min", 0);
  PSA_END.set_val("max", 13.1936);
  PSA_END.set_val("step", 0.013333);
  PSA_END.set_val("units", "μs");
  add_definition(PSA_END);

  SettingMeta PSA_START(rc + "PSA_START", SettingType::floating, "PSA_START");
  PSA_START.set_val("address", 11);
  PSA_START.set_val("min", 0);
  PSA_START.set_val("max", 13.1936);
  PSA_START.set_val("step", 0.013333);
  PSA_START.set_val("units", "μs");
  add_definition(PSA_START);

  SettingMeta PSM_GAIN_AVG(rc + "PSM_GAIN_AVG", SettingType::floating, "PSM_GAIN_AVG");
  PSM_GAIN_AVG.set_val("address", 38);
  PSM_GAIN_AVG.set_flag("readonly");
  add_definition(PSM_GAIN_AVG);

  SettingMeta PSM_GAIN_AVG_LEN(rc + "PSM_GAIN_AVG_LEN", SettingType::floating, "PSM_GAIN_AVG_LEN");
  PSM_GAIN_AVG_LEN.set_val("address", 39);
  PSM_GAIN_AVG_LEN.set_flag("readonly");
  add_definition(PSM_GAIN_AVG_LEN);

  // is this correct?
  SettingMeta PSM_GAIN_CORR(rc + "PSM_GAIN_CORR", SettingType::floating, "PSM_GAIN_CORR");
  PSM_GAIN_CORR.set_val("address", 42);
  PSM_GAIN_CORR.set_val("min", 0);
  PSM_GAIN_CORR.set_val("max", 1);
  PSM_GAIN_CORR.set_val("step", 0.1);
  add_definition(PSM_GAIN_CORR);

  SettingMeta PSM_TEMP_AVG(rc + "PSM_TEMP_AVG", SettingType::floating, "PSM_TEMP_AVG");
  PSM_TEMP_AVG.set_val("address", 40);
  PSM_TEMP_AVG.set_flag("readonly");
  PSM_TEMP_AVG.set_val("units", "°C");
  add_definition(PSM_TEMP_AVG);

  SettingMeta PSM_TEMP_AVG_LEN(rc + "PSM_TEMP_AVG_LEN", SettingType::floating, "PSM_TEMP_AVG_LEN");
  PSM_TEMP_AVG_LEN.set_val("address", 41);
  PSM_TEMP_AVG_LEN.set_val("min", 0);
  PSM_TEMP_AVG_LEN.set_val("max", 1);
  PSM_TEMP_AVG_LEN.set_val("step", 0.1);
  add_definition(PSM_TEMP_AVG_LEN);

  SettingMeta SFDT(rc + "SFDT", SettingType::floating, "SFDT");
  SFDT.set_val("address", 34);
  SFDT.set_flag("readonly");
  SFDT.set_val("units", "s");
  add_definition(SFDT);

  SettingMeta TAU(rc + "TAU", SettingType::floating, "TAU");
  TAU.set_val("address", 15);
  TAU.set_val("min", 0); // "1.5e-05"
  TAU.set_val("max", 65535);
  TAU.set_val("step", 1);
  TAU.set_val("units", "μs");
  add_definition(TAU);

  SettingMeta TRACE_DELAY(rc + "TRACE_DELAY", SettingType::floating, "TRACE_DELAY");
  TRACE_DELAY.set_val("address", 10);
  TRACE_DELAY.set_val("min", 0);
  TRACE_DELAY.set_val("max", 13.1936);
  TRACE_DELAY.set_val("step", 0.013333);
  TRACE_DELAY.set_val("units", "μs");
  add_definition(TRACE_DELAY);

  SettingMeta TRACE_LENGTH(rc + "TRACE_LENGTH", SettingType::floating, "TRACE_LENGTH");
  TRACE_LENGTH.set_val("address", 9);
  TRACE_LENGTH.set_val("min", 0);
  TRACE_LENGTH.set_val("max", 13.1936);
  TRACE_LENGTH.set_val("step", 0.013333);
  TRACE_LENGTH.set_val("units", "μs");
  add_definition(TRACE_LENGTH);

  SettingMeta TRIGGER_FLATTOP(rc + "TRIGGER_FLATTOP", SettingType::floating, "TRIGGER_FLATTOP");
  TRIGGER_FLATTOP.set_val("address", 5);
  TRIGGER_FLATTOP.set_val("min", 0);
  TRIGGER_FLATTOP.set_val("max", 500);
  TRIGGER_FLATTOP.set_val("step", 0.05);
  TRIGGER_FLATTOP.set_val("units", "μs");
  add_definition(TRIGGER_FLATTOP);

  SettingMeta TRIGGER_RISETIME(rc + "TRIGGER_RISETIME", SettingType::floating, "TRIGGER_RISETIME");
  TRIGGER_RISETIME.set_val("address", 4);
  TRIGGER_RISETIME.set_val("min", 0);
  TRIGGER_RISETIME.set_val("max", 500);
  TRIGGER_RISETIME.set_val("step", 0.05);
  TRIGGER_RISETIME.set_val("units", "μs");
  add_definition(TRIGGER_RISETIME);

  SettingMeta TRIGGER_THRESHOLD(rc + "TRIGGER_THRESHOLD", SettingType::integer, "TRIGGER_THRESHOLD");
  TRIGGER_THRESHOLD.set_val("address", 6);
  TRIGGER_THRESHOLD.set_val("min", 0);
  TRIGGER_THRESHOLD.set_val("max", 4095);
  TRIGGER_THRESHOLD.set_val("step", 1);
  add_definition(TRIGGER_THRESHOLD);

  SettingMeta VGAIN(rc + "VGAIN", SettingType::floating, "VGAIN");
  VGAIN.set_val("address", 7);
  VGAIN.set_val("min", 0);
  VGAIN.set_val("max", 65535);
  VGAIN.set_val("step", 0.000076);
  VGAIN.set_val("units", "V/V");
  VGAIN.set_flag("gain");
  add_definition(VGAIN);

  SettingMeta VOFFSET(rc + "VOFFSET", SettingType::floating, "VOFFSET");
  VOFFSET.set_val("address", 8);
  VOFFSET.set_val("min", -2.5);
  VOFFSET.set_val("max", 2.5);
  VOFFSET.set_val("step", 0.000076);
  VOFFSET.set_val("units", "V");
  add_definition(VOFFSET);

  SettingMeta XDT(rc + "XDT", SettingType::floating, "XDT");
  XDT.set_val("address", 17);
  XDT.set_val("min", 0.053333);
  XDT.set_val("max", 873.77333);
  XDT.set_val("step", 0.013333);
  XDT.set_val("units", "μs");
  add_definition(XDT);

  int32_t i{0};
  SettingMeta MODULE_CHANNEL(r + "module/channel", SettingType::stem);
  MODULE_CHANNEL.set_flag("saveworthy");
  MODULE_CHANNEL.set_enum(0, rc + "CHANNEL_CSRA");
  MODULE_CHANNEL.set_enum(1, rc + "CHANNEL_CSRB");
  MODULE_CHANNEL.set_enum(2, rc + "ENERGY_RISETIME");
  MODULE_CHANNEL.set_enum(3, rc + "ENERGY_FLATTOP");
  MODULE_CHANNEL.set_enum(4, rc + "TRIGGER_RISETIME");
  MODULE_CHANNEL.set_enum(5, rc + "TRIGGER_FLATTOP");
  MODULE_CHANNEL.set_enum(6, rc + "TRIGGER_THRESHOLD");
  MODULE_CHANNEL.set_enum(7, rc + "VGAIN");
  MODULE_CHANNEL.set_enum(8, rc + "VOFFSET");
  MODULE_CHANNEL.set_enum(9, rc + "TRACE_LENGTH");
  MODULE_CHANNEL.set_enum(10, rc + "TRACE_DELAY");
  MODULE_CHANNEL.set_enum(11, rc + "PSA_START");
  MODULE_CHANNEL.set_enum(12, rc + "PSA_END");
  MODULE_CHANNEL.set_enum(13, rc + "EMIN");
  MODULE_CHANNEL.set_enum(14, rc + "BINFACTOR");
  MODULE_CHANNEL.set_enum(15, rc + "TAU");
  MODULE_CHANNEL.set_enum(16, rc + "BLCUT");
  MODULE_CHANNEL.set_enum(17, rc + "XDT");
  MODULE_CHANNEL.set_enum(18, rc + "BASELINE_PERCENT");
  MODULE_CHANNEL.set_enum(19, rc + "CFD_THRESHOLD");
  MODULE_CHANNEL.set_enum(20, rc + "INTEGRATOR");
  MODULE_CHANNEL.set_enum(21, rc + "CHANNEL_CSRC");
  MODULE_CHANNEL.set_enum(22, rc + "GATE_WINDOW");
  MODULE_CHANNEL.set_enum(23, rc + "GATE_DELAY");
  MODULE_CHANNEL.set_enum(24, rc + "COINC_DELAY");
  MODULE_CHANNEL.set_enum(25, rc + "BLAVG");
  MODULE_CHANNEL.set_enum(26, rc + "LIVE_TIME");
  MODULE_CHANNEL.set_enum(27, rc + "INPUT_COUNT_RATE");
  MODULE_CHANNEL.set_enum(28, rc + "FAST_PEAKS");
  MODULE_CHANNEL.set_enum(29, rc + "OUTPUT_COUNT_RATE");
  MODULE_CHANNEL.set_enum(30, rc + "NOUT");
  MODULE_CHANNEL.set_enum(31, rc + "GATE_RATE");
  MODULE_CHANNEL.set_enum(32, rc + "GATE_COUNTS");
  MODULE_CHANNEL.set_enum(33, rc + "FTDT");
  MODULE_CHANNEL.set_enum(34, rc + "SFDT");
  MODULE_CHANNEL.set_enum(35, rc + "GDT");
  MODULE_CHANNEL.set_enum(36, rc + "CURRENT_ICR");
  MODULE_CHANNEL.set_enum(37, rc + "CURRENT_OORF");
  MODULE_CHANNEL.set_enum(38, rc + "PSM_GAIN_AVG");
  MODULE_CHANNEL.set_enum(39, rc + "PSM_GAIN_AVG_LEN");
  MODULE_CHANNEL.set_enum(40, rc + "PSM_TEMP_AVG");
  MODULE_CHANNEL.set_enum(41, rc + "PSM_TEMP_AVG_LEN");
  MODULE_CHANNEL.set_enum(42, rc + "PSM_GAIN_CORR");
  add_definition(MODULE_CHANNEL);

  SettingMeta ACTUAL_COINCIDENCE_WAIT(rm + "ACTUAL_COINCIDENCE_WAIT", SettingType::floating, "ACTUAL_COINCIDENCE_WAIT");
  ACTUAL_COINCIDENCE_WAIT.set_val("address", 6);
  ACTUAL_COINCIDENCE_WAIT.set_val("min", 13.333333);
  ACTUAL_COINCIDENCE_WAIT.set_val("max", 87374663.6448);
  ACTUAL_COINCIDENCE_WAIT.set_val("step", 13.333333);
  ACTUAL_COINCIDENCE_WAIT.set_val("units", "ns");
  add_definition(ACTUAL_COINCIDENCE_WAIT);


  /*
      <SettingMeta id="Pixie4/System/module/BOARD_VERSION" type="floating" name="BOARD_VERSION" address="24" writable="false" step="0" minimum="0" maximum="0" />

          <SettingMeta id="Pixie4/System/module/BUFFER_HEAD_LENGTH" type="integer" name="BUFFER_HEAD_LENGTH" address="16" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/CHANNEL_HEAD_LENGTH" type="integer" name="CHANNEL_HEAD_LENGTH" address="18" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/COINCIDENCE_PATTERN" type="binary" name="COINCIDENCE_PATTERN" address="5" writable="true" word_size="16">
      <flag bit="0" description="oooo" />
      <flag bit="1" description="ooo+" />
      <flag bit="2" description="oo+o" />
      <flag bit="3" description="oo++" />
      <flag bit="4" description="o+oo" />
      <flag bit="5" description="o+o+" />
      <flag bit="6" description="o++o" />
      <flag bit="7" description="o+++" />
      <flag bit="8" description="+ooo" />
      <flag bit="9" description="+oo+" />
      <flag bit="10" description="+o+o" />
      <flag bit="11" description="+o++" />
      <flag bit="12" description="++oo" />
      <flag bit="13" description="++o+" />
      <flag bit="14" description="+++o" />
      <flag bit="15" description="++++" />
      </SettingMeta>
      <SettingMeta id="Pixie4/System/module/DBLBUFCSR" type="binary" name="DBLBUFCSR" address="14" writable="true" word_size="16">
      <flag bit="0" description="Enable double buffer mode" />
      <flag bit="1" description="Host has read buffer" />
      <flag bit="3" description="Host should read first block, else should read second block" />
      </SettingMeta>
      <SettingMeta id="Pixie4/System/module/DSP_BUILD" type="floating" name="DSP_BUILD" address="27" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/DSP_RELEASE" type="floating" name="DSP_RELEASE" address="26" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/EVENT_HEAD_LENGTH" type="integer" name="EVENT_HEAD_LENGTH" address="17" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/EVENT_RATE" type="floating" name="EVENT_RATE" address="22" writable="false" step="0" minimum="0" maximum="0" unit="cps" />
      <SettingMeta id="Pixie4/System/module/FILTER_RANGE" type="int_menu" name="FILTER_RANGE" address="11" writable="true">
      <menu_item item_value="1" item_text="2 bits" />
      <menu_item item_value="2" item_text="4 bits" />
      <menu_item item_value="3" item_text="8 bits" />
      <menu_item item_value="4" item_text="16 bits" />
      <menu_item item_value="5" item_text="32 bits" />
      <menu_item item_value="6" item_text="64 bits" />
      </SettingMeta>
      <SettingMeta id="Pixie4/System/module/FIPPI_ID" type="floating" name="FIPPI_ID" address="28" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/IN_SYNCH" type="boolean" name="IN_SYNCH" address="9" writable="true" description="modules synchronized; clear to have system reset clocks at start of daq" />
      <SettingMeta id="Pixie4/System/module/MAX_EVENTS" type="integer" name="MAX_EVENTS" address="4" writable="true" step="1" minimum="0" maximum="100" />
      <SettingMeta id="Pixie4/System/module/MIN_COINCIDENCE_WAIT" type="integer" name="MIN_COINCIDENCE_WAIT" address="7" writable="false" step="0" minimum="0" maximum="0" unit="ticks" />
      <SettingMeta id="Pixie4/System/module/MODULEPATTERN" type="binary" name="MODULEPATTERN" address="12" writable="true" word_size="16">
      <flag bit="4" description="Gate event acceptance on FRONT panel input" />
      <flag bit="5" description="Gate event acceptance on LOCAL coincidence test" />
      <flag bit="6" description="Gate event acceptance on backplane STATUS line" />
      <flag bit="7" description="Gate event acceptance on GLOBAL coincidence test" />
      </SettingMeta>
      <SettingMeta id="Pixie4/System/module/MODULE_CSRA" type="binary" name="MODULE_CSRA" address="1" writable="true" word_size="16">
      <flag bit="1" description="acquire 32 buffers and write to EM, else only one at a time" />
      <flag bit="2" description="backpane trigger distribution: see manual..." />
      <flag bit="3" description="Pixie native MCA: bin sums to addback spectrum" />
      <flag bit="4" description="Pixie native MCA: individual spectra only singles" />
      <flag bit="5" description="front panel DSP-OUT distributed as veto signal to backplane" />
      <flag bit="6" description="chan3 hit status contributes to backplane STATUS" />
      <flag bit="7" description="polarity of front panel pulse counter" />
      <flag bit="9" description="write NNSHAREPATTERN to left neighbor (should be PDM) during ControlTask5" />
      <flag bit="10" description="Pixie-500 only: time stamps as 2ns, else 13.3ns" />
      <flag bit="12" description="drive low TOKEN backplane if local coincidence fails" />
      <flag bit="13" description="send hit patten to slot2 using PXI STAR trigger for each event" />
      <flag bit="14" description="front panel DSP-OUT as input to STATUS on backplane (wire-OR)" />
      <flag bit="15" description="backpane trigger distribution: see manual..." />
      </SettingMeta>
      <SettingMeta id="Pixie4/System/module/MODULE_CSRB" type="binary" name="MODULE_CSRB" address="2" writable="true" word_size="16">
      <flag bit="0" description="Execute user code routines programmed by user dsp" />
      </SettingMeta>
      <SettingMeta id="Pixie4/System/module/MODULE_CSRC" type="binary" name="MODULE_CSRC" address="15" writable="true" word_size="16" />
      <SettingMeta id="Pixie4/System/module/MODULE_FORMAT" type="binary" name="MODULE_FORMAT" address="3" writable="false" word_size="16" description="not used" />
      <SettingMeta id="Pixie4/System/module/MODULE_NUMBER" type="integer" name="MODULE_NUMBER" address="0" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/NNSHAREPATTERN" type="integer" name="NNSHAREPATTERN" address="13" writable="true" step="1" minimum="0" maximum="65535" description="User-defined control word for PXI-PDM" />
      <SettingMeta id="Pixie4/System/module/NUMBER_EVENTS" type="integer" name="NUMBER_EVENTS" address="20" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/OUTPUT_BUFFER_LENGTH" type="integer" name="OUTPUT_BUFFER_LENGTH" address="19" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/PDM_MASKA" type="integer" name="PDM_MASKA" address="31" writable="true" step="1" minimum="0" maximum="65535" />
      <SettingMeta id="Pixie4/System/module/PDM_MASKB" type="integer" name="PDM_MASKB" address="32" writable="true" step="1" minimum="0" maximum="65535" />
      <SettingMeta id="Pixie4/System/module/PDM_MASKC" type="integer" name="PDM_MASKC" address="33" writable="true" step="1" minimum="0" maximum="65535" />
      <SettingMeta id="Pixie4/System/module/RUN_TIME" type="floating" name="RUN_TIME" address="21" writable="false" step="0" minimum="0" maximum="0" unit="s" />
      <SettingMeta id="Pixie4/System/module/RUN_TYPE" type="int_menu" name="RUN_TYPE" address="10" writable="true">
      <menu_item item_value="0" item_text="Slow control run" />
      <menu_item item_value="256" item_text="Traces" />
      <menu_item item_value="257" item_text="Full" />
      <menu_item item_value="258" item_text="PSA only" />
      <menu_item item_value="259" item_text="Compressed" />
      <menu_item item_value="769" item_text="Pixie MCA" />
      </SettingMeta>
      <SettingMeta id="Pixie4/System/module/SERIAL_NUMBER" type="floating" name="SERIAL_NUMBER" address="25" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/SYNCH_WAIT" type="boolean" name="SYNCH_WAIT" address="8" writable="true" description="wait for all modules ready before starting daq" />
      <SettingMeta id="Pixie4/System/module/SYSTEM_ID" type="floating" name="SYSTEM_ID" address="29" writable="false" step="0" minimum="0" maximum="0" />
      <SettingMeta id="Pixie4/System/module/TOTAL_TIME" type="floating" name="TOTAL_TIME" address="23" writable="false" step="0" minimum="0" maximum="0" unit="s" />
      <SettingMeta id="Pixie4/System/module/XET_DELAY" type="integer" name="XET_DELAY" address="30" writable="true" step="1" minimum="0" maximum="65535" description="delay for generated event trigger from front panel to backplane" />
*/

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
