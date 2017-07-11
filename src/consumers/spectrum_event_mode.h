#pragma once

#include "spectrum.h"
#include "coincidence.h"

class SpectrumEventMode : virtual public Spectrum
{
public:
  SpectrumEventMode();

protected:
  bool _initialize() override;
  void _init_from_file(std::string name) override;

  void _push_event(const Event&) override;
  void _push_stats(const Status&) override;
  void _flush() override;

  bool channel_relevant(int16_t channel) const override;

  virtual bool event_relevant(const Event&) const;
  virtual bool validate_coincidence(const Coincidence&) const;
  virtual void add_coincidence(const Coincidence&) = 0;

protected:
  //cached parameters
  Pattern pattern_coinc_, pattern_anti_, pattern_add_;
  std::vector<int32_t> cutoff_logic_;
  std::vector<double>  delay_ns_;
  double max_delay_ {0};
  double coinc_window_ {0};

  //from status manifest
  std::vector<int> energy_idx_;

  //data:
  std::list<Coincidence> backlog;
  PreciseFloat total_coincidences_;
};
