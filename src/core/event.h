#pragma once

#include "event_model.h"
#include <vector>

namespace DAQuiri {

class Event
{
private:
  int16_t       source_channel_ {-1};
  TimeStamp     timestamp_;
  std::vector<DigitizedVal>          values_;
  std::vector<std::vector<uint16_t>> traces_;

public:
  inline Event() {}

  inline Event(int16_t sourcechan, const EventModel &model)
    : source_channel_(sourcechan)
    , timestamp_(0, model.timebase)
    , values_ (model.values)
  {
    for (auto t : model.traces)
    {
      size_t product = 1;
      for (auto d : t)
        product *= d;
      traces_.push_back(std::vector<uint16_t>(product, 0));
    }
  }

  //Accessors
  inline const int16_t& channel() const
  {
    return source_channel_;
  }

  inline const TimeStamp& timestamp() const
  {
    return timestamp_;
  }

  inline size_t value_count() const
  {
    return values_.size();
  }

  inline size_t trace_count() const
  {
    return traces_.size();
  }

  inline DigitizedVal value(size_t idx) const
  {
    if (idx >= values_.size())
      throw std::out_of_range("Event: bad value index");
    return values_.at(idx);
  }

  inline const std::vector<uint16_t>& trace(size_t idx) const
  {
    if (idx >= traces_.size())
      throw std::out_of_range("Event: bad trace index");
    return traces_.at(idx);
  }

  //Setters
  inline void set_native_time(uint64_t t)
  {
    timestamp_.set_native(t);
  }

  inline void set_timestamp(const TimeStamp& ts)
  {
    timestamp_ = ts;
  }

  inline void delay_ns(double ns)
  {
    timestamp_.delay(ns);
  }

  inline void set_value(size_t idx, uint16_t val)
  {
    if (idx >= values_.size())
      throw std::out_of_range("Event: bad value index");
    values_[idx].set_val(val);
  }

  inline void set_trace(size_t idx, const std::vector<uint16_t> &trc)
  {
    if (idx >= traces_.size())
      throw std::out_of_range("Event: bad trace index");
    auto& t = traces_[idx];
    size_t len = std::min(trc.size(), t.size());
    for (size_t i=0; i < len; ++i)
      t[i] = trc.at(i);
  }

  //Comparators
  inline bool operator==(const Event other) const
  {
    if (source_channel_ != other.source_channel_) return false;
    if (timestamp_ != other.timestamp_) return false;
    if (values_ != other.values_) return false;
    if (traces_ != other.traces_) return false;
    return true;
  }

  inline bool operator!=(const Event other) const
  {
    return !operator==(other);
  }

  inline std::string debug() const
  {
    std::stringstream ss;
    ss << "[ch" << source_channel_ << "|t" << timestamp_.debug();
    if (traces_.size())
      ss << "|ntraces=" << traces_.size();
    for (auto &v : values_)
      ss << " " << v.debug();
    ss << "]";
    return ss.str();
  }
};

}
