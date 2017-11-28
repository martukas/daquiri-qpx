#pragma once

#include <string>
#include <cmath>
#include "precise_float.h"

#include "json.hpp"
using namespace nlohmann;

namespace DAQuiri {

class TimeBase
{
private:
  PreciseFloat multiplier_ {1};
  PreciseFloat divider_ {1};

public:
  inline TimeBase() {}

  inline TimeBase(uint32_t multiplier, uint32_t divider)
  {
    if (!multiplier)
      throw std::domain_error("multiplier must be non-0");
    if (!divider)
      throw std::domain_error("divider must be non-0");
    multiplier_ = multiplier;
    divider_ =  divider;
    auto GCD = gcd(multiplier_, divider_);
    if (GCD > 1)
    {
      multiplier_ /= GCD;
      divider_ /= GCD;
    }
  }

  inline PreciseFloat multiplier() const
  {
    return multiplier_;
  }

  inline PreciseFloat divider() const
  {
    return divider_;
  }

  inline PreciseFloat to_nanosec(PreciseFloat native) const
  {
    return native * multiplier_ / divider_;
  }

  inline PreciseFloat to_microsec(PreciseFloat native) const
  {
    return native * multiplier_ / divider_ * 0.001;
  }

  inline PreciseFloat to_millisec(PreciseFloat native) const
  {
    return native * multiplier_ / divider_ * 0.000001;
  }

  inline PreciseFloat to_sec(PreciseFloat native) const
  {
    return native * multiplier_ / divider_ * 0.000000001;
  }

  inline PreciseFloat to_native(PreciseFloat ns)
  {
    return ns * divider_ / multiplier_;
  }

  inline bool operator==(const TimeBase& other) const
  {
    return ((divider_ == other.divider_)
            && (multiplier_ == other.multiplier_));
  }

  inline bool operator!=(const TimeBase& other) const
  {
    return (!operator ==(other));
  }

  static inline TimeBase common(const TimeBase& a, const TimeBase& b)
  {
    if (a == b)
      return a;
    return TimeBase(gcd(a.multiplier_, b.multiplier_),
                    lcm(a.divider_, b.divider_));
  }

  inline std::string debug() const
  {
    std::stringstream ss;
    if ((multiplier_ != 1) && (divider_ != 1))
      ss << "x(" << multiplier_ << "/" << divider_ << ")";
    return ss.str();
  }

private:
  static inline uint32_t lcm(uint32_t a, uint32_t b)
  {
    uint32_t m(a), n(b);
    while(m!=n)
    {
      if(m < n)
        m = m + a;
      else
        n = n + b;
    }
    return m;
  }

  static inline uint32_t gcd(uint32_t a, uint32_t b)
  {
    return b == 0 ? a : gcd(b, a % b);
  }
};

inline void to_json(json& j, const TimeBase& t)
{
  j["multiplier"] = t.multiplier();
  j["divider"] = t.divider();
}

inline void from_json(const json& j, TimeBase& t)
{
  t = TimeBase(j["multiplier"], j["divider"]);
}

}