#pragma once

#include <core/fitting/uncertain.h>
#include <string>

namespace DAQuiri {

class Radiation {
 public:
  Radiation() = default;
  Radiation(UncertainDouble en, UncertainDouble ab)
    : energy(en), abundance(ab) {}

  bool shallow_equals(const Radiation& other) const {return (energy == other.energy);}
  bool operator!= (const Radiation& other) const {return !(this->operator==(other));}
  bool operator== (const Radiation& other) const {
    if (energy != other.energy) return false;
    if (abundance != other.abundance) return false;
    return true;
  }

  UncertainDouble energy;
  UncertainDouble abundance;
  std::string comment;
  bool marked {false};
};

void to_json(nlohmann::json& j, const Radiation& s);
void from_json(const nlohmann::json& j, Radiation& s);

}
