#pragma once

#include <core/plugin/container.h>
#include <core/gamma/nuclides/transition.h>

namespace DAQuiri {

class Isotope {
 public:
  Isotope () {}
  Isotope (std::string nm) : Isotope() {name = nm;}

  bool shallow_equals(const Isotope& other) const {return (name == other.name);}
  bool operator!= (const Isotope& other) const {return !(this->operator==(other));}
  bool operator== (const Isotope& other) const {
    if (name != other.name) return false;
    if (half_life != other.half_life) return false;
    if (gammas != other.gammas) return false;
    return true;
  }

  std::string name;
  std::string comment;
  UncertainDouble half_life;
  std::string half_life_comment;
  Container<Radiation> gammas;
};

void to_json(nlohmann::json& j, const Isotope& s);
void from_json(const nlohmann::json& j, Isotope& s);

}
