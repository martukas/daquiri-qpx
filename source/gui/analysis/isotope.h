#pragma once

#include <core/plugin/container.h>
#include <string>

namespace RadTypes {

class Radiation {
 public:
  Radiation() : energy(0), abundance(0), marked(false) {}
  Radiation(double en, double ab) : energy(en), abundance(ab), marked(false) {}

  bool shallow_equals(const Radiation& other) const {return (energy == other.energy);}
  bool operator!= (const Radiation& other) const {return !(this->operator==(other));}
  bool operator== (const Radiation& other) const {
    if (energy != other.energy) return false;
    if (abundance != other.abundance) return false;
    return true;
  }

  double energy;
  double abundance;
  bool marked;
};

void to_json(nlohmann::json& j, const Radiation& s);
void from_json(const nlohmann::json& j, Radiation& s);

class Isotope {
 public:
  Isotope () : half_life(0.0) {}
  Isotope (std::string nm) : Isotope() {name = nm;}

  bool shallow_equals(const Isotope& other) const {return (name == other.name);}
  bool operator!= (const Isotope& other) const {return !(this->operator==(other));}
  bool operator== (const Isotope& other) const {
    if (name != other.name) return false;
    if (half_life != other.half_life) return false;
    if (gamma_constant != other.gamma_constant) return false;
    if (beta != other.beta) return false;
    if (gammas != other.gammas) return false;
    return true;
  }

  std::string name;
  double half_life;
  std::string gamma_constant;
  Radiation beta;
  Container<Radiation> gammas;
};

void to_json(nlohmann::json& j, const Isotope& s);
void from_json(const nlohmann::json& j, Isotope& s);

}
