#include <gui/analysis/isotope.h>

namespace RadTypes {

void to_json(nlohmann::json& j, const Radiation& s)
{
  j["energy"] = s.energy;
  j["abundance"] = s.abundance;
}

void from_json(const nlohmann::json& j, Radiation& s)
{
  s.energy = j["energy"];
  s.abundance = j["abundance"];
}

void to_json(nlohmann::json& j, const Isotope& s)
{
  j["name"] = s.name;
  j["half_life"] = s.half_life;
  j["gamma_constant"] = s.gamma_constant;
  j["beta"] = s.beta;
  if (s.gammas.size())
    j["gammas"] = s.gammas;
}

void from_json(const nlohmann::json& j, Isotope& s)
{
  s.name = j["name"];
  s.half_life = j["half_life"];
  s.gamma_constant = j["gamma_constant"];
  s.beta = j["beta"];
  if (j.count("gammas"))
    s.gammas = j["gammas"];
}


}
