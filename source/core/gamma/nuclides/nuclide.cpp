#include <core/gamma/nuclides/nuclide.h>

namespace DAQuiri {

void to_json(nlohmann::json& j, const Isotope& s)
{
  j["name"] = s.name;
  j["half_life"] = s.half_life;
  if (s.gammas.size())
    j["gammas"] = s.gammas;
}

void from_json(const nlohmann::json& j, Isotope& s)
{
  s.name = j["name"];
  s.half_life = j["half_life"];
  if (j.count("gammas"))
    s.gammas = j["gammas"];
}


}
