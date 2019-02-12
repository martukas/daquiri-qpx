#include <core/fitting/hypermet/GaussianPrecalc.h>

namespace DAQuiri
{

std::string to_string(const Side& s)
{
  if (s == Side::right)
    return "right";
  else
    return "left";
}

Side to_side(const std::string& s)
{
  if (s == "right")
    return Side::right;
  else
    return Side::left;
}


}
