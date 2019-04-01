#pragma once

#include <cinttypes>
#include <string>

namespace DAQuiri
{

struct PrecalcVals
{
  double width;
  double ampl;
  double half_ampl;
  double spread;

  double width_grad;
  double pos_grad;
  double amp_grad;

  int32_t i_width;
  int32_t i_pos;
  int32_t i_amp;
};

enum class Side
{
  left,
  right
};

inline double flip(const Side& side, double spread)
{
  if (side == Side::right)
    return -spread;
  return spread;
}

std::string to_string(const Side& s);
Side to_side(const std::string& s);

}
