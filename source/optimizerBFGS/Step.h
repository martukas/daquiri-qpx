#pragma once

#include <optimizerBFGS/Value.h>

namespace Hypermet
{

struct Step
{
  Step() = default;
  Step(Side s) : side(s) {}

  bool override{false};
  bool enabled{true};
  Value amplitude;
  Side side {Side::left};

  inline double flip(double spread) const
  {
    if (side == Side::right)
      return -spread;
    return spread;
  }

  double eval_with(const PrecalcVals& pre, double ampl) const;

  double eval(const PrecalcVals& pre) const;

  double eval_grad(const PrecalcVals& pre,
                   std::vector<double>& grads,
                   size_t i_width, size_t i_pos, size_t i_amp) const;

  void update_indices(int32_t& i);
  void put(std::vector<double>& fit) const;
  void get(const std::vector<double>& fit);
  void get_uncerts(const std::vector<double>& diagonals, double chisq_norm);
};

}
