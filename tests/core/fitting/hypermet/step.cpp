#include "gtest_color_print.h"

#include <dlib/matrix/matrix_mat.h>

#include <core/fitting/hypermet/Step.h>

class Step : public TestBase
{
  DAQuiri::Step step;

  virtual void SetUp()
  {
    step.amplitude.bound(0.000001, 0.05);
    step.amplitude.to_fit = true;
  }

  double diff()
  {
    DAQuiri::PrecalcVals pre;
    pre.width = 2;
    pre.ampl = 6;
    pre.half_ampl = 3;
    pre.spread = 10;
//    const auto& fe = std::bind(&Region::eval, region_, std::placeholders::_1);
//    double length(derivative(rosen)(starting_point) - rosen_derivative(starting_point));
  }
};

TEST_F(Step, GamValueTrue)
{

}

