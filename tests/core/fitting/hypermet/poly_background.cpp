#include "gtest_color_print.h"
#include "visualize_vector.h"

#include <core/fitting/hypermet/PolyBackground.h>

class PolyBackground : public TestBase
{
};


TEST_F(PolyBackground, val_grad)
{
  DAQuiri::PolyBackground pb;
  pb.x_offset = 10;
  pb.base.val(70);
  pb.slope.val(3);
  pb.curve.val(5);
  int32_t k{0};
  pb.update_indices(k);

  std::vector<double> x;
  std::vector<double> base, slope, curve;
  Eigen::VectorXd grad(3);
  for (size_t i=0; i < 30; ++i)
  {
    x.push_back(i);
    pb.eval_grad(i, grad);
    base.push_back(grad[0]);
    slope.push_back(grad[1]);
    curve.push_back(grad[2]);
  }

  auto y = pb.eval(x);

  MESSAGE() << "Val:\n" << visualize(x, y, 80) << "\n";
  MESSAGE() << "Grad(base):\n" << visualize(x, base, 80) << "\n";
  MESSAGE() << "Grad(slope):\n" << visualize(x, slope, 80) << "\n";
  MESSAGE() << "Grad(curve):\n" << visualize(x, curve, 80) << "\n";
}
