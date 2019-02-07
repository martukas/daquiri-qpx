#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

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

  MESSAGE() << "PolyBackground:\n" << pb.to_string(" ") << "\n";
  MESSAGE() << "val(chan):\n" << visualize(x, y, 80) << "\n";
  MESSAGE() << "base.grad(chan):\n" << visualize(x, base, 80) << "\n";
  MESSAGE() << "slope.grad(chan):\n" << visualize(x, slope, 80) << "\n";
  MESSAGE() << "curve.grad(chan):\n" << visualize(x, curve, 80) << "\n";
}
