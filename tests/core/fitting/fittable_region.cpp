#include "hypermet/function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/fittable_region.h>
#include <core/fitting/hypermet/Value.h>

class FittableLine : public DAQuiri::FittableRegion
{
 public:
  DAQuiri::Value val;

  double eval(double chan) const override
  {
    return val.val() * chan;
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return val.val_from(fit) * chan;
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    double ret = val.val_from(fit) * chan;
    if (val.index() >= 0)
      grads[val.index()] = val.grad_from(fit) * chan;
    return ret;
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    val.put(ret);
    return ret;
  }
};

class FittableRegion : public FunctionTest
{
 protected:
  FittableLine fl;

  virtual void SetUp()
  {
    fl.val.bound(-10, 10);
    fl.val.val(2.5);
  }
};

TEST_F(FittableRegion, DefaultConstructed)
{
  EXPECT_EQ(fl.variable_count, 0);
  //fl.val.update_index(fl.variable_count);
}
