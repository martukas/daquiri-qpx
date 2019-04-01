#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/gamma/fit_evaluation.h>

class FitEvaluation : public TestBase
{
};


TEST_F(FitEvaluation, DefaultConstructed)
{
  DAQuiri::FitEvaluation fe;
  EXPECT_TRUE(fe.empty());
}
