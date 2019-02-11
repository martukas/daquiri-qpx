#include "gtest_color_print.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Peak.h>

class Peak : public TestBase
{
 protected:

  virtual void SetUp()
  {
    peak = peak.gaussian_only();
    peak.amplitude.bound(0, 500);
    peak.amplitude.val(400);
    peak.position.bound(5, 15);
    peak.position.val(10);
    peak.width_override = true;
    peak.width_.val(3.2);
  }

  DAQuiri::WeightedData generate_data()
  {
    std::vector<double> channels;
    std::vector<double> y;
    for (size_t i = 0; i < 20; ++i)
    {
      channels.push_back(i);
      y.push_back(peak.eval(i).all());
    }
    return DAQuiri::WeightedData(channels, y);
  }

  void grad_at_point(DAQuiri::Value& variable, size_t channel_idx)
  {
    auto wdata = generate_data();

    int32_t idx {0};
    peak.update_indices(idx);
    EXPECT_EQ(idx, 3);

    double degrees_freedom = wdata.data.size() - idx;
    EXPECT_EQ(degrees_freedom, 17);

    auto chosen_point = wdata.data.at(channel_idx);
    EXPECT_EQ(chosen_point.x, channel_idx);
    EXPECT_EQ(chosen_point.y, peak.eval(chosen_point.x).all());

    size_t chosen_var_idx = variable.index();

    std::vector<double> x_val;
    std::vector<double> val_val;
    std::vector<double> chi_sq_norm;
    std::vector<double> gradient;

    Eigen::VectorXd chan_gradients;
    for (double val_x = -4; val_x < 4; val_x += 0.2)
    {
      x_val.push_back(val_x);
      val_val.push_back(variable.val_at(val_x));

      variable.x(val_x);
      chan_gradients.setConstant(idx, 0.0);
      double FTotal = peak.eval_grad(chosen_point.x, chan_gradients).all();
      double t3 = -2.0 * (chosen_point.y - FTotal) / square(chosen_point.weight_phillips_marlow);

      gradient.push_back(chan_gradients[chosen_var_idx] * t3);
      chi_sq_norm.push_back(
          square((chosen_point.y - FTotal) / chosen_point.weight_phillips_marlow)
              / degrees_freedom);
    }

    MESSAGE() << "chi_sq_norm(x):\n" << visualize(x_val, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(x):\n" << visualize(x_val, gradient, 100) << "\n";
    MESSAGE() << "sq_norm(val):\n" << visualize(val_val, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(val):\n" << visualize(val_val, gradient, 100) << "\n";
  }

  void grad(DAQuiri::Value& variable)
  {
    auto wdata = generate_data();

    int32_t idx {0};
    peak.update_indices(idx);
    ASSERT_EQ(idx, 3);

    double degrees_freedom = wdata.data.size() - idx;
    ASSERT_EQ(degrees_freedom, 17);

    size_t chosen_var_idx = variable.index();

    std::vector<double> x_val;
    std::vector<double> val_val;
    std::vector<double> chi_sq_norm;
    std::vector<double> gradient;

    Eigen::VectorXd gradients;
    Eigen::VectorXd chan_gradients;
    for (double val_x = -4; val_x < 4; val_x += 0.2)
    {
      x_val.push_back(val_x);
      val_val.push_back(variable.val_at(val_x));
      variable.x(val_x);

      gradients.setConstant(idx, 0.0);
      Eigen::VectorXd chan_gradients;
      double Chisq{0.0};
      for (const auto& data : wdata.data)
      {
        chan_gradients.setConstant(idx, 0.0);

        double FTotal = peak.eval_grad(data.x, chan_gradients).all();

        double t3 = -2.0 * (data.y - FTotal) / square(data.weight_phillips_marlow);

        for (size_t var = 0; var < static_cast<size_t>(idx); ++var)
          gradients[var] += chan_gradients[var] * t3;
        Chisq += square((data.y - FTotal) / data.weight_phillips_marlow);
      }

      gradient.push_back(gradients[chosen_var_idx]);
      chi_sq_norm.push_back(Chisq / degrees_freedom);
    }


    MESSAGE() << "chi_sq_norm(x):\n" << visualize(x_val, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(x):\n" << visualize(x_val, gradient, 100) << "\n";
    MESSAGE() << "chi_sq_norm(val):\n" << visualize(val_val, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(val):\n" << visualize(val_val, gradient, 100) << "\n";
  }

  DAQuiri::Peak peak;
};

TEST_F(Peak, CheckPeakSetup)
{
  MESSAGE() << "Peak:\n" << peak.to_string(" ") << "\n";
}

TEST_F(Peak, EvalGaussianOnly)
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 20; ++i)
  {
    channels.push_back(i);
    y.push_back(peak.eval(i).all());
  }
  MESSAGE() << "peak(channel):\n" << visualize(channels, y, 100) << "\n";
}

TEST_F(Peak, GradAmplitudeAtOnePoint)
{
  grad_at_point(peak.amplitude, 10);
}

TEST_F(Peak, GradAmplitudeOnly)
{
  grad(peak.amplitude);
}

TEST_F(Peak, GradWidthAtOnePoint)
{
  grad_at_point(peak.width_, 5);
}

TEST_F(Peak, GradWidthOnly)
{
  grad(peak.width_);
}

TEST_F(Peak, GradPositionAtOnePoint)
{
  grad_at_point(peak.position, 10);
}

TEST_F(Peak, GradPositionOnly)
{
  grad(peak.position);
}
