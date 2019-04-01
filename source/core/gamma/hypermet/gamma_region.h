#pragma once

#include <core/fitting/data_model/data_model.h>

#include <core/gamma/hypermet/peak.h>
#include <core/gamma/hypermet/poly_background.h>

#include <set>


namespace DAQuiri
{

//enum class Type : uint16_t
//{
//  Normal = 0,
//  Annihilation = 1,
//  Boron = 2,
//  GeTriangle = 3,
//  IntPeakShape = 4
//};
//Type region_type{Type::Normal};


class Region : public FittableRegion
{
  std::uniform_real_distribution<double> x_dist {-M_PI_2, M_PI_2};

 public:
  PolyBackground background;
  Peak default_peak_;
  std::map<double, Peak> peaks_;
  SUM4Edge LB_, RB_;

 public:
  Region();
  Region(const WeightedData& data, uint16_t background_samples);
  void replace_data(const WeightedData& data, const SUM4Edge& lb, const SUM4Edge& rb);
  void replace_data(const WeightedData& data, uint16_t left_samples, uint16_t right_samples);
  void replace_data(const WeightedData& data);
  void adjust_LB(const SUM4Edge& lb);
  void adjust_RB(const SUM4Edge& rb);

  double left() const;
  double right() const;
  bool empty() const;
  bool dirty() const;

  bool add_peak(double l, double r, double amp_hint = 0.0);
  bool adjust_sum4(double peakID, double left, double right);
  void auto_sum4();
  bool replace_hypermet(double peakID, Peak hyp);
  bool remove_peak(double peakID);
  bool remove_peaks(const std::set<double>& ids);

  // Fittable implementation
  void update_indices() override;
  Eigen::VectorXd variables() const override;
  double eval(double chan) const override;
  double eval_at(double chan, const Eigen::VectorXd& fit) const override;
  double eval_grad_at(double chan, const Eigen::VectorXd& fit, Eigen::VectorXd& grads) const override;
  void save_fit(const FitResult& result) override;

  bool perturb(std::mt19937& rng) override;
  bool sane() const override;

  void save_fit(const Eigen::VectorXd& variables);
  void reindex_peaks();

  std::string to_string(std::string prepend = "") const override;
  friend void to_json(nlohmann::json& j, const Region& s);
  friend void from_json(const nlohmann::json& j, Region& s);

 private:
  bool dirty_{false};

  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak

  void init_background();
  void cull_peaks();
};

}
