#pragma once

#include <core/fitting/hypermet/Hypermet.h>
#include <core/fitting/hypermet/PolyBackground.h>
#include <core/fitting/finder.h>

#include <core/fitting/BFGS/Fittable.h>

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


class Region : public Fittable
{
 public:
  PolyBackground background;
  Hypermet default_peak_;
  std::vector<Hypermet> peaks_;

 public:
  Region(const SpectrumData& data);
  double chi_sq_normalized() const;

  // Fitting related
  void map_fit();
  void save_fit_uncerts(const FitResult& result);

  // Fittable implementation
  std::vector<double> variables() const override;
  double degrees_of_freedom() const  override;
  double chi_sq(const std::vector<double>& fit) const override;
  double grad_chi_sq(const std::vector<double>& fit,
                     std::vector<double>& gradients) const override;

 private:
  int32_t var_count_{0};
  SpectrumData data_;

  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak

  size_t fit_var_count() const;
  void save_fit(const std::vector<double>& variables);
  double chi_sq() const;
  double grad_chi_sq(std::vector<double>& gradients) const;
};

}
