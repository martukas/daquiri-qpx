#pragma once

#include <optimizerBFGS/Region.h>
#include <cstdint>
#include <vector>

namespace Hypermet
{

class BFGS
{
 public:
  bool Cancelfit{false};

  int32_t nfuncprofile;
  int32_t DiffType{1};

  void BFGSMin(Region& RegObj, double tolf, size_t& iter);

 private:
  double Sign(double a, double b);
  double BrentDeriv(const Region& region, double a, double b, double c, double tol, double& xmin,
                    const std::vector<double>& gx, const std::vector<double>& gh);
  void Bracket(const Region& region, double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const std::vector<double>& gx, const std::vector<double>& gh);
  double fgv(const Region& region, double lambda, std::vector<double> gx, std::vector<double> g);
  double dfgv(const Region& region, double lambda, std::vector<double> gx, std::vector<double> g);
  void LinMin(const Region& region, std::vector<double>& x, std::vector<double> h, double& fmin);
};

}
