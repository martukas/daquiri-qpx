#pragma once

#include <cstdint>
#include <vector>

struct Region {
  std::vector<double> Vector;
  double DegreeOfFreedom;

  static double CalcChiSq(const std::vector<double>& v);
  static void GradChiSq(const std::vector<double>& v, std::vector<double>& v2, double&);
};

class BFGS {
 public:
  bool Cancelfit {false};

  int32_t nfuncprofile;
  int32_t DiffType {1};

  double erfc(double x);

  void BFGSMin(Region& RegObj, double tolf, size_t& iter);

 private:
  double Sign(double a, double b);
  double BrentDeriv(double a, double b, double c, double tol, double& xmin,
                    const std::vector<double>& gx, const std::vector<double>& gh);
  void Bracket(double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const std::vector<double>& gx, const std::vector<double>& gh);
  double fgv(double lambda, std::vector<double> gx, std::vector<double> g);
  double dfgv(double lambda, std::vector<double> gx, std::vector<double> g);
  void LinMin(std::vector<double>& x, std::vector<double> h, double& fmin);

};