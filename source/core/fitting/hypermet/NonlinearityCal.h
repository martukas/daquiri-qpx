#pragma once

#include <cstdint>
#include <string>

#include <core/util/eigen_fix.h>

namespace DAQuiri
{

class NonlinearityCal
{
  // This class loads and calculates the nonlinearity from an existing fit, as done by Hypermet-PC.
  // See also in:
  // B. Fazekas, Zs. Revay, J. ?st?r, T. Belgya, G. Moln?r, A. Simonits:
  // A new method for determination of gamma-ray spectrometer non-linearity,
  // Nucl. Instr. Meth. A 422 (1999) 469-473
 private:
  std::vector<double> apol_, bpol_, norm_factors_, coefficients_;
  Eigen::SparseMatrix<double> variance_;

  double n_c0, n_c1;
  double n_bl0, n_bl1;
  bool initialized_{false};

 public:
  void load(std::string flnm);
  void close();
  bool initialized() const;
  double val(double position) const;
  double sigma(double position) const;
  void set_base_points(double ch1, double ch2);

 private:
  double n_ortpol(size_t n, double X) const;
  double nonlin1(double Position);
};

}
