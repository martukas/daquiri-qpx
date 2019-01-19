#pragma once

#include <cstdint>
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

namespace Hypermet
{

class NonlinearityCal
{
  // This class loads and calculates the nonlinearity from an existing fit, as done by Hypermet-PC.
  // See also in:
  // B. Fazekas, Zs. Revay, J. ?st?r, T. Belgya, G. Moln?r, A. Simonits:
  // A new method for determination of gamma-ray spectrometer non-linearity,
  // Nucl. Instr. Meth. A 422 (1999) 469-473
 private:
  float n_maxdeg;
  double n_apol[7], n_bpol[6], n_normfact[8], n_poly_coeff[7];
  Eigen::SparseMatrix<double> n_VarMatrix;

  double n_c0, n_c1;
  double n_bl0, n_bl1;
  bool n_init_done{false};

 public:
  void Init(std::string flnm);
  void Close();
  bool InitDone() const;
  double Value(double Position) const;
  double Sigma(double Position) const;
  void SetBasePoints(double ch1, double ch2);

 private:
  double n_ortpol(size_t n, double X) const;
  double nonlin1(double Position);
};

}
