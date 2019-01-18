#pragma once

#include <cstdint>
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

class EfficiencyCal
{
  // This class loads and calculates the efficiency from an existing fit, as done by Hypermet-PC.
  // See also in:
  // G.L. Molnar, Zs. Revay, T. Belgya:
  // Wide energy range efficiency calibration method for Ge detectors,
  // Nucl. Instr. Meth. A 489 (2002) 140?159
 private:
  float e_maxdeg;
  double e_apol[7], e_bpol[6], e_normfact[8], e_poly_coeff[8];

  Eigen::SparseMatrix<double> e_VarMatrix{8, 8};

  double e_c0, e_c1;
  bool e_init_done{false};

 public:
  void Init(std::string flnm);
  void Close();
  bool InitDone() const;
  double Value(double Energy) const;
  double SigmaRel(double ByVal) const;

 private:
  double e_ortpol(size_t n, double X) const;
};
