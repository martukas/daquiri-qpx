#pragma once

#include <core/fitting/optimizer.h>

#include "TF1.h"
#include "TH1.h"
#include "TGraphErrors.h"

namespace DAQuiri
{

class WrapperROOT
{
 public:
  //General:
  static std::string def(uint16_t num, const Parameter&);
  static void set(TF1* f, uint16_t num, const Parameter&);
  // \todo use undertain type
  static std::pair<double, double> get(TF1* f, uint16_t num);

  static TH1D* fill_hist(const std::vector<double>& x,
                         const std::vector<double>& y);

  //CoefFunction:
  static uint16_t set_params(const CoefFunction&, TF1* f, uint16_t start);
  static uint16_t get_params(CoefFunction&, TF1* f, uint16_t start);
  static std::string def_of(CoefFunction* t);
  static std::string define_chain(const CoefFunction& func, const std::string& element);
  static std::string def_Polynomial(const CoefFunction& p);
  static std::string def_PolyLog(const CoefFunction& p);
  static std::string def_SqrtPoly(const CoefFunction& p);
  static std::string def_LogInverse(const CoefFunction& p);
  static std::string def_Effit(const CoefFunction& p);

  //Gaussian:
  static std::string def(const Gaussian& gaussian, uint16_t start = 0);
  static std::string def(const Gaussian& gaussian, uint16_t a, uint16_t c, uint16_t w);

  static void set_params(const Gaussian& gaussian, TF1* f, uint16_t start = 0);
  static void set_params(const Gaussian& gaussian, TF1* f, uint16_t a, uint16_t c, uint16_t w);
  static void get_params(Gaussian& gaussian, TF1* f, uint16_t start = 0);
  static void get_params(Gaussian& gaussian, TF1* f, uint16_t a, uint16_t c, uint16_t w);

  static std::vector<Gaussian> fit_multi_variw(const std::vector<double>& x,
                                               const std::vector<double>& y,
                                               std::vector<Gaussian> old,
                                               Polynomial& background,
                                               FitSettings settings
  );

  static std::vector<Gaussian> fit_multi_commonw(const std::vector<double>& x,
                                                 const std::vector<double>& y,
                                                 std::vector<Gaussian> old,
                                                 Polynomial& background,
                                                 FitSettings settings
  );

  //Hypermet:
  static std::string def(const Hypermet& h, uint16_t start = 0);
  static std::string def(const Hypermet& h, uint16_t width, uint16_t others_start);
  static std::string def_skew(const Parameter& ampl, const Parameter& slope,
                              const std::string& w, const std::string& xc,
                              const std::string& xcw, uint16_t idx);

  static void set_params(const Hypermet& hyp, TF1* f, uint16_t start = 0);
  static void set_params(const Hypermet& hyp, TF1* f, uint16_t w, uint16_t others_start);
  static void get_params(Hypermet& hyp, TF1* f, uint16_t start = 0);
  static void get_params(Hypermet& hyp, TF1* f, uint16_t w, uint16_t start);

  static std::vector<Hypermet> fit_multi_variw(const std::vector<double>& x,
                                               const std::vector<double>& y,
                                               std::vector<Hypermet> old,
                                               Polynomial& background,
                                               FitSettings settings
  );

  static std::vector<Hypermet> fit_multi_commonw(const std::vector<double>& x,
                                                 const std::vector<double>& y,
                                                 std::vector<Hypermet> old,
                                                 Polynomial& background,
                                                 FitSettings settings
  );

};

}
