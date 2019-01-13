#include <optimizerROOT/optimizer_ROOT.h>
#include <optimizerROOT/wrapper_ROOT.h>

//#include <core/util/custom_logger.h>

namespace DAQuiri
{

void OptimizerROOT::fit(std::shared_ptr<CoefFunction> func,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& x_sigma,
                        const std::vector<double>& y_sigma)
{
  if (!func || (x.size() != y.size()))
    return;

  TGraphErrors* g1;
  if ((x.size() == x_sigma.size()) && (x_sigma.size() == y_sigma.size()))
    g1 = new TGraphErrors(y.size(),
                          x.data(), y.data(),
                          x_sigma.data(), y_sigma.data());
  else
    g1 = new TGraphErrors(y.size(),
                          x.data(), y.data(),
                          nullptr, nullptr);

  std::string def = WrapperROOT::def_of(func.get());

  //DBG("Def = {}", def);

  TF1* f1 = new TF1("f1", def.c_str());
  WrapperROOT::set_params(*func, f1, 0);

  g1->Fit("f1", "QEMN");

  WrapperROOT::get_params(*func, f1, 0);
  func->chi2(f1->GetChisquare());

  f1->Delete();
  g1->Delete();
}

void OptimizerROOT::fit(Gaussian& gaussian,
                        const std::vector<double>& x,
                        const std::vector<double>& y)
{
  TH1D* h1 = WrapperROOT::fill_hist(x, y);
  if (!h1)
    return;

  initial_sanity(gaussian, x.front(), x.back(),
                 h1->GetMinimum(), h1->GetMaximum());

  TF1* f1 = new TF1("f1", WrapperROOT::def(gaussian).c_str());

  WrapperROOT::set_params(gaussian, f1);

  h1->Fit("f1", "QEMN");

  WrapperROOT::get_params(gaussian, f1);
  gaussian.set_chi2(f1->GetChisquare());

  f1->Delete();
  h1->Delete();
}

std::vector<Gaussian> OptimizerROOT::fit_multiplet(const std::vector<double>& x,
                                                   const std::vector<double>& y,
                                                   std::vector<Gaussian> old,
                                                   Polynomial& background,
                                                   FitSettings settings)
{
  bool use_w_common = (settings.width_common &&
      settings.cali_fwhm_.valid() &&
      settings.cali_nrg_.valid());

  if (use_w_common)
    return WrapperROOT::fit_multi_commonw(x, y, old, background, settings);
  else
    return WrapperROOT::fit_multi_variw(x, y, old, background, settings);
}

std::vector<Hypermet> OptimizerROOT::fit_multiplet(const std::vector<double>& x,
                                                   const std::vector<double>& y,
                                                   std::vector<Hypermet> old,
                                                   Polynomial& background,
                                                   FitSettings settings)
{
  bool use_w_common = (settings.width_common &&
      settings.cali_fwhm_.valid() &&
      settings.cali_nrg_.valid());

  if (use_w_common)
    return WrapperROOT::fit_multi_commonw(x, y, old, background, settings);
  else
    return WrapperROOT::fit_multi_variw(x, y, old, background, settings);
}

}
