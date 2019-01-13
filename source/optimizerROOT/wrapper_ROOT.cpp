#include <optimizerROOT/wrapper_ROOT.h>

//#include <core/util/custom_logger.h>

namespace DAQuiri
{

std::string WrapperROOT::def_of(CoefFunction* t)
{
  if (t->type() == "Polynomial")
    return def_Polynomial(*t);
  else if (t->type() == "PolyLog")
    return def_PolyLog(*t);
  else if (t->type() == "LogInverse")
    return def_LogInverse(*t);
  else if (t->type() == "SqrtPoly")
    return def_SqrtPoly(*t);
  else if (t->type() == "Effit")
    return def_Effit(*t);
  return "";
}

std::string WrapperROOT::def_Polynomial(const CoefFunction& p)
{
  std::string xc = "x";
  if (std::isfinite(p.x_offset().value()) && p.x_offset().value())
    xc = "(x-" + std::to_string(p.x_offset().value()) + ")";

  return define_chain(p, xc);
}

std::string WrapperROOT::def_PolyLog(const CoefFunction& p)
{
  std::string xc = "TMath::Log(x)";
  if (std::isfinite(p.x_offset().value()) && p.x_offset().value())
    xc = "TMath::Log(x-" + std::to_string(p.x_offset().value()) + ")";

  return "TMath::Exp(" + define_chain(p, xc) + ")";
}

std::string WrapperROOT::def_SqrtPoly(const CoefFunction& p)
{
  std::string xc = "x";
  if (std::isfinite(p.x_offset().value()) && p.x_offset().value())
    xc = "(x-" + std::to_string(p.x_offset().value()) + ")";

  return "TMath::Sqrt(" + define_chain(p, xc) + ")";
}

std::string WrapperROOT::def_LogInverse(const CoefFunction& p)
{
  std::string xc = "(1.0/x)";
  if (std::isfinite(p.x_offset().value()) && p.x_offset().value())
    xc = "1.0/(x-" + std::to_string(p.x_offset().value()) + ")";

  return "TMath::Exp(" + define_chain(p, xc) + ")";
}

std::string WrapperROOT::def_Effit(const CoefFunction& p)
{
  std::string definition = "((d + e*xb + f*(xb^2))^(-20))^(-0.05) where xb=ln(x/1000)";
  //   x= ((A + B*ln(x/100.0) + C*(ln(x/100.0)^2))^(-G) + (D + E*ln(x/1000.0) + F*(ln(x/1000.0)^2))^(1-G))^(1-(1/G))
  //pow(pow(A + B*xa + C*xa*xa,-G) + pow(D + E*xb + F*xb*xb,-G), -1.0/G);
//  if (G==0)
//      G = 20;
  return definition;
}

std::string WrapperROOT::def(uint16_t num, const Parameter& param)
{
  if (!param.fixed())
    return "[" + std::to_string(num) + "]";
  else if (std::isfinite(param.value()))
    return std::to_string(param.value());
  else
    return "0";
}

void WrapperROOT::set(TF1* f, uint16_t num, const Parameter& param)
{
  f->SetParameter(num, param.value());
  f->SetParLimits(num, param.lower(), param.upper());
}

std::pair<double, double> WrapperROOT::get(TF1* f, uint16_t num)
{
  return {f->GetParameter(num), f->GetParError(num)};
}

uint16_t WrapperROOT::set_params(const CoefFunction& func, TF1* f, uint16_t start)
{
  for (const auto& p : func.coeffs())
    if (!p.second.fixed())
    {
      set(f, start, p.second);
      start++;
    }
  return start;
}

uint16_t WrapperROOT::get_params(CoefFunction& func, TF1* f, uint16_t start)
{
  for (auto p : func.coeffs())
    if (!p.second.fixed())
    {
      // \todo use uncertainty
      func.set_coeff(p.first, get(f, start).first);
      start++;
    }
  return start;
}

std::string WrapperROOT::define_chain(const CoefFunction& func,
                                        const std::string& element)
{
  std::string definition;
  int i = 0;
  for (auto& c : func.coeffs())
  {
    if (i > 0)
      definition += "+";
    if (!c.second.fixed())
    {
      definition += def(i, c.second);
      i++;
    }
    else
      definition += std::to_string(c.second.value());
    if (c.first > 0)
      definition += "*" + element;
    if (c.first > 1)
      definition += "^" + std::to_string(c.first);
  }
  return definition;
}

std::string WrapperROOT::def(const Gaussian& gaussian, uint16_t start)
{
  return def(gaussian, start, start + 1, start + 2);
}

std::string WrapperROOT::def(const Gaussian& gaussian, uint16_t a, uint16_t c, uint16_t w)
{
  return def(a, gaussian.height()) +
      "*TMath::Exp(-((x-" + def(c, gaussian.center()) + ")/" + def(w, gaussian.hwhm()) + ")^2)";
}

void WrapperROOT::set_params(const Gaussian& gaussian, TF1* f, uint16_t start)
{
  set_params(gaussian, f, start, start + 1, start + 2);
}

void WrapperROOT::set_params(const Gaussian& gaussian, TF1* f, uint16_t a, uint16_t c, uint16_t w)
{
  set(f, a, gaussian.height());
  set(f, c, gaussian.center());
  set(f, w, gaussian.hwhm());
}

void WrapperROOT::get_params(Gaussian& gaussian, TF1* f, uint16_t start)
{
  get_params(gaussian, f, start, start + 1, start + 2);
}

void WrapperROOT::get_params(Gaussian& gaussian, TF1* f, uint16_t a, uint16_t c, uint16_t w)
{
  // \todo use uncertainties
  gaussian.set_height(get(f, a).first);
  gaussian.set_center(get(f, c).first);
  gaussian.set_hwhm(get(f, w).first);
}

TH1D* WrapperROOT::fill_hist(const std::vector<double>& x,
                               const std::vector<double>& y)
{
  if ((x.size() < 1) || (x.size() != y.size()))
    return nullptr;

  TH1D* h1 = new TH1D("h1", "h1", x.size(), x.front(), x.back());
  int i = 1;
  for (const auto& h : y)
  {
    h1->SetBinContent(i, h);
    h1->SetBinError(i, sqrt(h));
    i++;
  }
  return h1;
}

std::vector<Gaussian> WrapperROOT::fit_multi_variw(const std::vector<double>& x,
                                                     const std::vector<double>& y,
                                                     std::vector<Gaussian> old,
                                                     Polynomial& background,
                                                     FitSettings settings)
{
  if (old.empty())
    return old;

  TH1D* h1 = fill_hist(x, y);
  if (!h1)
    return old;

  for (auto& gaussian : old)
  {
    Optimizer::sanity_check(gaussian, x.front(), x.back(), h1->GetMinimum(), h1->GetMaximum());
    Optimizer::constrain_center(gaussian, settings.lateral_slack);

    double width_expected = gaussian.hwhm().value();
    if (settings.cali_fwhm_.valid() && settings.cali_nrg_.valid())
      width_expected = settings.bin_to_width(gaussian.center().value()) / (2 * sqrt(log(2)));

    gaussian.constrain_hwhm(width_expected * settings.width_common_bounds.lower(),
                            width_expected * settings.width_common_bounds.upper());
  }

  std::string definition = def_of(&background);
  size_t backgroundparams = background.coeffs().size(); //valid coef count?
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = backgroundparams + i * 3;
    definition += "+" + def(old.at(i), num);
  }

  TF1* f1 = new TF1("f1", definition.c_str());

  set_params(background, f1, 0);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = backgroundparams + i * 3;
    set_params(old[i], f1, num);
  }

  h1->Fit("f1", "QEMN");

  get_params(background, f1, 0);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = backgroundparams + i * 3;
    get_params(old[i], f1, num);
    old[i].set_chi2(f1->GetChisquare());
  }
  background.chi2(f1->GetChisquare());

  f1->Delete();
  h1->Delete();

  return old;
}

std::vector<Gaussian> WrapperROOT::fit_multi_commonw(const std::vector<double>& x,
                                                       const std::vector<double>& y,
                                                       std::vector<Gaussian> old,
                                                       Polynomial& background,
                                                       FitSettings settings)
{
  if (old.empty())
    return old;

  TH1D* h1 = fill_hist(x, y);
  if (!h1)
    return old;

  Parameter w_common{0.0};
  double width_expected = settings.bin_to_width((x.front() + x.back()) / 2) / (2 * sqrt(log(2)));
  w_common.constrain(width_expected * settings.width_common_bounds.lower(),
                     width_expected * settings.width_common_bounds.upper());

  for (auto& gaussian : old)
  {
    gaussian.set_hwhm(w_common);
    Optimizer::sanity_check(gaussian, x.front(), x.back(), h1->GetMinimum(), h1->GetMaximum());
    Optimizer::constrain_center(gaussian, settings.lateral_slack);
  }

  std::string definition = def_of(&background);
  size_t w_index = background.coeffs().size(); //valid coef count?
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = 1 + w_index + i * 2;
    definition += "+" + def(old.at(i), num, num + 1, w_index);
  }

  TF1* f1 = new TF1("f1", definition.c_str());

  set_params(background, f1, 0);
  set(f1, w_index, w_common);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = 1 + w_index + i * 2;
    set_params(old[i], f1, num, num + 1, w_index);
  }

  h1->Fit("f1", "QEMN");

  get_params(background, f1, 0);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = 1 + w_index + i * 2;
    get_params(old[i], f1, num, num + 1, w_index);
    old[i].set_chi2(f1->GetChisquare());
  }
  background.chi2(f1->GetChisquare());

  f1->Delete();
  h1->Delete();

  return old;
}

std::string WrapperROOT::def(const Hypermet& h, uint16_t start)
{
  return def(h, start, start + 1);
}

std::string WrapperROOT::def_skew(const Parameter& ampl, const Parameter& slope,
                                    const std::string& w, const std::string& xc,
                                    const std::string& xcw, uint16_t idx)
{
  // \todo make use of this function
//  if (!ampl.enabled())
//    return "";
  std::string skewh = def(idx, ampl);
  std::string skews = "/" + def(idx + 1, slope);
  return skewh + "*TMath::Exp((0.5*" + w + skews + ")^2 + (" + xc + skews + "))*TMath::Erfc((0.5*" + w + skews + ") + "
      + xcw + ")";
}

std::string WrapperROOT::def(const Hypermet& hyp, uint16_t width, uint16_t i)
{
  std::string h = def(i, hyp.height());
  std::string w = def(width, hyp.width());
  std::string xc = "(x-" + def(i + 1, hyp.center()) + ")";
  std::string xcw = xc + "/" + w;
  std::string lskew = "0";
  std::string rskew = "0";
  std::string tail = "0";
  std::string step = "0";

  std::string lskewh = def(i + 2, hyp.Lskew_amplitude());
  std::string lskews = "/" + def(i + 3, hyp.Lskew_slope());
  lskew =
      lskewh + "*TMath::Exp((0.5*" + w + lskews + ")^2 + (" + xc + lskews + "))*TMath::Erfc((0.5*" + w + lskews + ") + "
          + xcw + ")";

  std::string rskewh = def(i + 4, hyp.Rskew_amplitude());
  std::string rskews = "/" + def(i + 5, hyp.Rskew_slope());
  rskew =
      rskewh + "*TMath::Exp((0.5*" + w + rskews + ")^2 - (" + xc + rskews + "))*TMath::Erfc((0.5*" + w + rskews + ") - "
          + xcw + ")";

  std::string tailh = def(i + 6, hyp.tail_amplitude());
  std::string tails = "/" + def(i + 7, hyp.tail_slope());
  tail = tailh + "*TMath::Exp((0.5*" + w + tails + ")^2 + (" + xc + tails + "))*TMath::Erfc((0.5*" + w + tails + ") + "
      + xcw + ")";

  std::string steph = def(i + 8, hyp.step_amplitude());
  step = steph + "*TMath::Erfc(" + xc + "/" + w + ")";

  return h + "*( TMath::Exp(-(" + xcw + ")^2) + 0.5 * (" + lskew + " + " + rskew + " + " + tail + " + " + step + " ) )";
}

void WrapperROOT::set_params(const Hypermet& hyp, TF1* f, uint16_t start)
{
  set_params(hyp, f, start, start + 1);
}

void WrapperROOT::set_params(const Hypermet& hyp, TF1* f, uint16_t width, uint16_t others_start)
{
  set(f, width, hyp.width());
  set(f, others_start, hyp.height());
  set(f, others_start + 1, hyp.center());
  set(f, others_start + 2, hyp.Lskew_amplitude());
  set(f, others_start + 3, hyp.Lskew_slope());
  set(f, others_start + 4, hyp.Rskew_amplitude());
  set(f, others_start + 5, hyp.Rskew_slope());
  set(f, others_start + 6, hyp.tail_amplitude());
  set(f, others_start + 7, hyp.tail_slope());
  set(f, others_start + 8, hyp.step_amplitude());
}

void WrapperROOT::get_params(Hypermet& hyp, TF1* f, uint16_t start)
{
  get_params(hyp, f, start, start + 1);
}

void WrapperROOT::get_params(Hypermet& hyp, TF1* f, uint16_t w, uint16_t others_start)
{
  // \todo use uncertainties
  hyp.set_width(get(f, w).first);
  hyp.set_height(get(f, others_start).first);
  hyp.set_center(get(f, others_start + 1).first);
  hyp.set_Lskew_amplitude(get(f, others_start + 2).first);
  hyp.set_Lskew_slope(get(f, others_start + 3).first);
  hyp.set_Rskew_amplitude(get(f, others_start + 4).first);
  hyp.set_Rskew_slope(get(f, others_start + 5).first);
  hyp.set_tail_amplitude(get(f, others_start + 6).first);
  hyp.set_tail_slope(get(f, others_start + 7).first);
  hyp.set_step_amplitude(get(f, others_start + 8).first);
}

std::vector<Hypermet> WrapperROOT::fit_multi_variw(const std::vector<double>& x,
                                                     const std::vector<double>& y,
                                                     std::vector<Hypermet> old,
                                                     Polynomial& background,
                                                     FitSettings settings)
{
  if (old.empty())
    return old;

  TH1D* h1 = fill_hist(x, y);
  if (!h1)
    return old;

  for (auto& hyp : old)
  {
    Optimizer::sanity_check(hyp, x.front(), x.back(), h1->GetMinimum(), h1->GetMaximum());
    Optimizer::constrain_center(hyp, settings.lateral_slack);

    double width_expected = hyp.width().value();
    if (settings.cali_fwhm_.valid() && settings.cali_nrg_.valid())
      width_expected = settings.bin_to_width(hyp.center().value()) / (2 * sqrt(log(2)));

    hyp.constrain_width(width_expected * settings.width_common_bounds.lower(),
                        width_expected * settings.width_common_bounds.upper());
  }

  std::string definition = def_of(&background);
  uint16_t backgroundparams = background.coeffs().size(); //valid coef count?
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = backgroundparams + i * 10;
    definition += "+" + def(old.at(i), num);
  }

  TF1* f1 = new TF1("f1", definition.c_str());

  set_params(background, f1, 0);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = backgroundparams + i * 10;
    set_params(old[i], f1, num);
  }

  h1->Fit("f1", "QEMN");

  get_params(background, f1, 0);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = backgroundparams + i * 10;
    get_params(old[i], f1, num);
    old[i].set_chi2(f1->GetChisquare());
  }
  background.chi2(f1->GetChisquare());

  f1->Delete();
  h1->Delete();

  return old;
}

std::vector<Hypermet> WrapperROOT::fit_multi_commonw(const std::vector<double>& x,
                                                       const std::vector<double>& y,
                                                       std::vector<Hypermet> old,
                                                       Polynomial& background,
                                                       FitSettings settings)
{
  if (old.empty())
    return old;

  TH1D* h1 = fill_hist(x, y);
  if (!h1)
    return old;

  Parameter w_common{0.0};
  double width_expected = settings.bin_to_width((x.front() + x.back()) / 2) / (2 * sqrt(log(2)));
  w_common.constrain(width_expected * settings.width_common_bounds.lower(),
                     width_expected * settings.width_common_bounds.upper());

  for (auto& hyp : old)
  {
    hyp.set_width(w_common);
    Optimizer::sanity_check(hyp, x.front(), x.back(), h1->GetMinimum(), h1->GetMaximum());
    Optimizer::constrain_center(hyp, settings.lateral_slack);
  }

  std::string definition = def_of(&background);
  uint16_t w_index = background.coeffs().size(); //valid coef count?
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = 1 + w_index + i * 9;
    definition += "+" + def(old.at(i), w_index, num);
  }

  TF1* f1 = new TF1("f1", definition.c_str());

  set_params(background, f1, 0);
  set(f1, w_index, w_common);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = 1 + w_index + i * 9;
    set_params(old[i], f1, w_index, num);
  }

  h1->Fit("f1", "QEMN");

  get_params(background, f1, 0);
  for (size_t i = 0; i < old.size(); ++i)
  {
    uint16_t num = 1 + w_index + i * 9;
    get_params(old[i], f1, w_index, num);
    old[i].set_chi2(f1->GetChisquare());
  }
  background.chi2(f1->GetChisquare());

  f1->Delete();
  h1->Delete();

  return old;
}

}
