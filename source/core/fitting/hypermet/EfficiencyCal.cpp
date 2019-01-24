#include <core/fitting/hypermet/EfficiencyCal.h>
#include <core/util/more_math.h>

#include <fstream>

#include <core/util/custom_logger.h>

namespace Hypermet
{

void EfficiencyCal::load(std::string flnm)
{
  try
  {
    double junk;
    bool flag{true};

    std::ifstream file(flnm, std::ios::binary);

    float e_maxdeg;
    file.read(reinterpret_cast<char*>(&e_maxdeg), sizeof(e_maxdeg));

//            junk = (e_maxdeg And Not 15) >> 4
//            e_maxdeg = e_maxdeg And 15
//            If (e_maxdeg > 8 Or e_maxdeg < 1 Or junk > 1) Then
//                'illegal dimensions
//                Throw New Exception("Illegal dimensions!")
//            End If
//            FileGet(1, e_c0)
//            FileGet(1, e_c1)
//            ReDim e_apol(e_maxdeg - 1)
//            ReDim e_bpol(e_maxdeg - 2)
//            ReDim e_normfact(e_maxdeg)
//            ReDim e_poly_coeff(e_maxdeg)
//            FileGet(1, e_apol)
//            FileGet(1, e_bpol)
//            FileGet(1, e_normfact)
//
//            FileGet(1, e_poly_coeff, , False)
    std::vector<double> junk2(e_maxdeg);

    variance_.resize(e_maxdeg, e_maxdeg);
    for (size_t i = 0; i <= e_maxdeg; ++i)
    {
      //FileGet(1, junk2, , False)
      for (size_t k = 0; k <= e_maxdeg; ++k)
        variance_.coeffRef(i, k) = junk2[k];
    }

    initialized_ = true;

  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    initialized_ = false;
    std::throw_with_nested(std::runtime_error("Efficiency error"));
  }
}

void EfficiencyCal::close()
{
  initialized_ = false;
}

bool EfficiencyCal::initialized() const
{
  return initialized_;
}

double EfficiencyCal::e_ortpol(size_t n, double X) const
{
  try
  {
    std::vector<double> pol(n + 1);
    pol[0] = 1;
    if (n != 0)
      pol[1] = (X - apol_[0]);

    for (size_t i = 2; i <= n; ++i)
      pol[i] = (X - apol_[i - 1]) * pol[i - 1] - bpol_[i - 2] * pol[i - 2];
    return pol[n];
  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    return 0;
  }
}

double EfficiencyCal::val(double energy) const
{
  try
  {
    if (!initialized_)
      return 1;

    double retval{0};
    double X{e_c0 + std::log(energy) * e_c1};
    for (size_t i = 0; i < coefficients_.size(); ++i)
      retval += coefficients_[i] * norm_factors_[i] * e_ortpol(i, X);
    return std::exp(retval);
  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    return 0;
  }
}

double EfficiencyCal::sigma_rel(double val) const
{
  try
  {
    if (!initialized_)
      return 0;

    double X{e_c0 + std::log(val) * e_c1};
    double sigeff{0};
    for (size_t i = 0; i < coefficients_.size(); ++i)
    {
      double w{norm_factors_[i] * e_ortpol(i, X)};
      sigeff += square(w) * variance_.coeff(i, i);
      for (size_t j = 0; j < i; ++j)
        sigeff += 2 * w * norm_factors_[j] * e_ortpol(j, X) * variance_.coeff(i, j);
    }
    if (sigeff >= 0.0)
      return std::sqrt(sigeff);
    else
      return 0.0;
  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    return 0;
  }
}

}
