#include <optimizerBFGS/EfficiencyCal.h>
#include <optimizerBFGS/more_math.h>

#include <fstream>

#include <core/util/custom_logger.h>

namespace Hypermet
{

void EfficiencyCal::Init(std::string flnm)
{
  try
  {
    double junk;
    bool flag{true};

    std::ifstream file(flnm, std::ios::binary);

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

    e_VarMatrix.resize(e_maxdeg, e_maxdeg);
    for (size_t i = 0; i <= e_maxdeg; ++i)
    {
      //FileGet(1, junk2, , False)
      for (size_t k = 0; k <= e_maxdeg; ++k)
        e_VarMatrix.coeffRef(i, k) = junk2[k];
    }

    e_init_done = true;

  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    e_init_done = false;
    std::throw_with_nested(std::runtime_error("Efficiency error"));
  }
}

void EfficiencyCal::Close()
{
  e_init_done = false;
}

bool EfficiencyCal::InitDone() const
{
  return e_init_done;
}

double EfficiencyCal::e_ortpol(size_t n, double X) const
{
  try
  {
    std::vector<double> pol(n + 1);
    pol[0] = 1;
    if (n != 0)
      pol[1] = (X - e_apol[0]);

    for (size_t i = 2; i <= n; ++i)
      pol[i] = (X - e_apol[i - 1]) * pol[i - 1] - e_bpol[i - 2] * pol[i - 2];
    return pol[n];
  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    return 0;
  }
}

double EfficiencyCal::val(double Energy) const
{
  try
  {
    if (!e_init_done)
      return 1;

    double retval{0};
    double X{e_c0 + std::log(Energy) * e_c1};
    for (size_t i = 0; i <= e_maxdeg; ++i)
      retval += e_poly_coeff[i] * e_normfact[i] * e_ortpol(i, X);
    return std::exp(retval);
  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    return 0;
  }
}

double EfficiencyCal::SigmaRel(double Energy) const
{
  try
  {
    if (!e_init_done)
      return 0;

    double X{e_c0 + std::log(Energy) * e_c1};
    double sigeff{0};
    for (size_t i = 0; i <= e_maxdeg; ++i)
    {
      double w{e_normfact[i] * e_ortpol(i, X)};
      sigeff += square(w) * e_VarMatrix.coeff(i, i);
      for (size_t j = 0; j < i; ++j)
        sigeff += 2 * w * e_normfact[j] * e_ortpol(j, X) * e_VarMatrix.coeff(i, j);
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
