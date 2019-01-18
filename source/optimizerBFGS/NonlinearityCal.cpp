#include <optimizerBFGS/NonlinearityCal.h>
#include <optimizerBFGS/more_math.h>

#include <fstream>

#include <core/util/custom_logger.h>

void NonlinearityCal::Init(std::string flnm)
{
  try
  {
    double junk;
    bool flag{true};

    std::ifstream file(flnm, std::ios::binary);

    file.read(reinterpret_cast<char*>(&n_maxdeg), sizeof(n_maxdeg));

//            FileGet(1, n_maxdeg)
//            FileGet(1, n_c0)
//            FileGet(1, n_c1)
//            ReDim n_apol(n_maxdeg - 1)
//            ReDim n_bpol(n_maxdeg - 2)
//            ReDim n_normfact(n_maxdeg)
//            ReDim n_poly_coeff(n_maxdeg - 2)
//            FileGet(1, n_apol)
//            FileGet(1, n_bpol)
//            FileGet(1, n_normfact)
//            FileGet(1, junk)
//            FileGet(1, junk)
//            FileGet(1, n_poly_coeff, , False)
//            Dim junk2(n_maxdeg + 1) As Double
//            FileGet(1, junk2, , False)
//            FileGet(1, junk2, , False)
    std::vector<double> junk2(n_maxdeg);

    n_VarMatrix.resize(n_maxdeg - 2, n_maxdeg - 2);
    for (size_t i = 0; i <= (n_maxdeg - 2); ++i)
    {
      //FileGet(1, junk2, , False)
      for (size_t k = 0; k <= (n_maxdeg - 2); ++k)
        n_VarMatrix.coeffRef(i, k) = junk2[k];
    }

    n_init_done = true;
  }
  catch (...)
  {
    ERR("Nonlinearity error: {exception}");
    n_init_done = false;
  }


//        try
//        {
//            SetBasePoints(Spectrum.Calibration.EnergyCal(0).Channel,
//                Spectrum.Calibration.EnergyCal(1).Channel);
//        }
//        catch (...)
//        {
//            ERR("Nonlinearity setbasepoint error: {exception}");
//            n_bl0 = 0;
//            n_bl1 = 0;
//        }
}

void NonlinearityCal::Close()
{
  n_init_done = false;
}

bool NonlinearityCal::InitDone() const
{
  return n_init_done;
}

double NonlinearityCal::n_ortpol(size_t n, double X) const
{
  try
  {
    std::vector<double> pol(n + 1);
    pol[0] = 1;
    pol[1] = (X - n_apol[0]);

    for (size_t i = 2; i <= n; ++i)
      pol[i] = (X - n_apol[i - 1]) * pol[i - 1] - n_bpol[i - 2] * pol[i - 2];
    pol[n];
  }
  catch (...)
  {
    ERR("Nonlinearity error: {exception}");
    return 0;
  }
}

double NonlinearityCal::Value(double Position) const
{
  try
  {
    if (!n_init_done)
      return 0;
    double X = n_c0 + (Position + 1) * n_c1;
    double retval = n_bl0 + n_bl1 * X;
    for (size_t i = 0; i <= n_maxdeg - 2; ++i)
      retval += n_poly_coeff[i] * n_normfact[i + 2] * n_ortpol(i + 2, X);
    return retval;
  }
  catch (...)
  {
    ERR("Nonlinearity error: {exception}");
    return 0;
  }
}

double NonlinearityCal::Sigma(double Position)
{
  try
  {
    if (!n_init_done)
      return 0;
    double X = n_c0 + (Position + 1) * n_c1;
    double siglin = 0;
    for (size_t i = 0; i <= n_maxdeg - 2; ++i)
    {
      double w = n_normfact[i + 2] * n_ortpol(i + 2, X);
      siglin += square(w) * n_VarMatrix.coeff(i, i);
      for (size_t j = 0; j <= i - 1; ++i)
        siglin += 2 * w * n_normfact[j + 2] *
            n_ortpol(j + 2, X) *
            n_VarMatrix.coeff(i, j);
    }
    if (siglin >= 0)
      siglin = std::sqrt(siglin);
    else
      siglin = 0;
    return siglin;
  }
  catch (...)
  {
    ERR("Nonlinearity error: {exception}");
    return 0;
  }
}

double NonlinearityCal::nonlin1(double Position)
{
  try
  {
    if (!n_init_done)
      return 0;
    double X = n_c0 + (Position + 1) * n_c1;
    double retval = 0;
    for (size_t i = 0; i <= n_maxdeg - 2; ++i)
      retval += n_poly_coeff[i] * n_normfact[i + 2] * n_ortpol(i + 2, X);
    return retval;
  }
  catch (...)
  {
    ERR("Nonlinearity error: {exception}");
    return 0;
  }
}

void NonlinearityCal::SetBasePoints(double ch1, double ch2)
{
  if (!n_init_done)
    throw (std::runtime_error("No nonlinearity correction loaded"));

  try
  {
    double m1 = nonlin1(ch1);
    double m2 = nonlin1(ch2);
    n_bl1 = (m2 - m1) / (n_c1 * (ch1 - ch2));
    n_bl0 = -m1 - n_bl1 * (n_c0 + n_c1 * ch1);
  }
  catch (...)
  {
    ERR("Nonlinearity error: {exception}");
    return;
  }
}
