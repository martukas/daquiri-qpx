#include <optimizerBFGS/SpectrumIO.h>
#include <fstream>

#include <core/util/custom_logger.h>

void CEff::Init(std::string flnm)
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

void CEff::Close()
{
  e_init_done = false;
}

bool CEff::InitDone() const
{
  return e_init_done;
}

double CEff::e_ortpol(size_t n, double X) const
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

double CEff::Value(double Energy) const
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

double CEff::SigmaRel(double Energy)
{
  try
  {
    double w;
    int32_t j;

    if (!e_init_done)
      return 0;

    double X{e_c0 + std::log(Energy) * e_c1};
    double sigeff{0};
    for (size_t i = 0; i <= e_maxdeg; ++i)
    {
      w = e_normfact[i] * e_ortpol(i, X);
      sigeff += w * w * e_VarMatrix.coeff(i, i);
      for (size_t j = 0; j < i; ++j)
        sigeff += 2 * w * e_normfact[j] * e_ortpol(j, X) * e_VarMatrix.coeff(i, j);
    }
    if (sigeff >= 0)
      sigeff = std::sqrt(sigeff);
    else
      sigeff = 0;
    return sigeff;
  }
  catch (...)
  {
    ERR("Efficiency error: {exception}");
    return 0;
  }
}

void CNonlin::Init(std::string flnm)
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

void CNonlin::Close()
{
  n_init_done = false;
}

bool CNonlin::InitDone() const
{
  return n_init_done;
}

double CNonlin::n_ortpol(size_t n, double X) const
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

double CNonlin::Value(double Position) const
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

double CNonlin::Sigma(double Position)
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
      siglin += w * w * n_VarMatrix.coeff(i, i);
      for (size_t j = 0; j <= i - 1; ++i)
        siglin += 2 * w * n_normfact[j + 2] * n_ortpol(j + 2, X) * n_VarMatrix.coeff(i, j);
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

double CNonlin::nonlin1(double Position)
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

void CNonlin::SetBasePoints(double ch1, double ch2)
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

uint8_t CCalibration::CalOrder() const
{
  return _CalOrder;
}
void CCalibration::CalOrder(uint8_t Value)
{
  if (Value <= 2)
  {
    _CalOrder = Value;
    EnergyCal.resize(Value);
    WidthCal.resize(Value);
  }
}

void CCalibration::New()
{
  //EnergyCal(0).Channel = 0
  //EnergyCal(0).Value = 0
  //EnergyCal(1).Channel = 1
  //EnergyCal(1).Value = 1
}

double CCalibration::EnergyConst() const
{
  double c0 = -1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[1].Value + 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[0].Value + EnergyCal[0].Value;
  return c0;
}

double CCalibration::EnergySlope() const
{
  double c1 = (1 / (EnergyCal[1].Channel - EnergyCal[0].Channel)
      * EnergyCal[1].Value - 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) * EnergyCal[0].Value);
  return c1;
}

double CCalibration::ChannelToEnergy(double Channel) const
{
  Channel += Nonlinearity.Value(Channel);
  double c0 = -1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[1].Value + 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[0].Value + EnergyCal[0].Value;
  //(1 / (Ch2 - Ch1) * e2 - 1 / (Ch2 - Ch1) * e1)
  double c1 = (1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) * (EnergyCal[1].Value - EnergyCal[0].Value));
  return c0 + c1 * Channel;
}

double CCalibration::EnergyToChannel(double Energy) const
{
  double c0 = -1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Channel * EnergyCal[1].Value + 1 /
      (EnergyCal[1].Channel - EnergyCal[0].Channel) * EnergyCal[0].Channel *
      EnergyCal[0].Value + EnergyCal[0].Value;
  double
      c1 = (1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[1].Value - 1 / (EnergyCal[1].Channel - EnergyCal[0].Channel) *
      EnergyCal[0].Value);
  Energy -= c0;
  Energy /= c1;
  Energy -= Nonlinearity.Value(Energy);
  return Energy;
}

double CCalibration::Width(double Channel) const
{
  //fwhm_b = (-fwhm2 ^ 2 + fwhm1 ^ 2) / (-Energy(chfwhm2) + Energy(chfwhm1))
  //fwhm_a = -(Energy(chfwhm2) * fwhm1 ^ 2 - fwhm2 ^ 2 * Energy(chfwhm1)) / (-Energy(chfwhm2) + Energy(chfwhm1))
}

double CSpectrum::Weight(size_t i) const
{
  double k0 = Channel[i];

  if (k0 >= 25)
    return std::sqrt(k0);
  else
  {
    double k1 = 1;
    if ((i > 0) && (i < Channel.size()))
      k1 = Channel[i - 1] + Channel[i] + Channel[i + 1] / 3.0;
    return std::max(std::sqrt(k1), 1.0);
  }
}

// template type for Val;
int8_t CSpectrum::Sign(double Val)
{
  if (Val < 0)
    return -1;
  if (Val > 0)
    return 1;
  return 0;
}

double CSpectrum::DeadTime(double TrueTime, double LiveTime)
{
  if (TrueTime > 0.0)
    return (TrueTime - LiveTime) / TrueTime * 100.0;
  return 0.0;
}

double CSpectrum::Rate(double LiveTime, double SumCounts)
{
  if (LiveTime > 0.0)
    return SumCounts / LiveTime * 100;
  return 0.0;
}

size_t CSpectrum::mystery_function(double Val)
{
  bool Ready = false;
  if ((Val >= 0) && (Val < pow(2, 36)))
    Ready = true;

  while (!Ready)
  {
    size_t exponent = std::log(std::abs(Val)) / std::log(2);
    Val = Val - Sign(Val) * pow(2, exponent);
    if ((Val >= 0) && (Val < pow(2, 36)))
      Ready = true;
  }
  return Val;
}

