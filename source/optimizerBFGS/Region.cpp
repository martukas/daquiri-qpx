#include <optimizerBFGS/Region.h>
#include <optimizerBFGS/more_math.h>

#include <core/util/custom_logger.h>

namespace Hypermet
{

Region::Region(CSpectrum& spe, double FromChannel, double ToChannel)
    : spectrum(spe)
      , FirstChannel(std::min(FromChannel, ToChannel))
      , LastChannel(std::max(FromChannel, ToChannel))
{
  try
  {
    DEL.max(4);
    DEL.min(0.8);

    AST.max(1.5);
    AST.min(0.02);
    AST.to_fit = true;

    BST.max(0.5);
    BST.min(0.2);
    BST.to_fit = true;

    ALT.max(0.15);
    ALT.min(0.0001);
    ALT.to_fit = true;

    BLT.max(50);
    BLT.min(2.5);
    BLT.to_fit = true;

    ART.max(0.9);
    ART.min(0.01);
    ART.to_fit = true;

    BRT.max(1.5);
    BRT.min(0.3);
    BRT.to_fit = true;

    SIG.max(0.05);
    SIG.min(0.000001);
    SIG.to_fit = true;

    BSL.to_fit = true;
    BCV.to_fit = true;

    //  SearchPeaks();

  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("Object Region failed!"));
  }
}

bool Region::LeftTail() const
{
  //Get/Set Left Long Tail Status
  return TailFlag;
}

void Region::LeftTail(bool enable)
{
  TailFlag = enable;
  ALT.to_fit = enable;
  BLT.to_fit = enable;
}

bool Region::RightTail() const
{
  return RightTailFlag;
}
void Region::RightTail(bool enable)
{
  RightTailFlag = enable;
  ART.to_fit = enable;
  BRT.to_fit = enable;
}

bool Region::Slope() const
{
  return SlopeFlag;
}

void Region::Slope(bool enable)
{
  SlopeFlag = enable;
  BSL.to_fit = enable;
}

bool Region::Curve() const
{
  return CurveFlag;
}

void Region::Curve(bool enable)
{
  CurveFlag = enable;
  BCV.to_fit = enable;
}

bool Region::StepBkg() const
{
  //Get/Set Step Status
  return StepFlag;
}

void Region::StepBkg(bool enable)
{
  StepFlag = enable;
  SIG.to_fit = enable;
}

int32_t Region::L(int32_t i, int32_t j, int32_t m)
{
  //M = FWHM
  int32_t ret{0};

  if (j - m <= i && i <= j - 1)
    ret = -1;
  if (j <= i && i <= j + m - 1)
    ret = 2;
  if (j + m <= i && i <= j + 2 * m - 1)
    ret = -1;

  return ret;
}

void Region::SearchPeaks(uint8_t Threshold)
{
  int32_t m;
  try
  {
    m = static_cast<int32_t>(1.6551 * DEL.val());
  }
  catch (...)
  {
    m = 3;
  }

  size_t i;
  try
  {
    for (size_t j = FirstChannel; j <= LastChannel; ++j)
    {
      double val = 0;
      double Var = 0;
      for (i = j - m; i <= (j + 2 * m - 1); ++j)
        if (i > 1)
          val = val + L(i, j, m) * spectrum.Channel[i];
      Var += square(L(i, j, m)) * spectrum.Channel[i];

      //Conv(j - FirstChannel) = val / std::sqrt(Var)
      //if(((Conv(j - FirstChannel - 2) < Conv(j - FirstChannel - 1)) && _
      //(Conv(j - FirstChannel) < Conv(j - FirstChannel - 1)) && _
      //(Conv(j - FirstChannel - 1) > Threshold))) {
      //AddPeak(j - 1, j - 2, j, std::sqrt(spectrum.Channel[j]))
    }
  }
  catch (...)
  {
    ERR("Search Peaks failed!");
  }
}

void Region::AddPeak(double Position, double Min, double Max, double Gamma)
{
  try
  {
    peaks.emplace_back();

    peaks.back().GAM.x_index = peaks[std::max(peaks.size() - 2, size_t(0))].GAM.x_index + 2;
    peaks.back().position.x_index = peaks[std::max(peaks.size() - 2, size_t(0))].position.x_index + 2;

    peaks.back().position.min(Min);
    peaks.back().position.max(Max);
    peaks.back().position.val(Position);
    peaks.back().position.uncert_value = 0;

    peaks.back().GAM.val(Gamma);
    peaks.back().GAM.uncert_value = 0;

  }
  catch (...)
  {
    ERR("Add Peak failed!");
  }
}

void Region::DeletePeak(size_t index)
{
  try
  {
    if (index >= peaks.size())
    {
      ERR("Can't delete the peak! (invalid index)");
      return;
    }
    if (peaks.size() == 1)
    {
      ERR("Can't delete the only peak!");
      return;
    }
    for (size_t i = index; i < peaks.size() - 1; ++i)
    {
      peaks[i].position = peaks[i + 1].position;
      peaks[i].GAM = peaks[i + 1].GAM;
    } //i
    peaks.resize(peaks.size() - 1);
  }
  catch (...)
  {
    ERR("Delete Peak failed!");
  }
}

double Region::PeakArea(size_t PeakIndex) const
{
  return peaks[PeakIndex].GAM.val() * DEL.val() * (std::sqrt(M_PI) +
      AST.val() * BST.val() * std::exp(-0.25 / square(BST.val())) +
      ART.val() * BRT.val() * std::exp(-0.25 / square(BRT.val())));
}

double Region::UncPeakArea(size_t PeakIndex) const
{
  double t = PeakArea(PeakIndex);
  //, i, j As Integer
  //Dim cs As Double = ChisqNorm() * 0.5
  //for( i = 0 To FitVars - 1
  //for( j = 0 To i - 1
  //t += Gradient(i) * Gradient(j) * Hessinv.coeff(i, j) * cs
  //} //j
  //} //i
  //
  //Dim Bracket As Double = (std::sqrt(M_PI) + AST.val() * BST.val() * std::exp(-0.25 / BST.val() ^ 2) + ART.val() * BRT.val() * std::exp(-0.25 / BRT.val() ^ 2))
  //(dGAM/dX*dArea/dGAM)^2*Var(X)*Chisq
  //t = (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) ^ 2 * Hessinv.coeff(DEL.x_index, DEL.x_index) * cs
  //(dGAM/dX*dArea/dGAM)*(dDEL/dY*dArea/dDEL)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) * (DEL.GradAt(DEL.X) * Peak(PeakIndex).GAM.val() * Bracket) * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, DEL.x_index) * cs
  //if(AST.to_fit = true) {
  ////(dGAM/dX*dArea/dGAM)*(dAST/dY*dArea/dAST)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) * (AST.GradAt(AST.X) * Peak(PeakIndex).GAM.val() * DEL.val() * BST.val() * std::exp(-0.25 / BST.val() ^ 2)) * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, AST.x_index) * cs
  //}
  //if(BST.to_fit = true) {
  //(dGAM/dX*dArea/dGAM)*(dBST/dY*dArea/dBST)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.val() * Bracket) * (BST.GradAt(BST.X) * Peak(PeakIndex).GAM.val() * DEL.val() * AST.val() * (1 + 0.5 / BST.val() ^ 2) * std::exp(-0.25 / BST.val() ^ 2)) * Hessinv.coeff(Peak(PeakIndex).GAM.x_index, AST.x_index) * cs
  //}
  return std::sqrt(t) * std::max(1.0, ChisqNorm());
}

double Region::PeakAreaEff(size_t PeakIndex, const Calibration& cal)
{
  double
      eff = cal.efficiency.val(cal.channel_to_energy(peaks[PeakIndex].position.val()));
  //eff = 1 if uninitialized
  return (peaks[PeakIndex].GAM.val() * DEL.val() * (std::sqrt(M_PI) +
      AST.val() * BST.val() * std::exp(-0.25 / square(BST.val())) +
      ART.val() * BRT.val() * std::exp(-0.25 / square(BRT.val())))) / eff;
}

double Region::UncPeakAreaEff(size_t PeakIndex, const Calibration& cal)
{
  double t = PeakArea(PeakIndex);
  double
      eff = cal.efficiency.val(cal.channel_to_energy(peaks[PeakIndex].position.val()));
  double sigrel_eff =
      cal.efficiency.SigmaRel(cal.channel_to_energy(peaks[PeakIndex].position.val()));
  //sigrel_eff = 0 if uninitialized
  return (square(std::sqrt(std::sqrt(t) / t)) + square(sigrel_eff)) *
      (t / eff) * std::max(1.0, ChisqNorm());
}

size_t Region::FitVars() const
{
  size_t n = 2;    //BLN,DEL: always on!
  if (AST.to_fit)
    n += 1; //AST,BST
  if (BST.to_fit)
    n += 1;
  if (TailFlag)
    n += 2; //ALT,BLT
  if (RightTailFlag)
    n += 2; //ART,BRT
  if (StepBkg())
    n += 1;
  if (SlopeFlag)
    n += 1;
  if (CurveFlag)
    n += 1;
  n += 2 * peaks.size(); //GAM, POS
  return n;
}

void Region::setupFit()
{
  auto vars = FitVars(); // - 1 ?
  Vector.resize(vars);
  Gradient.resize(vars);
  //ReDim ChisqGradient(FitVars - 1)
  Hessinv.resize(vars, vars);

  int32_t shift = 0;

  Vector[0] = BLN.x();
  BLN.x_index = 0;
  Vector[1] = DEL.x();
  DEL.x_index = 1;

  if (AST.to_fit)
  {
    Vector[2] = AST.x();
    AST.x_index = 2;
    shift += 1;
  }

  if (BST.to_fit)
  {
    Vector[2 + shift] = BST.x();
    BST.x_index = 2 + shift;
    shift += 1;
  }

  if (TailFlag)
  {
    Vector[2 + shift] = ALT.x();
    ALT.x_index = 2 + shift;
    Vector[3 + shift] = BLT.x();
    BLT.x_index = 3 + shift;
    shift += 2;
  }

  if (RightTailFlag)
  {
    Vector[2 + shift] = ART.x();
    ART.x_index = 2 + shift;
    Vector[3 + shift] = BRT.x();
    BRT.x_index = 3 + shift;
    shift += 2;
  }

  if (StepBkg())
  {
    Vector[2 + shift] = SIG.x();
    SIG.x_index = 2 + shift;
    shift += 1;
  }

  if (Slope())
  {
    Vector[2 + shift] = BSL.x();
    BSL.x_index = 2 + shift;
    shift += 1;
  }

  if (Curve())
  {
    Vector[2 + shift] = BCV.x();
    BCV.x_index = 2 + shift;
    shift += 1;
  }

  for (auto& p : peaks)
  {
    Vector[2 + shift] = p.GAM.x();
    p.GAM.x_index = 2 + shift;
    Vector[3 + shift] = p.position.x();
    p.position.x_index = 3 + shift;
    shift += 2;
  }
}

void Region::storeFit()
{
  double Chisq_norm = std::max(ChisqNorm(), 1.0) * 0.5;
  int32_t shift{0};

  double df = DegreeOfFreedom();
  for (size_t i = 0; i < FitVars(); ++i)
    for (size_t j = 0; j < FitVars(); ++j)
      Hessinv.coeffRef(i, j) *= df;

  BLN.x(Vector[0]);
  BLN.uncert_value = std::sqrt(std::abs(Hessinv.coeff(0, 0) *
      BLN.grad_at(square(Vector[0])) * Chisq_norm));
  DEL.x(Vector[1]);
  DEL.uncert_value = std::sqrt(std::abs(Hessinv.coeff(1, 1) *
      DEL.grad_at(square(Vector[1])) * Chisq_norm));

  if (AST.to_fit)
  {
    AST.x(Vector[2]);
    AST.uncert_value = std::sqrt(std::abs(Hessinv.coeff(2, 2) *
        AST.grad_at(square(Vector[2])) * Chisq_norm));
    shift += 1;
  }

  if (BST.to_fit)
  {
    BST.x(Vector[2 + shift]);
    BST.uncert_value = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                              * BST.grad_at(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (TailFlag)
  {
    ALT.x(Vector[2 + shift]);
    ALT.uncert_value = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                              * ALT.grad_at(square(Vector[2 + shift])) * Chisq_norm));
    BLT.x(Vector[3 + shift]);
    BLT.uncert_value = std::sqrt(std::abs(Hessinv.coeff(3 + shift, 3 + shift)
                                              * BLT.grad_at(square(Vector[3 + shift])) * Chisq_norm));
    shift += 2;
  }

  if (RightTailFlag)
  {
    ART.x(Vector[2 + shift]);
    ART.uncert_value = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                              * ART.grad_at(square(Vector[2 + shift])) * Chisq_norm));
    BRT.x(Vector[3 + shift]);
    BRT.uncert_value = std::sqrt(std::abs(Hessinv.coeff(3 + shift, 3 + shift)
                                              * BRT.grad_at(square(Vector[3 + shift])) * Chisq_norm));
    shift += 2;
  }

  if (StepBkg())
  {
    SIG.x(Vector[2 + shift]);
    SIG.uncert_value = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                              * SIG.grad_at(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (Slope())
  {
    BSL.x(Vector[2 + shift]);
    BSL.uncert_value = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                              * BSL.grad_at(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (Curve())
  {
    BCV.x(Vector[2 + shift]);
    BCV.uncert_value = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                              * BCV.grad_at(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  for (auto& p : peaks)
  {
    p.GAM.x(Vector[2 + shift]);
    p.GAM.uncert_value =
        std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift) *
            p.GAM.grad_at(square(Vector[2 + shift])) * Chisq_norm));
    p.position.x(Vector[3 + shift]);
    p.position.uncert_value =
        std::sqrt(std::abs(Hessinv.coeff(3 + shift, 3 + shift) *
            p.position.grad_at(square(Vector[3 + shift])) * Chisq_norm));
    shift += 2;
  }
}

void Region::FuncValue(double E, std::vector<double>& ret) const
{
  //returns the value of the fitted curve and the background at Energy E
  //Dim FTotal, FBkg0, FBkg, FPeak As Double
  ret.resize(peaks.size() + 2);
  double _DE, _GAM, _DEL, _BST, _BRT, _BLT;
  try
  {
    //val(1):bkg
    ret[1] = BLN.val();
    if (Slope())
      ret[1] += BSL.val() * (E - FirstChannel);
    if (Curve())
      ret[1] += BCV.val() * square(E - FirstChannel);

    _BRT = BRT.val();
    _BST = BST.val();
    _DEL = DEL.val();
    _BLT = BLT.val();

    for (size_t i = 0; i < peaks.size(); ++i)
    {
      _DE = E - peaks[i].position.val();
      _GAM = peaks[i].GAM.val();

      if (LeftTail())
        ret[1] += _GAM * 0.5 * ALT.val() *
            std::exp(_DE / (_BLT * _DEL)) *
            std::erfc(_DE / _DEL + 0.5 / _BLT);
      if (StepBkg())
        ret[1] += SIG.val() * 0.5 * _GAM *
            std::erfc(peaks[i].step_type() * _DE / _DEL);

      ret[i + 2] = _GAM * std::exp(-1.0 * (square(_DE) / square(_DEL)));
      ret[i + 2] += _GAM * 0.5 * AST.val() *
          std::exp(_DE / (_BST * _DEL)) *
          std::erfc(_DE / _DEL + 0.5 / _BST);
      if (RightTail())
        ret[i + 2] += _GAM * 0.5 * ART.val() *
            std::exp(-1.0 * _DE / (_BRT * _DEL)) *
            std::erfc(-1.0 * _DE / _DEL + 0.5 / _BRT);
    }
    //val(0):FTotal
    ret[0] = ret[1];
    for (size_t i = 0; i < peaks.size(); ++i)
      ret[0] += ret[i + 2];
  }
  catch (...)
  {
    ERR("Expection");
  }
}

double Region::CalcChiSq(const std::vector<double>& XVector) const
{
  //Calculates the normalized Chi-square over a region
  try
  {
    Chisq = 0;

    double _AST = AST.to_fit ? AST.val_at(XVector[AST.x_index]) : AST.val();
    double _BST = BST.to_fit ? BST.val_at(XVector[BST.x_index]) : BST.val();
    double _DEL = DEL.val_at(XVector[DEL.x_index]);
    double _ALT = LeftTail() ? ALT.val_at(XVector[ALT.x_index]) : 0.0;
    double _BLT = LeftTail() ? BLT.val_at(XVector[BLT.x_index]) : 0.0;
    double _SIG = StepBkg() ? SIG.val_at(XVector[SIG.x_index]) : 0.0;
    double _ART = RightTail() ? ART.val_at(XVector[ART.x_index]) : 0.0;
    double _BRT = RightTail() ? BRT.val_at(XVector[BRT.x_index]) : 0.0;

    for (size_t j = FirstChannel; j <= LastChannel; ++j)
    {
      // Background
      double FTotal = BLN.val_at(XVector[BLN.x_index]);
      if (Slope())
        FTotal += BSL.val_at(XVector[BSL.x_index]) * (j - FirstChannel);
      if (Curve())
        FTotal += BCV.val_at(XVector[BCV.x_index]) * square(j - FirstChannel);

      for (auto& p : peaks)
      {
        double DE = j - p.position.val_at(XVector[p.position.x_index]);
        double _GAM = p.GAM.val_at(XVector[p.GAM.x_index]);
        if (LeftTail())
        {
          FTotal += _GAM * 0.5 * _ALT *
              std::exp(DE / (_BLT * _DEL)) *
              std::erfc(DE / _DEL + 0.5 / _BLT);
        }
        if (StepBkg())
        {
          FTotal += _GAM * 0.5 * _SIG *
              std::erfc(p.step_type() * DE / _DEL);
        }
        FTotal = std::max(FTotal, 0.0);
        //--- Peak components ---
        // Gaussian
        FTotal += _GAM * std::exp(-1.0 * square(DE / _DEL));
        // Short tail
        FTotal += _GAM * 0.5 * _AST *
            std::exp(DE / (_BST * _DEL)) *
            std::erfc(DE / _DEL + 0.5 / _BST);
        if (RightTail())
        {
          FTotal += _GAM * 0.5 * _ART *
              std::exp(-1.0 * DE / (_BRT * _DEL)) *
              std::erfc(0.5 / _BRT - DE / _DEL);
        }
      }
      Chisq += square((spectrum.Channel[j] - FTotal) /
          spectrum.Weight(j));
    } //Channel

    return Chisq;
  }
  catch (...)
  {
    std::throw_with_nested(std::runtime_error("CalcChiSq failed"));
  }
}

double Region::ChisqNorm() const
{
  return Chisq / ((LastChannel - FirstChannel) - FitVars());
}

size_t Region::DegreeOfFreedom() const
{
  return ((LastChannel - FirstChannel) - FitVars());
}

void Region::GradChiSq(const std::vector<double>& XVector,
                       std::vector<double>& XGradient, double& Chisq) const
{
  //Calculates the Chi-square and its gradient

  /*if(DiffType = 2)
  {
      Call dfunc2(reg, XVector, XGradient, Chisq)
      Exit Sub
  }

  if(DiffType = 3)
  {
      Call dfunc3(reg, XVector, XGradient, Chisq)
      Exit Sub
  }*/

  //Dim XGradient2(XGradient.GetLength(0) - 1) As Double, Chisq2 As Double
  //dfunc2(reg, XVector, XGradient2, Chisq2)
  try
  {
    // zero-out arrays
    std::vector<double> XXGradient(XGradient.size(), 0.0);
    XGradient.assign(XGradient.size(), 0.0);

    double t2;

    Chisq = 0;

    double _AST = AST.to_fit ? AST.val_at(XVector[AST.x_index]) : AST.val();
    double _BST = BST.to_fit ? BST.val_at(XVector[BST.x_index]) : BST.val();
    double _DEL = DEL.val_at(XVector[DEL.x_index]);
    double _ALT = LeftTail() ? ALT.val_at(XVector[ALT.x_index]) : 0.0;
    double _BLT = LeftTail() ? BLT.val_at(XVector[BLT.x_index]) : 0.0;
    double _SIG = StepBkg() ? SIG.val_at(XVector[SIG.x_index]) : 0.0;
    double _ART = RightTail() ? ART.val_at(XVector[ART.x_index]) : 0.0;
    double _BRT = RightTail() ? BRT.val_at(XVector[BRT.x_index]) : 0.0;

    for (size_t j = FirstChannel; j <= LastChannel; ++j)
    {
      //--- Poly Background ---
      double FTotal = BLN.val_at(XVector[BLN.x_index]);
      XXGradient[BLN.x_index] = BLN.grad_at(XVector[BLN.x_index]);
      if (Slope())
      {
        FTotal += BSL.val_at(XVector[BSL.x_index]) *
            (j - FirstChannel);
        XXGradient[BSL.x_index] = (j - FirstChannel);
      }

      if (Curve())
      {
        FTotal += BCV.val_at(XVector[BCV.x_index]) *
            square(j - FirstChannel);
        XXGradient[BCV.x_index] = square(j - FirstChannel);
      }

      for (auto& p : peaks)
      {

        double DE = j - p.position.val_at(XVector[p.position.x_index]);
        double _GAM = p.GAM.val_at(XVector[p.GAM.x_index]);
        double t1 = DE / _DEL;
        //---Left Tail---
        if (LeftTail())
        {
          double _LongTail = _GAM * 0.5 * _ALT * std::exp(t1 / _BLT) *
              std::erfc(t1 + 0.5 / _BLT);

          FTotal += _LongTail;

          //t2 = (_GAM * _ALT * std::exp(t1 / _BLT) / M_PI ^ (0.5) * std::exp(-(1.0 / (2.0 * _BLT) + t1) ^ 2) * t1 / _DEL)
          t2 = (_GAM * _ALT * std::exp(t1 / _BLT) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * _BLT) + t1)) / _DEL);
          XXGradient[DEL.x_index] += DEL.grad_at(XVector[DEL.x_index])
              * (-1.0 * t1 / (_DEL * _BLT) * _LongTail + t2 * t1);
          XXGradient[p.position.x_index] += -1.0 / (_BLT * _DEL) *
              _LongTail + t2;
          XXGradient[p.GAM.x_index] += _LongTail / _GAM;

          XXGradient[ALT.x_index] += _LongTail / _ALT *
              ALT.grad_at(XVector[ALT.x_index]);
          XXGradient[BLT.x_index] += BLT.grad_at(XVector[BLT.x_index])
              * ((-1.0 * t1 / square(_BLT)) *
                  _LongTail + (_DEL / (2.0 * square(_BLT)) * t2));

        }
        //---Step---
        if (StepBkg())
        {
          double _StepBkg = _SIG * 0.5 * _GAM *
              std::erfc(p.step_type() * t1);
          FTotal += _StepBkg;

          XXGradient[DEL.x_index] += DEL.grad_at(XVector[DEL.x_index]) *
              (_GAM * _SIG * p.step_type() / std::sqrt(M_PI) *
                  std::exp(-DE / _DEL * t1) * t1 / _DEL);
          XXGradient[p.GAM.x_index] += _StepBkg / _GAM;
          XXGradient[SIG.x_index] += _StepBkg / _SIG *
              SIG.grad_at(XVector[SIG.x_index]);
        }
        FTotal = std::max(FTotal, 0.0);

        //---Gaussian---
        double _Gauss = _GAM * std::exp(-1.0 * square(t1));
        FTotal += _Gauss;

        XXGradient[DEL.x_index] += DEL.grad_at(XVector[DEL.x_index]) *
            (2.0 * square(t1) / _DEL * _Gauss);

        XXGradient[p.position.x_index] += 2.0 * t1 / _DEL * _Gauss;
        XXGradient[p.GAM.x_index] += _Gauss / _GAM;

        //---Short Tail---

        double _ShortTail = _GAM * 0.5 * _AST * std::exp(t1 / _BST) *
            std::erfc(t1 + 0.5 / _BST);
        FTotal += _ShortTail;

        //t2 = (_GAM * _AST * std::exp(t1 / _BST) / M_PI ^ (0.5) * std::exp(-1.0 * (1.0 / (2.0 * _BST) + t1) ^ 2) * t1 / _DEL)
        t2 = (_GAM * _AST * std::exp(t1 / _BST) / std::sqrt(M_PI) *
            std::exp(-1.0 * square(1.0 / (2.0 * _BST) + t1)) / _DEL);
        XXGradient[DEL.x_index] += DEL.grad_at(XVector[DEL.x_index]) *
            (-1.0 * t1 / (_DEL * _BST) * _ShortTail + t2 * t1);

        XXGradient[p.position.x_index] += -1.0 / (_BST * _DEL) *
            _ShortTail + t2;
        XXGradient[p.GAM.x_index] += _ShortTail / _GAM;

        if (AST.to_fit)
          XXGradient[AST.x_index] += _ShortTail / _AST *
              AST.grad_at(XVector[AST.x_index]);
        if (BST.to_fit)
          XXGradient[BST.x_index] += BST.grad_at(XVector[BST.x_index]) *
              ((-1.0 * t1 / square(_BST)) *
                  _ShortTail + (_DEL / (2.0 * square(_BST)) * t2));

        //---Right Tail---
        if (RightTail())
        {
          double _RightTail = _GAM * 0.5 * _ART *
              std::exp(-1.0 * t1 / _BRT) *
              std::erfc(0.5 / _BRT - t1);
          FTotal += _RightTail;

          //t2 = (_GAM * _ART * std::exp(-1.0 * t1 / _BRT) / M_PI ^ (0.5) * std::exp(-(1.0 / (2.0 * _BRT) - t1) ^ 2) * t1 / _DEL)
          t2 = (_GAM * _ART * std::exp(-1.0 * t1 / _BRT) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * _BRT) - t1)) / _DEL);
          XXGradient[DEL.x_index] += DEL.grad_at(XVector[DEL.x_index]) *
              ((t1 / (_DEL * _BRT) * _RightTail - t2 * t1));

          XXGradient[p.position.x_index] += 1.0 / (_BRT * _DEL) *
              _RightTail - t2;
          XXGradient[p.GAM.x_index] +=
              _RightTail / _GAM;

          XXGradient[ART.x_index] += _RightTail / _ART *
              ART.grad_at(XVector[ART.x_index]);
          XXGradient[BRT.x_index] += BRT.grad_at(XVector[BRT.x_index])
              * ((t1 / square(_BRT)) * _RightTail + (_DEL / (2.0 *
                  square(_BRT)) * t2));
        }

        //XXGradient(DEL.x_index) *= DEL.GradAt(XVector(DEL.x_index))
        XXGradient[p.GAM.x_index] *= p.GAM.grad_at(XVector[p.GAM.x_index]);
        XXGradient[p.position.x_index] *= p.position.grad_at(XVector[p.position.x_index]);
      } //Peak

      double t3 = -2 * (spectrum.Channel[j] - FTotal) /
          square(spectrum.Weight(j));

      for (size_t k = 0; k < FitVars(); ++k)
      {
        XGradient[k] += XXGradient[k] * t3;
        XXGradient[k] = 0.0;
      }
      Chisq += square((spectrum.Channel[j] - FTotal) / spectrum.Weight(j));
    } //j //Channel
    //Chisq /= df
  }
  catch (...)
  {
    ERR("Error in GradChiSq");
  }
}

}
