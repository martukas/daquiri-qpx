#include <optimizerBFGS/Region.h>

#include <core/util/custom_logger.h>

Region::Region(CSpectrum& spe, double FromChannel, double ToChannel)
: spectrum(spe)
, FirstChannel(std::min(FromChannel, ToChannel))
, LastChannel(std::max(FromChannel, ToChannel))
{
  try
  {
    DEL.Max(4);
    DEL.Min(0.8);

    AST.Max(1.5);
    AST.Min(0.02);
    AST.ToFit = true;

    BST.Max(0.5);
    BST.Min(0.2);
    BST.ToFit = true;

    ALT.Max(0.15);
    ALT.Min(0.0001);
    ALT.ToFit = true;

    BLT.Max(50);
    BLT.Min(2.5);
    BLT.ToFit = true;

    ART.Max(0.9);
    ART.Min(0.01);
    ART.ToFit = true;

    BRT.Max(1.5);
    BRT.Min(0.3);
    BRT.ToFit = true;

    SIG.Max(0.05);
    SIG.Min(0.000001);
    SIG.ToFit = true;

    BSL.ToFit = true;
    BCV.ToFit = true;

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

void Region::LeftTail(bool Value)
{
  TailFlag = Value;
  ALT.ToFit = Value;
  BLT.ToFit = Value;
}

bool Region::RightTail() const
{
  return RightTailFlag;
}
void Region::RightTail(bool Value)
{
  RightTailFlag = Value;
  ART.ToFit = Value;
  BRT.ToFit = Value;
}

bool Region::Slope() const
{
  return SlopeFlag;
}

void Region::Slope(bool Value)
{
  SlopeFlag = Value;
  BSL.ToFit = Value;
}

bool Region::Curve() const
{
  return CurveFlag;
}

void Region::Curve(bool Value)
{
  CurveFlag = Value;
  BCV.ToFit = Value;
}

bool Region::StepBkg() const
{
  //Get/Set Step Status
  return StepFlag;
}

void Region::StepBkg(bool Value)
{
  StepFlag = Value;
  SIG.ToFit = Value;
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
    m = static_cast<int32_t>(1.6551 * DEL.Value());
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
      double Value = 0;
      double Var = 0;
      for (i = j - m; i <= (j + 2 * m - 1); ++j)
        if (i > 1)
          Value = Value + L(i, j, m) * spectrum.Channel[i];
      Var += square(L(i, j, m)) * spectrum.Channel[i];

      //Conv(j - FirstChannel) = Value / std::sqrt(Var)
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
    Peak.push_back(CPeak());

    Peak.back().GAM.XIndex = Peak[std::max(Peak.size() - 2, size_t(0))].GAM.XIndex + 2;
    Peak.back().POS.XIndex = Peak[std::max(Peak.size() - 2, size_t(0))].POS.XIndex + 2;

    Peak.back().POS.Min(Min);
    Peak.back().POS.Max(Max);
    Peak.back().POS.Value(Position);
    Peak.back().POS.UncValue = 0;

    Peak.back().GAM.Value(Gamma);
    Peak.back().GAM.UncValue = 0;

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
    if (index >= Peak.size())
    {
      ERR("Can't delete the peak! (invalid index)");
      return;
    }
    if (Peak.size() == 1)
    {
      ERR("Can't delete the only peak!");
      return;
    }
    for (size_t i = index; i < Peak.size() - 1; ++i)
    {
      Peak[i].POS = Peak[i + 1].POS;
      Peak[i].GAM = Peak[i + 1].GAM;
    } //i
    Peak.resize(Peak.size() - 1);
  }
  catch (...)
  {
    ERR("Delete Peak failed!");
  }
}

void Region::SortPeak()
{
  try
  {
    //Peak.s
  }
  catch (...)
  {
    ERR("Sort Peaks failed!");
  }
}

double Region::PeakArea(size_t PeakIndex) const
{
  return Peak[PeakIndex].GAM.Value() * DEL.Value() * (std::sqrt(M_PI) +
      AST.Value() * BST.Value() * std::exp(-0.25 / square(BST.Value())) +
      ART.Value() * BRT.Value() * std::exp(-0.25 / square(BRT.Value())));
}

double Region::UncPeakArea(size_t PeakIndex)
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
  //Dim Bracket As Double = (std::sqrt(M_PI) + AST.Value() * BST.Value() * std::exp(-0.25 / BST.Value() ^ 2) + ART.Value() * BRT.Value() * std::exp(-0.25 / BRT.Value() ^ 2))
  //(dGAM/dX*dArea/dGAM)^2*Var(X)*Chisq
  //t = (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.Value() * Bracket) ^ 2 * Hessinv.coeff(DEL.XIndex, DEL.XIndex) * cs
  //(dGAM/dX*dArea/dGAM)*(dDEL/dY*dArea/dDEL)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.Value() * Bracket) * (DEL.GradAt(DEL.X) * Peak(PeakIndex).GAM.Value() * Bracket) * Hessinv.coeff(Peak(PeakIndex).GAM.XIndex, DEL.XIndex) * cs
  //if(AST.ToFit = true) {
  ////(dGAM/dX*dArea/dGAM)*(dAST/dY*dArea/dAST)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.Value() * Bracket) * (AST.GradAt(AST.X) * Peak(PeakIndex).GAM.Value() * DEL.Value() * BST.Value() * std::exp(-0.25 / BST.Value() ^ 2)) * Hessinv.coeff(Peak(PeakIndex).GAM.XIndex, AST.XIndex) * cs
  //}
  //if(BST.ToFit = true) {
  //(dGAM/dX*dArea/dGAM)*(dBST/dY*dArea/dBST)*Covar(X,Y)*Chisq
  //t += (Peak(PeakIndex).GAM.GradAt(Peak(PeakIndex).GAM.X) * DEL.Value() * Bracket) * (BST.GradAt(BST.X) * Peak(PeakIndex).GAM.Value() * DEL.Value() * AST.Value() * (1 + 0.5 / BST.Value() ^ 2) * std::exp(-0.25 / BST.Value() ^ 2)) * Hessinv.coeff(Peak(PeakIndex).GAM.XIndex, AST.XIndex) * cs
  //}
  return std::sqrt(t) * std::max(1.0, ChisqNorm());
}

double Region::PeakAreaEff(size_t PeakIndex, const CCalibration& cal)
{
  double
      eff = cal.Efficiency.Value(cal.ChannelToEnergy(Peak[PeakIndex].POS.Value()));
  //eff = 1 if uninitialized
  return (Peak[PeakIndex].GAM.Value() * DEL.Value() * (std::sqrt(M_PI) +
      AST.Value() * BST.Value() * std::exp(-0.25 / square(BST.Value())) +
      ART.Value() * BRT.Value() * std::exp(-0.25 / square(BRT.Value())))) / eff;
}

double Region::UncPeakAreaEff(size_t PeakIndex, const CCalibration& cal)
{
  double t = PeakArea(PeakIndex);
  double
      eff = cal.Efficiency.Value(cal.ChannelToEnergy(Peak[PeakIndex].POS.Value()));
  double sigrel_eff =
      cal.Efficiency.SigmaRel(cal.ChannelToEnergy(Peak[PeakIndex].POS.Value()));
  //sigrel_eff = 0 if uninitialized
  return (square(std::sqrt(std::sqrt(t) / t)) + square(sigrel_eff)) *
      (t / eff) * std::max(1.0, ChisqNorm());
}

size_t Region::FitVars() const
{
  size_t n = 2;    //BLN,DEL: always on!
  if (AST.ToFit)
    n += 1; //AST,BST
  if (BST.ToFit)
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
  n += 2 * Peak.size(); //GAM, POS
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

  Vector[0] = BLN.X();
  BLN.XIndex = 0;
  Vector[1] = DEL.X();
  DEL.XIndex = 1;

  if (AST.ToFit)
  {
    Vector[2] = AST.X();
    AST.XIndex = 2;
    shift += 1;
  }

  if (BST.ToFit)
  {
    Vector[2 + shift] = BST.X();
    BST.XIndex = 2 + shift;
    shift += 1;
  }

  if (TailFlag)
  {
    Vector[2 + shift] = ALT.X();
    ALT.XIndex = 2 + shift;
    Vector[3 + shift] = BLT.X();
    BLT.XIndex = 3 + shift;
    shift += 2;
  }

  if (RightTailFlag)
  {
    Vector[2 + shift] = ART.X();
    ART.XIndex = 2 + shift;
    Vector[3 + shift] = BRT.X();
    BRT.XIndex = 3 + shift;
    shift += 2;
  }

  if (StepBkg())
  {
    Vector[2 + shift] = SIG.X();
    SIG.XIndex = 2 + shift;
    shift += 1;
  }

  if (Slope())
  {
    Vector[2 + shift] = BSL.X();
    BSL.XIndex = 2 + shift;
    shift += 1;
  }

  if (Curve())
  {
    Vector[2 + shift] = BCV.X();
    BCV.XIndex = 2 + shift;
    shift += 1;
  }

  for (auto& p : Peak)
  {
    Vector[2 + shift] = p.GAM.X();
    p.GAM.XIndex = 2 + shift;
    Vector[3 + shift] = p.POS.X();
    p.POS.XIndex = 3 + shift;
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

  BLN.X(Vector[0]);
  BLN.UncValue = std::sqrt(std::abs(Hessinv.coeff(0, 0) *
      BLN.GradAt(square(Vector[0])) * Chisq_norm));
  DEL.X(Vector[1]);
  DEL.UncValue = std::sqrt(std::abs(Hessinv.coeff(1, 1) *
      DEL.GradAt(square(Vector[1])) * Chisq_norm));

  if (AST.ToFit)
  {
    AST.X(Vector[2]);
    AST.UncValue = std::sqrt(std::abs(Hessinv.coeff(2, 2) *
        AST.GradAt(square(Vector[2])) * Chisq_norm));
    shift += 1;
  }

  if (BST.ToFit)
  {
    BST.X(Vector[2 + shift]);
    BST.UncValue = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                          * BST.GradAt(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (TailFlag)
  {
    ALT.X(Vector[2 + shift]);
    ALT.UncValue = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                          * ALT.GradAt(square(Vector[2 + shift])) * Chisq_norm));
    BLT.X(Vector[3 + shift]);
    BLT.UncValue = std::sqrt(std::abs(Hessinv.coeff(3 + shift, 3 + shift)
                                          * BLT.GradAt(square(Vector[3 + shift])) * Chisq_norm));
    shift += 2;
  }

  if (RightTailFlag)
  {
    ART.X(Vector[2 + shift]);
    ART.UncValue = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                          * ART.GradAt(square(Vector[2 + shift])) * Chisq_norm));
    BRT.X(Vector[3 + shift]);
    BRT.UncValue = std::sqrt(std::abs(Hessinv.coeff(3 + shift, 3 + shift)
                                          * BRT.GradAt(square(Vector[3 + shift])) * Chisq_norm));
    shift += 2;
  }

  if (StepBkg())
  {
    SIG.X(Vector[2 + shift]);
    SIG.UncValue = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                          * SIG.GradAt(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (Slope())
  {
    BSL.X(Vector[2 + shift]);
    BSL.UncValue = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                          * BSL.GradAt(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  if (Curve())
  {
    BCV.X(Vector[2 + shift]);
    BCV.UncValue = std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift)
                                          * BCV.GradAt(square(Vector[2 + shift])) * Chisq_norm));
    shift += 1;
  }

  for (auto& p : Peak)
  {
    p.GAM.X(Vector[2 + shift]);
    p.GAM.UncValue =
        std::sqrt(std::abs(Hessinv.coeff(2 + shift, 2 + shift) *
            p.GAM.GradAt(square(Vector[2 + shift])) * Chisq_norm));
    p.POS.X(Vector[3 + shift]);
    p.POS.UncValue =
        std::sqrt(std::abs(Hessinv.coeff(3 + shift, 3 + shift) *
            p.POS.GradAt(square(Vector[3 + shift])) * Chisq_norm));
    shift += 2;
  }
}

void Region::FuncValue(double E, std::vector<double>& Value)
{
  //returns the value of the fitted curve and the background at Energy E
  //Dim FTotal, FBkg0, FBkg, FPeak As Double
  Value.resize(Peak.size() + 2);
  double _DE, _GAM, _DEL, _BST, _BRT, _BLT;
  try
  {
    //Value(1):bkg
    Value[1] = BLN.Value();
    if (Slope())
      Value[1] += BSL.Value() * (E - FirstChannel);
    if (Curve())
      Value[1] += BCV.Value() * square(E - FirstChannel);

    _BRT = BRT.Value();
    _BST = BST.Value();
    _DEL = DEL.Value();
    _BLT = BLT.Value();

    for (size_t i = 0; i < Peak.size(); ++i)
    {
      _DE = E - Peak[i].POS.Value();
      _GAM = Peak[i].GAM.Value();

      if (LeftTail())
        Value[1] += _GAM * 0.5 * ALT.Value() *
            std::exp(_DE / (_BLT * _DEL)) *
            erfc(_DE / _DEL + 0.5 / _BLT);
      if (StepBkg())
        Value[1] += SIG.Value() * 0.5 * _GAM *
            erfc(Peak[i].StepType() * _DE / _DEL);

      Value[i + 2] = _GAM * std::exp(-1.0 * (square(_DE) / square(_DEL)));
      Value[i + 2] += _GAM * 0.5 * AST.Value() *
          std::exp(_DE / (_BST * _DEL)) *
          erfc(_DE / _DEL + 0.5 / _BST);
      if (RightTail())
        Value[i + 2] += _GAM * 0.5 * ART.Value() *
            std::exp(-1.0 * _DE / (_BRT * _DEL)) *
            erfc(-1.0 * _DE / _DEL + 0.5 / _BRT);
    }
    //Value(0):FTotal
    Value[0] = Value[1];
    for (size_t i = 0; i < Peak.size(); ++i)
      Value[0] += Value[i + 2];
  }
  catch (...)
  {
    ERR("Expection");
  }
}

double Region::CalcChiSq(const std::vector<double>& XVector) const
{
  //Calculates the normalized Chi-square over a region
  double DE;
  double _GAM, _DEL, _AST, _BST, _ART, _BRT, _ALT, _BLT, _SIG, FTotal;
  double _Gauss, _ShortTail, _LongTail, _RightTail, _StepBkg;
  try
  {
    Chisq = 0;

    if (AST.ToFit)
      _AST = AST.ValueAt(XVector[AST.XIndex]);
    else
      _AST = AST.Value();

    if (BST.ToFit)
      _BST = BST.ValueAt(XVector[BST.XIndex]);
    else
      _BST = BST.Value();

    _DEL = DEL.ValueAt(XVector[DEL.XIndex]);

    if (LeftTail())
    {
      _ALT = ALT.ValueAt(XVector[ALT.XIndex]);
      _BLT = BLT.ValueAt(XVector[BLT.XIndex]);
    }

    if (StepBkg())
      _SIG = SIG.ValueAt(XVector[SIG.XIndex]);

    if (RightTail())
    {
      _ART = ART.ValueAt(XVector[ART.XIndex]);
      _BRT = BRT.ValueAt(XVector[BRT.XIndex]);
    }

    for (size_t j = FirstChannel; j <= LastChannel; ++j)
    {

      FTotal = BLN.ValueAt(XVector[BLN.XIndex]);
      if (Slope())
        FTotal += BSL.ValueAt(XVector[BSL.XIndex]) *
            (j - FirstChannel);
      if (Curve())
        FTotal += BCV.ValueAt(XVector[BCV.XIndex]) *
            square(j - FirstChannel);

      for (auto& p : Peak)
      {
        DE = j - p.POS.ValueAt(XVector[p.POS.XIndex]);
        _GAM = p.GAM.ValueAt(XVector[p.GAM.XIndex]);
        if (LeftTail())
        {
          _LongTail = _GAM * 0.5 * _ALT *
              std::exp(DE / (_BLT * _DEL)) *
              erfc(DE / _DEL + 0.5 / _BLT);
          FTotal += _LongTail;
        }
        if (StepBkg())
        {
          _StepBkg = _SIG * 0.5 * _GAM *
              erfc(p.StepType() * DE / _DEL);
          FTotal += _StepBkg;
        }
        FTotal = std::max(FTotal, 0.0);
        //--- Peak components ---
        _Gauss = _GAM * std::exp(-1.0 * square(DE / _DEL));
        FTotal += _Gauss;
        _ShortTail = _GAM * 0.5 * _AST *
            std::exp(DE / (_BST * _DEL)) *
            erfc(DE / _DEL + 0.5 / _BST);
        FTotal += _ShortTail;
        if (RightTail())
        {
          _RightTail = _GAM * 0.5 * _ART *
              std::exp(-1.0 * DE / (_BRT * _DEL)) *
              erfc(0.5 / _BRT - DE / _DEL);
          FTotal += _RightTail;
        }
      } //i
      Chisq += square((spectrum.Channel[j] - FTotal) /
          spectrum.Weight(j));
    } //j //Channel

    return Chisq;
  }
  catch (...)
  {
    ERR("Error in func: ");
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
    std::vector<double> XXGradient(XGradient.size() - 1);

    double DE, t1, t2;
    double _GAM, _DEL, _AST, _BST, _ART, _BRT, _ALT, _BLT, _SIG, FTotal;
    double _Gauss, _ShortTail, _LongTail, _RightTail, _StepBkg;

    Chisq = 0;

    for (size_t k = 0; k <= XGradient.size(); ++k)
    {
      XGradient[k] = 0.0;
      XXGradient[k] = 0.0;
    } //k

    if (AST.ToFit)
      _AST = AST.ValueAt(XVector[AST.XIndex]);
    else
      _AST = AST.Value();

    if (BST.ToFit)
      _BST = BST.ValueAt(XVector[BST.XIndex]);
    else
      _BST = BST.Value();

    _DEL = DEL.ValueAt(XVector[DEL.XIndex]);

    if (LeftTail())
    {
      _ALT = ALT.ValueAt(XVector[ALT.XIndex]);
      _BLT = BLT.ValueAt(XVector[BLT.XIndex]);
    }

    if (StepBkg())
      _SIG = SIG.ValueAt(XVector[SIG.XIndex]);

    if (RightTail())
    {
      _ART = ART.ValueAt(XVector[ART.XIndex]);
      _BRT = BRT.ValueAt(XVector[BRT.XIndex]);
    }

    for (size_t j = FirstChannel; j <= LastChannel; ++j)
    {


      //--- Poly Background ---

      FTotal = BLN.ValueAt(XVector[BLN.XIndex]);
      XXGradient[BLN.XIndex] = BLN.GradAt(XVector[BLN.XIndex]);
      if (Slope())
      {
        FTotal += BSL.ValueAt(XVector[BSL.XIndex]) *
            (j - FirstChannel);
        XXGradient[BSL.XIndex] = (j - FirstChannel);
      }

      if (Curve())
      {
        FTotal += BCV.ValueAt(XVector[BCV.XIndex]) *
            square(j - FirstChannel);
        XXGradient[BCV.XIndex] = square(j - FirstChannel);
      }

      for (auto& p : Peak)
      {

        DE = j - p.POS.ValueAt(XVector[p.POS.XIndex]);
        _GAM = p.GAM.ValueAt(XVector[p.GAM.XIndex]);
        t1 = DE / _DEL;
        //---Left Tail---
        if (LeftTail())
        {
          _LongTail = _GAM * 0.5 * _ALT * std::exp(t1 / _BLT) *
              erfc(t1 + 0.5 / _BLT);

          FTotal += _LongTail;

          //t2 = (_GAM * _ALT * std::exp(t1 / _BLT) / M_PI ^ (0.5) * std::exp(-(1.0 / (2.0 * _BLT) + t1) ^ 2) * t1 / _DEL)
          t2 = (_GAM * _ALT * std::exp(t1 / _BLT) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * _BLT) + t1)) / _DEL);
          XXGradient[DEL.XIndex] += DEL.GradAt(XVector[DEL.XIndex])
              * (-1.0 * t1 / (_DEL * _BLT) * _LongTail + t2 * t1);
          XXGradient[p.POS.XIndex] += -1.0 / (_BLT * _DEL) *
              _LongTail + t2;
          XXGradient[p.GAM.XIndex] += _LongTail / _GAM;

          XXGradient[ALT.XIndex] += _LongTail / _ALT *
              ALT.GradAt(XVector[ALT.XIndex]);
          XXGradient[BLT.XIndex] += BLT.GradAt(XVector[BLT.XIndex])
              * ((-1.0 * t1 / square(_BLT)) *
                  _LongTail + (_DEL / (2.0 * square(_BLT)) * t2));

        }
        //---Step---
        if (StepBkg())
        {
          _StepBkg = _SIG * 0.5 * _GAM *
              erfc(p.StepType() * t1);
          FTotal += _StepBkg;

          XXGradient[DEL.XIndex] += DEL.GradAt(XVector[DEL.XIndex]) *
              (_GAM * _SIG * p.StepType() / std::sqrt(M_PI) *
                  std::exp(-DE / _DEL * t1) * t1 / _DEL);
          XXGradient[p.GAM.XIndex] += _StepBkg / _GAM;
          XXGradient[SIG.XIndex] += _StepBkg / _SIG *
              SIG.GradAt(XVector[SIG.XIndex]);
        }
        FTotal = std::max(FTotal, 0.0);

        //---Gaussian---
        _Gauss = _GAM * std::exp(-1.0 * square(t1));
        FTotal += _Gauss;

        XXGradient[DEL.XIndex] += DEL.GradAt(XVector[DEL.XIndex]) *
            (2.0 * square(t1) / _DEL * _Gauss);

        XXGradient[p.POS.XIndex] += 2.0 * t1 / _DEL * _Gauss;
        XXGradient[p.GAM.XIndex] += _Gauss / _GAM;

        //---Short Tail---

        _ShortTail = _GAM * 0.5 * _AST * std::exp(t1 / _BST) *
            erfc(t1 + 0.5 / _BST);
        FTotal += _ShortTail;

        //t2 = (_GAM * _AST * std::exp(t1 / _BST) / M_PI ^ (0.5) * std::exp(-1.0 * (1.0 / (2.0 * _BST) + t1) ^ 2) * t1 / _DEL)
        t2 = (_GAM * _AST * std::exp(t1 / _BST) / std::sqrt(M_PI) *
            std::exp(-1.0 * square(1.0 / (2.0 * _BST) + t1)) / _DEL);
        XXGradient[DEL.XIndex] += DEL.GradAt(XVector[DEL.XIndex]) *
            (-1.0 * t1 / (_DEL * _BST) * _ShortTail + t2 * t1);

        XXGradient[p.POS.XIndex] += -1.0 / (_BST * _DEL) *
            _ShortTail + t2;
        XXGradient[p.GAM.XIndex] += _ShortTail / _GAM;

        if (AST.ToFit)
          XXGradient[AST.XIndex] += _ShortTail / _AST *
              AST.GradAt(XVector[AST.XIndex]);
        if (BST.ToFit)
          XXGradient[BST.XIndex] += BST.GradAt(XVector[BST.XIndex]) *
              ((-1.0 * t1 / square(_BST)) *
                  _ShortTail + (_DEL / (2.0 * square(_BST)) * t2));

        //---Right Tail---
        if (RightTail())
        {

          _RightTail = _GAM * 0.5 * _ART *
              std::exp(-1.0 * t1 / _BRT) *
              erfc(0.5 / _BRT - t1);
          FTotal += _RightTail;

          //t2 = (_GAM * _ART * std::exp(-1.0 * t1 / _BRT) / M_PI ^ (0.5) * std::exp(-(1.0 / (2.0 * _BRT) - t1) ^ 2) * t1 / _DEL)
          t2 = (_GAM * _ART * std::exp(-1.0 * t1 / _BRT) / std::sqrt(M_PI) *
              std::exp(-square(1.0 / (2.0 * _BRT) - t1)) / _DEL);
          XXGradient[DEL.XIndex] += DEL.GradAt(XVector[DEL.XIndex]) *
              ((t1 / (_DEL * _BRT) * _RightTail - t2 * t1));

          XXGradient[p.POS.XIndex] += 1.0 / (_BRT * _DEL) *
              _RightTail - t2;
          XXGradient[p.GAM.XIndex] +=
              _RightTail / _GAM;

          XXGradient[ART.XIndex] += _RightTail / _ART *
              ART.GradAt(XVector[ART.XIndex]);
          XXGradient[BRT.XIndex] += BRT.GradAt(XVector[BRT.XIndex])
              * ((t1 / square(_BRT)) * _RightTail + (_DEL / (2.0 *
                  square(_BRT)) * t2));
        }

        //XXGradient(DEL.XIndex) *= DEL.GradAt(XVector(DEL.XIndex))
        XXGradient[p.GAM.XIndex] *= p.GAM.GradAt(XVector[p.GAM.XIndex]);
        XXGradient[p.POS.XIndex] *= p.POS.GradAt(XVector[p.POS.XIndex]);
      } //Peak

      double t3 = -2 * (spectrum.Channel[j] - FTotal) /
          square(spectrum.Weight(j));

      for (size_t k = 0; k < FitVars(); ++k)
      {
        XGradient[k] += XXGradient[k] * t3;
        XXGradient[k] = 0;
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

