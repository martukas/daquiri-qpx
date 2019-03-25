#include <gui/fitter/FormPeakInfo.h>
#include "ui_FormPeakInfo.h"
#include <gui/widgets/qt_util.h>

#include <core/util/custom_logger.h>

FormPeakInfo::FormPeakInfo(DAQuiri::Peak& hm,
                           const DAQuiri::FCalibration& cal, QWidget* parent)
  : QDialog(parent)
  , ui(new Ui::FormPeakInfo)
  , calib_(cal)
  , hm_(hm)
{
  ui->setupUi(this);
  this->setFixedSize(this->size());

  ui->labelCaliEnergy->setText(QS(calib_.cali_nrg_.debug()));
  ui->labelCaliFWHM->setText(QS(calib_.cali_fwhm_.debug()));

  auto position = hm_.peak_position();
  ui->labelCenter->setText(QString::number(position.value()));
  ui->labelCenterUncert->setText("+-" + QString::number(position.sigma()));
  ui->labelCenterPercent->setText(QS(position.error_percent_fancy()));

  auto energy = hm_.peak_energy(calib_.cali_nrg_);
  ui->labelEnergy->setText(QString::number(energy.value()));
  ui->labelEnergyUncert->setText("+-" + QString::number(energy.sigma()));
  ui->labelEnergyPercent->setText(QS(energy.error_percent_fancy()));

  UncertainDouble amp {hm_.amplitude.val(), hm_.amplitude.uncert()};
  ui->labelAmplitude->setText(QString::number(amp.value()));
  ui->labelAmplitudeUncert->setText("+-" + QString::number(amp.sigma()));
  ui->labelAmplitudePercent->setText(QS(amp.error_percent_fancy()));

  ui->checkWidthOverride->setChecked(hm_.width_override);
  ui->checkWidthFit->setChecked(hm_.width.to_fit);
  UncertainDouble width {hm_.width.val(), hm_.width.uncert()};
  ui->labelWidth->setText(QString::number(width.value()));
  ui->labelWidthUncert->setText("+-" + QString::number(width.sigma()));
  ui->labelWidthPercent->setText(QS(width.error_percent_fancy()));

  auto fwhm = hm_.fwhm();
  ui->labelFWHM->setText(QString::number(fwhm.value()));
  ui->labelFWHMUncert->setText("+-" + QString::number(fwhm.sigma()));
  ui->labelFWHMPercent->setText(QS(fwhm.error_percent_fancy()));

  auto fwhm_energy = hm_.fwhm_energy(calib_.cali_nrg_);
  ui->labelFWEnergy->setText(QString::number(fwhm_energy.value()));
  ui->labelFWEnergyUncert->setText("+-" + QString::number(fwhm_energy.sigma()));
  ui->labelFWEnergyPercent->setText(QS(fwhm_energy.error_percent_fancy()));

  UncertainDouble step {hm_.step.amplitude.val(), hm_.step.amplitude.uncert()};
  ui->labelStep->setText(QS(step.to_string()));
  ui->labelStepPercent->setText(QS(step.error_percent_fancy()));

  UncertainDouble tail_amp {hm_.long_tail.amplitude.val(), hm_.long_tail.amplitude.uncert()};
  ui->labelTailH->setText(QS(tail_amp.to_string()));
  ui->labelTailHPercent->setText(QS(tail_amp.error_percent_fancy()));
  UncertainDouble tail_slope {hm_.long_tail.slope.val(), hm_.long_tail.slope.uncert()};
  ui->labelTailS->setText(QS(tail_slope.to_string()));
  ui->labelTailSPercent->setText(QS(tail_slope.error_percent_fancy()));

  UncertainDouble lskew_amp {hm_.short_tail.amplitude.val(), hm_.short_tail.amplitude.uncert()};
  ui->labelLskewH->setText(QS(lskew_amp.to_string()));
  ui->labelLskewHPercent->setText(QS(lskew_amp.error_percent_fancy()));
  UncertainDouble lskew_slope {hm_.short_tail.slope.val(), hm_.short_tail.slope.uncert()};
  ui->labelLskewS->setText(QS(lskew_slope.to_string()));
  ui->labelLskewSPercent->setText(QS(lskew_slope.error_percent_fancy()));

  UncertainDouble rskew_amp {hm_.right_tail.amplitude.val(), hm_.right_tail.amplitude.uncert()};
  ui->labelRskewH->setText(QS(rskew_amp.to_string()));
  ui->labelRskewHPercent->setText(QS(rskew_amp.error_percent_fancy()));
  UncertainDouble rskew_slope {hm_.right_tail.slope.val(), hm_.right_tail.slope.uncert()};
  ui->labelRskewS->setText(QS(rskew_slope.to_string()));
  ui->labelRskewSPercent->setText(QS(rskew_slope.error_percent_fancy()));

  // \todo min+max instead of epsilon for these
  ui->doubleCenter->setValue(hm_.position.val());
  ui->doubleAmplitude->setValue(hm_.amplitude.val());
  ui->doubleWidth->setValue(hm_.width.val());

  ui->checkStepEnable->setChecked(hm_.step.enabled);
  ui->checkStepFixed->setChecked(!hm_.step.amplitude.to_fit);
  ui->doubleMinStep->setValue(hm_.step.amplitude.min());
  ui->doubleMaxStep->setValue(hm_.step.amplitude.max());
  ui->doubleInitStep->setValue(hm_.step.amplitude.val());

  ui->checkTailEnable->setChecked(hm_.long_tail.enabled);
  ui->checkTailFixed->setChecked(!hm_.long_tail.amplitude.to_fit); // \todo for slope
  ui->doubleMinTailAmp->setValue(hm_.long_tail.amplitude.min());
  ui->doubleMaxTailAmp->setValue(hm_.long_tail.amplitude.max());
  ui->doubleInitTailAmp->setValue(hm_.long_tail.amplitude.val());
  ui->doubleMinTailSlope->setValue(hm_.long_tail.slope.min());
  ui->doubleMaxTailSlope->setValue(hm_.long_tail.slope.max());
  ui->doubleInitTailSlope->setValue(hm_.long_tail.slope.val());

  ui->checkEnableLskew->setChecked(hm_.short_tail.enabled);
  ui->checkLskewFixed->setChecked(!hm_.short_tail.amplitude.to_fit); // \todo for slope
  ui->doubleMinLskewAmp->setValue(hm_.short_tail.amplitude.min());
  ui->doubleMaxLskewAmp->setValue(hm_.short_tail.amplitude.max());
  ui->doubleInitLskewAmp->setValue(hm_.short_tail.amplitude.val());
  ui->doubleMinLskewSlope->setValue(hm_.short_tail.slope.min());
  ui->doubleMaxLskewSlope->setValue(hm_.short_tail.slope.max());
  ui->doubleInitLskewSlope->setValue(hm_.short_tail.slope.val());

  ui->checkEnableRskew->setChecked(hm_.right_tail.enabled);
  ui->checkRskewFixed->setChecked(!hm_.right_tail.amplitude.to_fit); // \todo for slope
  ui->doubleMinRskewAmp->setValue(hm_.right_tail.amplitude.min());
  ui->doubleMaxRskewAmp->setValue(hm_.right_tail.amplitude.max());
  ui->doubleInitRskewAmp->setValue(hm_.right_tail.amplitude.val());
  ui->doubleMinRskewSlope->setValue(hm_.right_tail.slope.min());
  ui->doubleMaxRskewSlope->setValue(hm_.right_tail.slope.max());
  ui->doubleInitRskewSlope->setValue(hm_.right_tail.slope.val());

//  ui->checkWidthCommon->setChecked(hm_.width_common);
//  ui->doubleMinWidthCommon->setValue(hm_.width_common_bounds.lower());
//  ui->doubleMaxWidthCommon->setValue(hm_.width_common_bounds.upper());
//  ui->doubleMinWidthVariable->setValue(hm_.width_variable_bounds.lower());
//  ui->doubleMaxWidthVariable->setValue(hm_.width_variable_bounds.upper());

//  ui->doubleLateralSlack->setValue(hm_.lateral_slack);
//  ui->spinFitterMaxIterations->setValue(hm_.fitter_max_iter);
}

void FormPeakInfo::on_buttonBox_accepted()
{
  DAQuiri::Parameter p, s;

  hm_.position.val(ui->doubleCenter->value());
  hm_.amplitude.val(ui->doubleAmplitude->value());
  hm_.width.val(ui->doubleWidth->value());

  // \todo fixed?
  hm_.step.enabled = ui->checkStepEnable->isChecked();
  hm_.step.amplitude.val(ui->doubleInitStep->value());
  hm_.step.amplitude.bound(ui->doubleMinStep->value(),
                           ui->doubleMaxStep->value());

  hm_.long_tail.enabled = ui->checkTailEnable->isChecked();
  hm_.long_tail.amplitude.val(ui->doubleInitTailAmp->value());
  hm_.long_tail.amplitude.bound(ui->doubleMinTailAmp->value(),
                                ui->doubleMaxTailAmp->value());
  hm_.long_tail.slope.val(ui->doubleInitTailSlope->value());
  hm_.long_tail.slope.bound(ui->doubleMinTailSlope->value(),
                            ui->doubleMaxTailSlope->value());

  hm_.short_tail.enabled = ui->checkEnableLskew->isChecked();
  hm_.short_tail.amplitude.val(ui->doubleInitLskewAmp->value());
  hm_.short_tail.amplitude.bound(ui->doubleMinLskewAmp->value(),
                                ui->doubleMaxLskewAmp->value());
  hm_.short_tail.slope.val(ui->doubleInitLskewSlope->value());
  hm_.short_tail.slope.bound(ui->doubleMinLskewSlope->value(),
                            ui->doubleMaxLskewSlope->value());

  hm_.right_tail.enabled = ui->checkEnableRskew->isChecked();
  hm_.right_tail.amplitude.val(ui->doubleInitRskewAmp->value());
  hm_.right_tail.amplitude.bound(ui->doubleMinRskewAmp->value(),
                                 ui->doubleMaxRskewAmp->value());
  hm_.right_tail.slope.val(ui->doubleInitRskewSlope->value());
  hm_.right_tail.slope.bound(ui->doubleMinRskewSlope->value(),
                             ui->doubleMaxRskewSlope->value());

//  hm_.width_common = ui->checkWidthCommon->isChecked();
//  hm_.width_common_bounds.set(ui->doubleMinWidthCommon->value();
//  hm_.width_common_bounds.ubound = ui->doubleMaxWidthCommon->value();
//  hm_.width_variable_bounds.set(ui->doubleMinWidthVariable->value();
//  hm_.width_variable_bounds.ubound = ui->doubleMaxWidthVariable->value();

//  hm_.lateral_slack = ui->doubleLateralSlack->value();
//  hm_.fitter_max_iter = ui->spinFitterMaxIterations->value();

  accept();
}

void FormPeakInfo::enforce_bounds()
{
  ui->doubleInitStep->setMinimum(ui->doubleMinStep->value());
  ui->doubleMaxStep->setMinimum(ui->doubleMinStep->value());
  ui->doubleInitStep->setMaximum(ui->doubleMaxStep->value());
  ui->doubleMinStep->setMaximum(ui->doubleMaxStep->value());

  ui->doubleInitTailAmp->setMinimum(ui->doubleMinTailAmp->value());
  ui->doubleMaxTailAmp->setMinimum(ui->doubleMinTailAmp->value());
  ui->doubleInitTailAmp->setMaximum(ui->doubleMaxTailAmp->value());
  ui->doubleMinTailAmp->setMaximum(ui->doubleMaxTailAmp->value());

  ui->doubleInitTailSlope->setMinimum(ui->doubleMinTailSlope->value());
  ui->doubleMaxTailSlope->setMinimum(ui->doubleMinTailSlope->value());
  ui->doubleInitTailSlope->setMaximum(ui->doubleMaxTailSlope->value());
  ui->doubleMinTailSlope->setMaximum(ui->doubleMaxTailSlope->value());

  ui->doubleInitLskewAmp->setMinimum(ui->doubleMinLskewAmp->value());
  ui->doubleMaxLskewAmp->setMinimum(ui->doubleMinLskewAmp->value());
  ui->doubleInitLskewAmp->setMaximum(ui->doubleMaxLskewAmp->value());
  ui->doubleMinLskewAmp->setMaximum(ui->doubleMaxLskewAmp->value());

  ui->doubleInitLskewSlope->setMinimum(ui->doubleMinLskewSlope->value());
  ui->doubleMaxLskewSlope->setMinimum(ui->doubleMinLskewSlope->value());
  ui->doubleInitLskewSlope->setMaximum(ui->doubleMaxLskewSlope->value());
  ui->doubleMinLskewSlope->setMaximum(ui->doubleMaxLskewSlope->value());

  ui->doubleInitRskewAmp->setMinimum(ui->doubleMinRskewAmp->value());
  ui->doubleMaxRskewAmp->setMinimum(ui->doubleMinRskewAmp->value());
  ui->doubleInitRskewAmp->setMaximum(ui->doubleMaxRskewAmp->value());
  ui->doubleMinRskewAmp->setMaximum(ui->doubleMaxRskewAmp->value());

  ui->doubleInitRskewSlope->setMinimum(ui->doubleMinRskewSlope->value());
  ui->doubleMaxRskewSlope->setMinimum(ui->doubleMinRskewSlope->value());
  ui->doubleInitRskewSlope->setMaximum(ui->doubleMaxRskewSlope->value());
  ui->doubleMinRskewSlope->setMaximum(ui->doubleMaxRskewSlope->value());
}

FormPeakInfo::~FormPeakInfo()
{
  delete ui;
}

void FormPeakInfo::closeEvent(QCloseEvent* event)
{
  event->accept();
}

void FormPeakInfo::on_buttonBox_rejected()
{
  reject();
}

void FormPeakInfo::on_doubleMinRskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMaxRskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMinRskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMaxRskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMinLskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMaxLskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMinLskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMaxLskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMinTailSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMaxTailSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMinTailAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMaxTailAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMinStep_valueChanged(double)
{
  enforce_bounds();
}

void FormPeakInfo::on_doubleMaxStep_valueChanged(double)
{
  enforce_bounds();
}
