#include <gui/fitter/FormPeakInfo.h>
#include "ui_FormPeakInfo.h"
#include <gui/widgets/qt_util.h>

#include <core/util/custom_logger.h>

FormPeakInfo::FormPeakInfo(DAQuiri::Peak& hm, QWidget* parent)
  : QDialog(parent)
  , ui(new Ui::FormPeakInfo)
  , hm_(hm)
{
  ui->setupUi(this);
  this->setFixedSize(this->size());

//  ui->labelCaliEnergy->setText(QString::number(hm_.cali_nrg_)));
//  ui->labelCaliFWHM->setText(QString::number(hm_.cali_fwhm_)));

  // \todo reenable uncertainties

  ui->labelCenter->setText(QS(hm_.peak_position().to_string()));
  ui->labelCenterPercent->setText(QS(hm_.peak_position().error_percent_fancy()));
  UncertainDouble amp {hm_.amplitude.val(), hm_.amplitude.val_uncert_};
  ui->labelAmplitude->setText(QS(amp.to_string()));
  ui->labelAmplitudePercent->setText(QS(amp.error_percent_fancy()));
  ui->labelWidth->setText(QS(hm_.fwhm().to_string()));
  ui->labelWidthPercent->setText(QS(hm_.fwhm().error_percent_fancy()));

  UncertainDouble step {hm_.step.amplitude.val(), hm_.step.amplitude.val_uncert_};
  ui->labelStep->setText(QS(step.to_string()));
  ui->labelStepPercent->setText(QS(step.error_percent_fancy()));

  UncertainDouble tail_amp {hm_.long_tail.amplitude.val(), hm_.long_tail.amplitude.val_uncert_};
  ui->labelTailH->setText(QS(tail_amp.to_string()));
  ui->labelTailHPercent->setText(QS(tail_amp.error_percent_fancy()));
  UncertainDouble tail_slope {hm_.long_tail.slope.val(), hm_.long_tail.slope.val_uncert_};
  ui->labelTailS->setText(QS(tail_slope.to_string()));
  ui->labelTailSPercent->setText(QS(tail_slope.error_percent_fancy()));

  UncertainDouble lskew_amp {hm_.short_tail.amplitude.val(), hm_.short_tail.amplitude.val_uncert_};
  ui->labelLskewH->setText(QS(lskew_amp.to_string()));
  ui->labelLskewHPercent->setText(QS(lskew_amp.error_percent_fancy()));
  UncertainDouble lskew_slope {hm_.short_tail.slope.val(), hm_.short_tail.slope.val_uncert_};
  ui->labelLskewS->setText(QS(lskew_slope.to_string()));
  ui->labelLskewSPercent->setText(QS(lskew_slope.error_percent_fancy()));

  UncertainDouble rskew_amp {hm_.right_tail.amplitude.val(), hm_.right_tail.amplitude.val_uncert_};
  ui->labelRskewH->setText(QS(rskew_amp.to_string()));
  ui->labelRskewHPercent->setText(QS(rskew_amp.error_percent_fancy()));
  UncertainDouble rskew_slope {hm_.right_tail.slope.val(), hm_.right_tail.slope.val_uncert_};
  ui->labelRskewS->setText(QS(rskew_slope.to_string()));
  ui->labelRskewSPercent->setText(QS(rskew_slope.error_percent_fancy()));

  // \todo min+max instead of epsilon for these
  ui->doubleCenter->setValue(hm_.position.val());
  ui->doubleCenterEpsilon->setValue(hm_.position.val_uncert_);
  ui->doubleAmplitude->setValue(hm_.amplitude.val());
  ui->doubleAmplitudeEpsilon->setValue(hm_.amplitude.val_uncert_);
  ui->doubleWidth->setValue(hm_.width_.val());
  ui->doubleWidthEpsilon->setValue(hm_.width_.val_uncert_);

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
  hm_.width_.val(ui->doubleWidth->value());

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
