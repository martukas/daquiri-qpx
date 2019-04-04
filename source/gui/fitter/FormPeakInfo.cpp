#include <gui/fitter/FormPeakInfo.h>
#include "ui_FormPeakInfo.h"
#include <gui/fitter/uncertain_double_widget.h>
#include <gui/fitter/fit_parameter_widget.h>
#include <gui/fitter/peak_skew_widget.h>
#include <QCloseEvent>

#include <gui/widgets/qt_util.h>

#include <core/util/logger.h>

FormPeakInfo::FormPeakInfo(DAQuiri::Peak& hm,
                           const DAQuiri::FCalibration& cal, QWidget* parent)
  : QDialog(parent)
  , ui(new Ui::FormPeakInfo)
  , calib_(cal)
  , peak_(hm)
{
  peak_backup_ = peak_;

  double spin_width = 120;
  double unc_width = 90;

  ui->setupUi(this);
  this->setFixedSize(this->size());

  ui->labelCaliEnergy->setText(QS(calib_.cali_nrg_.debug()));
  ui->labelCaliFWHM->setText(QS(calib_.cali_fwhm_.debug()));

  ui->EnergyLabel->setText("Energy (" + QS(calib_.cali_nrg_.to().units) + ") ");
  ui->FWEnergyLabel->setText("FWHM (" + QS(calib_.cali_nrg_.to().units) + ") ");

  position_ = new FitParameterWidget(peak_.position, spin_width, unc_width);
  ui->horizontalPosition->addWidget(position_);
  connect(position_, SIGNAL(updated()), this, SLOT(update()));
  energy_ = new UncertainDoubleWidget(unc_width);
  ui->energyLayout->addWidget(energy_);

  UncertainDouble amp {peak_.amplitude.val(), peak_.amplitude.uncert()};
  ui->labelAmplitude->setText(QString::number(amp.value()));
  ui->labelAmplitudeUncert->setText("\u00B1" + QString::number(amp.sigma()));
  ui->labelAmplitudePercent->setText(QS(amp.error_percent_fancy()));
  ui->doubleAmplitude->setValue(peak_.amplitude.val());
  area_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalArea->addWidget(area_);

  ui->checkWidthOverride->setChecked(peak_.width_override);
  width_ = new FitParameterWidget(peak_.width, spin_width, unc_width);
  ui->horizontalWidth->addWidget(width_);
  connect(width_, SIGNAL(updated()), this, SLOT(update()));
  fwhm_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWHM->addWidget(fwhm_);
  fwhm_energy_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWEnergy->addWidget(fwhm_energy_);

  ui->checkStepEnable->setChecked(peak_.step.enabled);
  ui->checkStepOverride->setChecked(peak_.step.override);
  step_amp_ = new FitParameterWidget(peak_.step.amplitude, spin_width, unc_width);
  ui->horizontalStepAmp->addWidget(step_amp_);
  connect(step_amp_, SIGNAL(updated()), this, SLOT(update()));

  left_skew_ = new TailWidget("Left skew", peak_.short_tail, spin_width, unc_width);
  ui->horizontalLeftSkew->addWidget(left_skew_);
  connect(left_skew_, SIGNAL(updated()), this, SLOT(update()));

  right_skew_ = new TailWidget("Right skew", peak_.right_tail, spin_width, unc_width);
  ui->horizontalRightSkew->addWidget(right_skew_);
  connect(right_skew_, SIGNAL(updated()), this, SLOT(update()));

  tail_ = new TailWidget("Long tail", peak_.long_tail, spin_width, unc_width);
  ui->horizontalTail->addWidget(tail_);
  connect(tail_, SIGNAL(updated()), this, SLOT(update()));

  update();
}

void FormPeakInfo::update()
{
  peak_.position = position_->parameter();
  energy_->update(peak_.peak_energy(calib_.cali_nrg_));

  peak_.amplitude.val(ui->doubleAmplitude->value());
  area_->update(peak_.area());

  peak_.width_override = ui->checkWidthOverride->isChecked();
  peak_.width = width_->parameter();
  fwhm_->update(peak_.fwhm());
  fwhm_energy_->update(peak_.fwhm_energy(calib_.cali_nrg_));
  if (width_->changed())
    ui->checkWidthOverride->setChecked(true);

  peak_.step.enabled = ui->checkStepEnable->isChecked();
  peak_.step.override = ui->checkStepOverride->isChecked();
  peak_.step.amplitude = step_amp_->parameter();
  if (step_amp_->changed())
    ui->checkStepOverride->setChecked(true);

  peak_.short_tail = left_skew_->tail();
  peak_.right_tail = right_skew_->tail();
  peak_.long_tail = tail_->tail();
}


void FormPeakInfo::on_buttonBox_accepted()
{
  accept();
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
  peak_ = peak_backup_;
  reject();
}

