#include <gui/fitter/peak_dialog.h>
#include "ui_peak_dialog.h"
#include <gui/fitter/widgets/uncertain_double_widget.h>
#include <gui/fitter/widgets/positive_parameter_widget.h>
#include <gui/fitter/widgets/bounded_parameter_widget.h>
#include <gui/fitter/widgets/peak_skew_widget.h>
#include <gui/fitter/widgets/peak_step_widget.h>
#include <QCloseEvent>

#include <gui/widgets/qt_util.h>

#include <core/util/logger.h>

PeakDialog::PeakDialog(DAQuiri::Peak& peak,
                       const DAQuiri::FCalibration& calibration, QWidget* parent)
    : QDialog(parent)
      , ui(new Ui::PeakDialog)
      , calibration_(calibration)
      , peak_(peak)
{
  peak_backup_ = peak_;

  double spin_width = 120;
  double unc_width = 90;

  ui->setupUi(this);
  this->setFixedSize(this->size());

  ui->labelCaliEnergy->setText(QS(calibration_.cali_nrg_.debug()));
  ui->labelCaliFWHM->setText(QS(calibration_.cali_fwhm_.debug()));

  ui->EnergyLabel->setText("Energy (" + QS(calibration_.cali_nrg_.to().units) + ") ");
  ui->FWEnergyLabel->setText("FWHM (" + QS(calibration_.cali_nrg_.to().units) + ") ");

  position_ = new BoundedParameterWidget(peak_.position, spin_width, unc_width);
  ui->horizontalPosition->addWidget(position_);
  connect(position_, SIGNAL(updated()), this, SLOT(update()));
  energy_ = new UncertainDoubleWidget(unc_width);
  ui->energyLayout->addWidget(energy_);

  amplitude_ = new PositiveParameterWidget(peak_.amplitude, spin_width, unc_width);
  amplitude_->set_decimals(2);
  amplitude_->set_step(10);
  ui->horizontalAmplitude->addWidget(amplitude_);
  connect(amplitude_, SIGNAL(updated()), this, SLOT(update()));
  area_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalArea->addWidget(area_);

  ui->checkWidthOverride->setChecked(peak_.width_override);
  width_ = new BoundedParameterWidget(peak_.width, spin_width, unc_width);
  ui->horizontalWidth->addWidget(width_);
  connect(width_, SIGNAL(updated()), this, SLOT(update()));
  fwhm_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWHM->addWidget(fwhm_);
  fwhm_energy_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWEnergy->addWidget(fwhm_energy_);

  left_skew_ = new SkewWidget("Left skew", peak_.left_skew, spin_width, unc_width);
  ui->horizontalLeftSkew->addWidget(left_skew_);
  connect(left_skew_, SIGNAL(updated()), this, SLOT(update()));

  right_skew_ = new SkewWidget("Right skew", peak_.right_skew, spin_width, unc_width);
  ui->horizontalRightSkew->addWidget(right_skew_);
  connect(right_skew_, SIGNAL(updated()), this, SLOT(update()));

  step_ = new StepWidget("Step", peak_.step, spin_width, unc_width);
  ui->horizontalStep->addWidget(step_);
  connect(step_, SIGNAL(updated()), this, SLOT(update()));

  tail_ = new SkewWidget("Tail", peak_.tail, spin_width, unc_width);
  ui->horizontalTail->addWidget(tail_);
  connect(tail_, SIGNAL(updated()), this, SLOT(update()));

  update();
}

void PeakDialog::update()
{
  peak_.position = position_->parameter();
  energy_->update(peak_.peak_energy(calibration_.cali_nrg_));

  peak_.amplitude = amplitude_->parameter();
  area_->update(peak_.area());

  peak_.width = width_->parameter();
  fwhm_->update(peak_.fwhm());
  fwhm_energy_->update(peak_.fwhm_energy(calibration_.cali_nrg_));
  if (width_->changed())
    ui->checkWidthOverride->setChecked(true);
  peak_.width_override = ui->checkWidthOverride->isChecked();

  peak_.left_skew = left_skew_->skew();
  peak_.right_skew = right_skew_->skew();

  peak_.step = step_->step();
  peak_.tail = tail_->skew();
}

void PeakDialog::on_buttonBox_accepted()
{
  accept();
}

PeakDialog::~PeakDialog()
{
  delete ui;
}

void PeakDialog::closeEvent(QCloseEvent* event)
{
  event->accept();
}

void PeakDialog::on_buttonBox_rejected()
{
  peak_ = peak_backup_;
  reject();
}

