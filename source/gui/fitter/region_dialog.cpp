#include <gui/fitter/region_dialog.h>
#include "ui_region_dialog.h"
#include <gui/fitter/widgets/uncertain_double_widget.h>
#include <gui/fitter/widgets/positive_parameter_widget.h>
#include <gui/fitter/widgets/bounded_parameter_widget.h>
#include <gui/fitter/widgets/peak_skew_widget.h>
#include <gui/fitter/widgets/peak_step_widget.h>
#include <QCloseEvent>

#include <gui/widgets/qt_util.h>

#include <core/util/logger.h>

RegionDialog::RegionDialog(DAQuiri::Region& region,
                           const DAQuiri::FCalibration& calibration, QWidget* parent)
    : QDialog(parent)
      , ui(new Ui::RegionDialog)
      , calibration_(calibration)
      , region_(region)
{
  region_backup_ = region_;

  double spin_width = 120;
  double unc_width = 90;

  ui->setupUi(this);
  this->setFixedSize(this->size());

  ui->labelCaliEnergy->setText(QS(calibration_.cali_nrg_.debug()));
  ui->labelCaliFWHM->setText(QS(calibration_.cali_fwhm_.debug()));

  ui->FWEnergyLabel->setText("FWHM (" + QS(calibration_.cali_nrg_.to().units) + ") ");

  background_base_ = new BoundedParameterWidget(region_.background.base, spin_width, unc_width);
  ui->horizontalBackgroundBase->addWidget(background_base_);
  connect(background_base_, SIGNAL(updated()), this, SLOT(update()));

  background_slope_ = new BoundedParameterWidget(region_.background.slope, spin_width, unc_width);
  ui->horizontalBackgroundSlope->addWidget(background_slope_);
  connect(background_slope_, SIGNAL(updated()), this, SLOT(update()));

  background_curve_ = new BoundedParameterWidget(region_.background.curve, spin_width, unc_width);
  ui->horizontalBackgroundCurve->addWidget(background_curve_);
  connect(background_curve_, SIGNAL(updated()), this, SLOT(update()));

  width_ = new BoundedParameterWidget(region_.default_peak_.width, spin_width, unc_width);
  ui->horizontalWidth->addWidget(width_);
  connect(width_, SIGNAL(updated()), this, SLOT(update()));
  fwhm_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWHM->addWidget(fwhm_);
  fwhm_energy_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWEnergy->addWidget(fwhm_energy_);

  left_skew_ = new SkewWidget("Left skew", region_.default_peak_.left_skew, spin_width, unc_width);
  ui->horizontalLeftSkew->addWidget(left_skew_);
  connect(left_skew_, SIGNAL(updated()), this, SLOT(update()));

  right_skew_ = new SkewWidget("Right skew", region_.default_peak_.right_skew, spin_width, unc_width);
  ui->horizontalRightSkew->addWidget(right_skew_);
  connect(right_skew_, SIGNAL(updated()), this, SLOT(update()));

  step_ = new StepWidget("Step", region_.default_peak_.step, spin_width, unc_width);
  ui->horizontalStep->addWidget(step_);
  connect(step_, SIGNAL(updated()), this, SLOT(update()));

  tail_ = new SkewWidget("Tail", region_.default_peak_.tail, spin_width, unc_width);
  ui->horizontalTail->addWidget(tail_);
  connect(tail_, SIGNAL(updated()), this, SLOT(update()));

  update();
}

void RegionDialog::update()
{
  region_.background.base = background_base_->parameter();
  region_.background.slope = background_slope_->parameter();
  region_.background.curve = background_curve_->parameter();

  region_.default_peak_.width = width_->parameter();
  fwhm_->update(region_.default_peak_.fwhm());
  fwhm_energy_->update(region_.default_peak_.fwhm_energy(calibration_.cali_nrg_));
//  if (width_->changed())
//    ui->checkWidthOverride->setChecked(true);

  region_.default_peak_.left_skew = left_skew_->skew();
  region_.default_peak_.right_skew = right_skew_->skew();

  region_.default_peak_.step = step_->step();
  region_.default_peak_.tail = tail_->skew();
}

void RegionDialog::on_buttonBox_accepted()
{
  accept();
}

RegionDialog::~RegionDialog()
{
  delete ui;
}

void RegionDialog::closeEvent(QCloseEvent* event)
{
  event->accept();
}

void RegionDialog::on_buttonBox_rejected()
{
  region_ = region_backup_;
  reject();
}

