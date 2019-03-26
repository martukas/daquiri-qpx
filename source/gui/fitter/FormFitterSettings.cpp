#include <gui/fitter/FormFitterSettings.h>
#include "ui_FormFitterSettings.h"
#include <gui/fitter/uncertain_double_widget.h>
#include <gui/fitter/fit_parameter_widget.h>
#include <gui/fitter/peak_skew_widget.h>
#include <QCloseEvent>

#include <gui/widgets/qt_util.h>

#include <core/util/custom_logger.h>

FormFitterSettings::FormFitterSettings(DAQuiri::FitSettings& fs, QWidget* parent)
    : QDialog(parent)
      , ui(new Ui::FormFitterSettings)
      , fit_settings_(fs)
{
  backup_settings_ = fit_settings_;

  double spin_width = 120;
  double unc_width = 90;

  ui->setupUi(this);
  this->setFixedSize(this->size());

  ui->labelCaliEnergy->setText(QString::fromStdString(fit_settings_.calib.cali_nrg_.debug()));
  ui->labelCaliFWHM->setText(QString::fromStdString(fit_settings_.calib.cali_fwhm_.debug()));

  ui->spinFinderCutoffKeV->setValue(fit_settings_.finder_cutoff_kev);

  ui->spinKONwidth->setValue(fit_settings_.kon_settings.width);
  ui->doubleKONsigmaSpectrum->setValue(fit_settings_.kon_settings.sigma_spectrum);
  ui->doubleKONsigmaResid->setValue(fit_settings_.kon_settings.sigma_resid);
  ui->doubleKONEdgeWidthFactor->setValue(fit_settings_.kon_settings.edge_width_factor);

  ui->spinRegionMaxPeaks->setValue(fit_settings_.ROI_max_peaks);
  ui->doubleRegionExtendBackground->setValue(fit_settings_.ROI_extend_background);

  ui->spinEdgeSamples->setValue(fit_settings_.background_edge_samples);

  ui->checkResidAuto->setChecked(fit_settings_.resid_auto);
  ui->spinResidMaxIterations->setValue(fit_settings_.resid_max_iterations);
  ui->spinResidMinAmplitude->setValue(fit_settings_.resid_min_amplitude);
  ui->doubleResidTooClose->setValue(fit_settings_.resid_too_close);

  ui->checkSmallSimplify->setChecked(fit_settings_.small_simplify);
  ui->spinSmallMaxAmplitude->setValue(fit_settings_.small_max_amplitude);

  ui->checkWidthCommon->setChecked(fit_settings_.width_common);
  ui->doubleMinWidthCommon->setValue(fit_settings_.width_common_bounds.lower());
  ui->doubleMaxWidthCommon->setValue(fit_settings_.width_common_bounds.upper());

  ui->checkWidthAt511->setChecked(fit_settings_.width_at_511_variable);
  ui->spinWidthAt511Tolerance->setValue(fit_settings_.width_at_511_tolerance);

  ui->spinFitterMaxIterations->setValue(fit_settings_.fitter_max_iter);

//  on_checkGaussOnly_clicked();

  ui->FWEnergyLabel->setText("FWHM (" + QS(fit_settings_.calib.cali_nrg_.to().units) + ") ");

  width_ = new FitParameterWidget(fit_settings_.default_peak.width, spin_width, unc_width);
  ui->horizontalWidth->addWidget(width_);
  connect(width_, SIGNAL(updated()), this, SLOT(update()));
  fwhm_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWHM->addWidget(fwhm_);
  fwhm_energy_ = new UncertainDoubleWidget(unc_width);
  ui->horizontalFWEnergy->addWidget(fwhm_energy_);

  ui->checkStepEnable->setChecked(fit_settings_.default_peak.step.enabled);
  step_amp_ = new FitParameterWidget(fit_settings_.default_peak.step.amplitude, spin_width, unc_width);
  ui->horizontalStepAmp->addWidget(step_amp_);
  connect(step_amp_, SIGNAL(updated()), this, SLOT(update()));

  left_skew_ = new TailWidget("Left skew", fit_settings_.default_peak.short_tail, spin_width, unc_width);
  ui->horizontalLeftSkew->addWidget(left_skew_);
  connect(left_skew_, SIGNAL(updated()), this, SLOT(update()));

  right_skew_ = new TailWidget("Right skew", fit_settings_.default_peak.right_tail, spin_width, unc_width);
  ui->horizontalRightSkew->addWidget(right_skew_);
  connect(right_skew_, SIGNAL(updated()), this, SLOT(update()));

  tail_ = new TailWidget("Long tail", fit_settings_.default_peak.long_tail, spin_width, unc_width);
  ui->horizontalTail->addWidget(tail_);
  connect(tail_, SIGNAL(updated()), this, SLOT(update()));

  update();
}

void FormFitterSettings::on_buttonBox_accepted()
{
  fit_settings_.finder_cutoff_kev = ui->spinFinderCutoffKeV->value();

  fit_settings_.kon_settings.width = ui->spinKONwidth->value();
  fit_settings_.kon_settings.sigma_spectrum = ui->doubleKONsigmaSpectrum->value();
  fit_settings_.kon_settings.sigma_resid = ui->doubleKONsigmaResid->value();
  fit_settings_.kon_settings.edge_width_factor = ui->doubleKONEdgeWidthFactor->value();

  fit_settings_.ROI_max_peaks = ui->spinRegionMaxPeaks->value();
  fit_settings_.ROI_extend_background = ui->doubleRegionExtendBackground->value();

  fit_settings_.background_edge_samples = ui->spinEdgeSamples->value();

  fit_settings_.resid_auto = ui->checkResidAuto->isChecked();
  fit_settings_.resid_max_iterations = ui->spinResidMaxIterations->value();
  fit_settings_.resid_min_amplitude = ui->spinResidMinAmplitude->value();
  fit_settings_.resid_too_close = ui->doubleResidTooClose->value();

  fit_settings_.small_simplify = ui->checkSmallSimplify->isChecked();
  fit_settings_.small_max_amplitude = ui->spinSmallMaxAmplitude->value();

  fit_settings_.width_common = ui->checkWidthCommon->isChecked();
  fit_settings_.width_common_bounds.constrain(ui->doubleMinWidthCommon->value(),
                                              ui->doubleMaxWidthCommon->value());

  fit_settings_.width_at_511_variable = ui->checkWidthAt511->isChecked();
  fit_settings_.width_at_511_tolerance = ui->spinWidthAt511Tolerance->value();

  fit_settings_.fitter_max_iter = ui->spinFitterMaxIterations->value();

  accept();
}

FormFitterSettings::~FormFitterSettings()
{
  delete ui;
}

void FormFitterSettings::update()
{
  fit_settings_.default_peak.width = width_->parameter();
  fwhm_->update(fit_settings_.default_peak.fwhm());
  fwhm_energy_->update(fit_settings_.default_peak.fwhm_energy(fit_settings_.calib.cali_nrg_));

  fit_settings_.default_peak.step.enabled = ui->checkStepEnable->isChecked();
  fit_settings_.default_peak.step.amplitude = step_amp_->parameter();

  fit_settings_.default_peak.short_tail = left_skew_->tail();
  fit_settings_.default_peak.right_tail = right_skew_->tail();
  fit_settings_.default_peak.long_tail = tail_->tail();
}

void FormFitterSettings::closeEvent(QCloseEvent* event)
{
  event->accept();
}

void FormFitterSettings::on_buttonBox_rejected()
{
  fit_settings_ = backup_settings_;
  reject();
}

//void FormFitterSettings::on_checkGaussOnly_clicked()
//{
//  bool enabled = !ui->checkGaussOnly->isChecked();
//}
