#include <gui/fitter/FormFitterSettings.h>
#include "ui_FormFitterSettings.h"

#include <core/util/custom_logger.h>

FormFitterSettings::FormFitterSettings(DAQuiri::FitSettings &fs, QWidget *parent)
  : QDialog(parent)
  , ui(new Ui::FormFitterSettings)
  , fit_settings_(fs)
{
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

  on_checkGaussOnly_clicked();
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

void FormFitterSettings::enforce_bounds()
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

FormFitterSettings::~FormFitterSettings()
{
  delete ui;
}

void FormFitterSettings::closeEvent(QCloseEvent *event) {
  event->accept();
}


void FormFitterSettings::on_buttonBox_rejected()
{
  reject();
}

void FormFitterSettings::on_doubleMinRskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMaxRskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMinRskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMaxRskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMinLskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMaxLskewSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMinLskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMaxLskewAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMinTailSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMaxTailSlope_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMinTailAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMaxTailAmp_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMinStep_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_doubleMaxStep_valueChanged(double)
{
  enforce_bounds();
}

void FormFitterSettings::on_checkGaussOnly_clicked()
{
  bool enabled = !ui->checkGaussOnly->isChecked();

  ui->checkStepEnable->setEnabled(enabled);
  ui->doubleMinStep->setEnabled(enabled);
  ui->doubleMaxStep->setEnabled(enabled);
  ui->doubleInitStep->setEnabled(enabled);

  ui->checkTailEnable->setEnabled(enabled);
  ui->doubleMinTailAmp->setEnabled(enabled);
  ui->doubleMaxTailAmp->setEnabled(enabled);
  ui->doubleInitTailAmp->setEnabled(enabled);
  ui->doubleMinTailSlope->setEnabled(enabled);
  ui->doubleMaxTailSlope->setEnabled(enabled);
  ui->doubleInitTailSlope->setEnabled(enabled);

  ui->checkEnableLskew->setEnabled(enabled);
  ui->doubleMinLskewAmp->setEnabled(enabled);
  ui->doubleMaxLskewAmp->setEnabled(enabled);
  ui->doubleInitLskewAmp->setEnabled(enabled);
  ui->doubleMinLskewSlope->setEnabled(enabled);
  ui->doubleMaxLskewSlope->setEnabled(enabled);
  ui->doubleInitLskewSlope->setEnabled(enabled);

  ui->checkEnableRskew->setEnabled(enabled);
  ui->doubleMinRskewAmp->setEnabled(enabled);
  ui->doubleMaxRskewAmp->setEnabled(enabled);
  ui->doubleInitRskewAmp->setEnabled(enabled);
  ui->doubleMinRskewSlope->setEnabled(enabled);
  ui->doubleMaxRskewSlope->setEnabled(enabled);
  ui->doubleInitRskewSlope->setEnabled(enabled);
}
