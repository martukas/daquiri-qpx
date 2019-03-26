#include <gui/fitter/FormPeakInfo.h>
#include "ui_FormPeakInfo.h"
#include <gui/widgets/qt_util.h>

#include <core/util/custom_logger.h>


#include <QBoxLayout>
#include <QLabel>

UncertainDoubleWidget::UncertainDoubleWidget(double width, QWidget* parent)
    : QWidget(parent)
{
  QHBoxLayout *hl = new QHBoxLayout();

  val_ = make_framed_label(width);
  hl->addWidget(val_);

  sigma_ = make_framed_label(width);
  hl->addWidget(sigma_);

  percent_= make_framed_label(width);
  hl->addWidget(percent_);

  setLayout(hl);
  layout()->setContentsMargins(0, 0, 0, 0);
}

QLabel* UncertainDoubleWidget::make_framed_label(double width)
{
  auto ret = new QLabel();
  ret->setFixedHeight(25);
  ret->setFixedWidth(width);
  ret->setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
  return ret;
}

void UncertainDoubleWidget::update(const UncertainDouble& val)
{
  val_->setText(QString::number(val.value()));
  sigma_->setText("+-" + QString::number(val.sigma()));
  percent_->setText(QString::number(val.error_percent())+ "%");
}

FitParameterWidget::FitParameterWidget(const DAQuiri::Value& val,
                                       double spin_width, double label_width, QWidget* parent)
    : QWidget(parent)
      , original_(val)
      , parameter_(val)
{
  QColor red;
  red.setNamedColor("#E46D59");

  QHBoxLayout *hl = new QHBoxLayout();

  modified_ = new QLabel();
  hl->addWidget(modified_);

  to_fit_ = new QCheckBox();
  to_fit_->setText("To fit");
  hl->addWidget(to_fit_);

  min_ = new QDoubleSpinBox();
  min_->setFixedHeight(25);
  min_->setFixedWidth(spin_width);
  min_->setDecimals(6);
  hl->addWidget(min_);

  QPalette pal = min_->palette();
  pal.setColor(min_->backgroundRole(), red);
  min_->setPalette(pal);

  auto lt1 = new QLabel();
  lt1->setText(" \u2264 ");
  hl->addWidget(lt1);

  val_ = new QDoubleSpinBox();
  val_->setFixedHeight(25);
  val_->setFixedWidth(spin_width);
  val_->setDecimals(6);
  hl->addWidget(val_);

  auto lt2 = new QLabel();
  lt2->setText(" \u2264 ");
  hl->addWidget(lt2);

  max_ = new QDoubleSpinBox();
  max_->setFixedHeight(25);
  max_->setFixedWidth(spin_width);
  max_->setDecimals(6);
  hl->addWidget(max_);

  pal = max_->palette();
  pal.setColor(max_->backgroundRole(), red);
  max_->setPalette(pal);

  unc_val_ = new UncertainDoubleWidget(label_width);
  hl->addWidget(unc_val_);

  setLayout(hl);
  layout()->setContentsMargins(0, 0, 0, 0);

  update();

  connect(to_fit_, SIGNAL(clicked()), this, SLOT(to_fit_clicked()));
  connect(val_, SIGNAL(valueChanged(double)), this, SLOT(valChanged(double)));
  connect(min_, SIGNAL(valueChanged(double)), this, SLOT(minChanged(double)));
  connect(max_, SIGNAL(valueChanged(double)), this, SLOT(maxChanged(double)));
}

bool FitParameterWidget::changed() const
{
  if (parameter_.to_fit != original_.to_fit)
    return true;
  if (parameter_.min() != original_.min())
    return true;
  if (parameter_.max() != original_.max())
    return true;
  return  (parameter_.val() != original_.val());
}

void FitParameterWidget::update()
{
  modified_->setText(changed() ? "*" : "");

  to_fit_->setChecked(parameter_.to_fit);

  val_->setMinimum(parameter_.min());
  val_->setMaximum(parameter_.max());
  val_->setValue(parameter_.val());

  min_->setMinimum(parameter_.min() - std::abs(parameter_.min()));
  min_->setMaximum(parameter_.max());
  min_->setValue(parameter_.min());

  max_->setMinimum(parameter_.min());
  max_->setMaximum(parameter_.max() + std::abs(parameter_.max()));
  max_->setValue(parameter_.max());

  unc_val_->update(UncertainDouble(parameter_.val(), parameter_.uncert()));
}

void FitParameterWidget::to_fit_clicked()
{
  parameter_.to_fit = to_fit_->isChecked();
  update();
  emit updated();
}

void FitParameterWidget::valChanged(double v)
{
  parameter_.val(val_->value());
  update();
  emit updated();
}

void FitParameterWidget::minChanged(double v)
{
  parameter_.bound(v, max_->value());
  parameter_.val(val_->value());
  update();
  emit updated();
}

void FitParameterWidget::maxChanged(double v)
{
  parameter_.bound(min_->value(), v);
  parameter_.val(val_->value());
  update();
  emit updated();
}


TailWidget::TailWidget(QString name, const DAQuiri::Tail& t,
                       double spin_width, double label_width, QWidget* parent)
    : QWidget(parent)
      , original_(t)
      , tail_(t)
{
  QVBoxLayout *vl = new QVBoxLayout();

  QHBoxLayout *hl1 = new QHBoxLayout();
  auto name_label = new QLabel();
  name_label->setStyleSheet("font-weight: bold");
  name_label->setText(name);
  name_label->setMinimumWidth(100);
  name_label->sizePolicy().setHorizontalPolicy(QSizePolicy::Preferred);
  hl1->addWidget(name_label);
  modified_ = new QLabel();
  modified_->sizePolicy().setHorizontalPolicy(QSizePolicy::Preferred);
  hl1->addWidget(modified_);
  override_ = new QCheckBox();
  override_->setText("Override");
  override_->setMinimumWidth(90);
  override_->sizePolicy().setHorizontalPolicy(QSizePolicy::Minimum);
  hl1->addWidget(override_);
  enabled_ = new QCheckBox();
  enabled_->setText("Enabled");
  enabled_->sizePolicy().setHorizontalPolicy(QSizePolicy::Expanding);
  hl1->addWidget(enabled_);
  vl->addLayout(hl1);

  QHBoxLayout *hl2 = new QHBoxLayout();
  hl2->insertSpacing(0, 40);
  auto amp_label = new QLabel();
  amp_label->setText("Amplitude (relative)");
  hl2->addWidget(amp_label);
  amplitude_ = new FitParameterWidget(tail_.amplitude, spin_width, label_width);
  hl2->addWidget(amplitude_);
  vl->addLayout(hl2);

  QHBoxLayout *hl3 = new QHBoxLayout();
  hl3->insertSpacing(0, 40);
  auto slp_label = new QLabel();
  slp_label->setText("Slope");
  hl3->addWidget(slp_label);
  slope_ = new FitParameterWidget(tail_.slope, spin_width, label_width);
  hl3->addWidget(slope_);
  vl->addLayout(hl3);

  setLayout(vl);
  layout()->setContentsMargins(0, 0, 0, 0);

  update();

  connect(enabled_, SIGNAL(clicked()), this, SLOT(enabled_clicked()));
  connect(override_, SIGNAL(clicked()), this, SLOT(override_clicked()));
  connect(amplitude_, SIGNAL(updated()), this, SLOT(param_changed()));
  connect(slope_, SIGNAL(updated()), this, SLOT(param_changed()));
}

bool TailWidget::changed() const
{
  if (amplitude_->changed())
    return true;
  if (slope_->changed())
    return true;
  if (tail_.enabled != original_.enabled)
    return true;
  return  (tail_.override != original_.override);
}

void TailWidget::enabled_clicked()
{
  tail_.enabled = enabled_->isChecked();
  if (changed())
    tail_.override = true;
  update();
}

void TailWidget::override_clicked()
{
  update();
}

void TailWidget::param_changed()
{
  tail_.amplitude = amplitude_->parameter();
  tail_.slope = slope_->parameter();
  if (changed())
    tail_.override = true;
  update();
}

void TailWidget::update()
{
  modified_->setText(changed() ? "*" : "");
  enabled_->setChecked(tail_.enabled);
  override_->setChecked(tail_.override);
  emit updated();
}

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

