#include <gui/fitter/peak_skew_widget.h>

#include <gui/fitter/fit_parameter_widget.h>
#include <QBoxLayout>
#include <QLabel>
#include <QCheckBox>

TailWidget::TailWidget(QString name, const DAQuiri::Tail& t,
                       double spin_width, double label_width, QWidget* parent)
    : QWidget(parent)
      , original_(t)
      , tail_(t)
{
  QVBoxLayout *vl = new QVBoxLayout();

  QHBoxLayout *hl1 = new QHBoxLayout();
  hl1->setSpacing(0);
  modified_ = new QLabel();
  modified_ ->setMinimumWidth(0);
  modified_ ->setMaximumWidth(8);
  modified_->sizePolicy().setHorizontalPolicy(QSizePolicy::Maximum);
  hl1->addWidget(modified_);
  auto name_label = new QLabel();
  name_label->setStyleSheet("font-weight: bold");
  name_label->setText(name);
  name_label->setFixedWidth(100);
  name_label->sizePolicy().setHorizontalPolicy(QSizePolicy::Fixed);
  hl1->addWidget(name_label);
  override_ = new QCheckBox();
  override_->setText("Override");
  override_->setFixedWidth(90);
  override_->sizePolicy().setHorizontalPolicy(QSizePolicy::Fixed);
  hl1->addWidget(override_);
  enabled_ = new QCheckBox();
  enabled_->setText("Enabled");
  enabled_->sizePolicy().setHorizontalPolicy(QSizePolicy::MinimumExpanding);
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
