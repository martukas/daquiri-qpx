#include <gui/fitter/widgets/peak_step_widget.h>

#include <gui/fitter/widgets/bounded_parameter_widget.h>
#include <QBoxLayout>
#include <QLabel>
#include <QCheckBox>

StepWidget::StepWidget(QString name, const DAQuiri::Step& s,
                       double spin_width, double label_width, QWidget* parent)
    : QWidget(parent)
      , original_(s)
      , current_(s)
{
  QVBoxLayout* vl = new QVBoxLayout();

  QHBoxLayout* hl1 = new QHBoxLayout();
  hl1->setSpacing(0);
  modified_ = new QLabel();
  modified_->setMinimumWidth(0);
  modified_->setMaximumWidth(8);
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

  QHBoxLayout* hl2 = new QHBoxLayout();
  hl2->insertSpacing(0, 40);
  auto amp_label = new QLabel();
  amp_label->setText("Amplitude (relative)");
  hl2->addWidget(amp_label);
  amplitude_ = new BoundedParameterWidget(current_.amplitude, spin_width, label_width);
  hl2->addWidget(amplitude_);
  vl->addLayout(hl2);

  setLayout(vl);
  layout()->setContentsMargins(0, 0, 0, 0);

  update_values();

  connect(enabled_, SIGNAL(clicked()), this, SLOT(enabled_clicked()));
  connect(override_, SIGNAL(clicked()), this, SLOT(override_clicked()));
  connect(amplitude_, SIGNAL(updated()), this, SLOT(param_changed()));
}

bool StepWidget::changed() const
{
  if (amplitude_->changed())
    return true;
  if (current_.enabled != original_.enabled)
    return true;
  return (current_.override != original_.override);
}

void StepWidget::enabled_clicked()
{
  current_.enabled = enabled_->isChecked();
  if (changed())
    current_.override = true;
  update_values();
}

void StepWidget::override_clicked()
{
  update_values();
}

void StepWidget::param_changed()
{
  current_.amplitude = amplitude_->parameter();
  if (changed())
    current_.override = true;
  update_values();
}

void StepWidget::update_values()
{
  modified_->setText(changed() ? "*" : "");
  enabled_->setChecked(current_.enabled);
  override_->setChecked(current_.override);
  emit updated();
}
