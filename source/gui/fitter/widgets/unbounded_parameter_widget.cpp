#include <gui/fitter/widgets/unbounded_parameter_widget.h>

#include <gui/fitter/widgets/uncertain_double_widget.h>

#include <QBoxLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QCheckBox>

UnboundedParameterWidget::UnboundedParameterWidget(const DAQuiri::UnboundedParam& val,
                                       double spin_width, double label_width, QWidget* parent)
    : QWidget(parent)
      , original_(val)
      , current_(val)
{
  QHBoxLayout *hl = new QHBoxLayout();

  modified_ = new QLabel();
  hl->addWidget(modified_);

  to_fit_ = new QCheckBox();
  to_fit_->setText("To fit");
  hl->addWidget(to_fit_);

  val_ = new QDoubleSpinBox();
  val_->setFixedHeight(25);
  val_->setFixedWidth(spin_width);
  val_->setMinimum(-std::numeric_limits<double>::max());
  val_->setMaximum(std::numeric_limits<double>::max());
  val_->setDecimals(6);
  hl->addWidget(val_);

  unc_val_ = new UncertainDoubleWidget(label_width);
  hl->addWidget(unc_val_);

  setLayout(hl);
  layout()->setContentsMargins(0, 0, 0, 0);

  update_values();

  connect(to_fit_, SIGNAL(clicked()), this, SLOT(to_fit_clicked()));
  connect(val_, SIGNAL(valueChanged(double)), this, SLOT(valChanged(double)));
}

bool UnboundedParameterWidget::changed() const
{
  if (current_.to_fit != original_.to_fit)
    return true;
  return  (current_.val() != original_.val());
}

void UnboundedParameterWidget::set_decimals(int decimals)
{
  val_->setDecimals(decimals);
}

void UnboundedParameterWidget::set_step(double step)
{
  val_->setSingleStep(step);
}

void UnboundedParameterWidget::update_values()
{
  modified_->setText(changed() ? "*" : "");

  to_fit_->setChecked(current_.to_fit);

  val_->setValue(current_.val());

  unc_val_->update(UncertainDouble(current_.val(), current_.uncert()));
}

void UnboundedParameterWidget::to_fit_clicked()
{
  current_.to_fit = to_fit_->isChecked();
  update_values();
  emit updated();
}

void UnboundedParameterWidget::valChanged(double v)
{
  current_.val(val_->value());
  update_values();
  emit updated();
}
