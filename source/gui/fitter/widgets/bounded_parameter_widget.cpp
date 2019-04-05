#include <gui/fitter/widgets/bounded_parameter_widget.h>

#include <gui/fitter/widgets/uncertain_double_widget.h>

#include <QBoxLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QCheckBox>

BoundedParameterWidget::BoundedParameterWidget(const DAQuiri::SineBoundedParam& val,
                                       double spin_width, double label_width, QWidget* parent)
    : QWidget(parent)
      , original_(val)
      , current_(val)
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

  update_values();

  connect(to_fit_, SIGNAL(clicked()), this, SLOT(to_fit_clicked()));
  connect(val_, SIGNAL(valueChanged(double)), this, SLOT(valChanged(double)));
  connect(min_, SIGNAL(valueChanged(double)), this, SLOT(minChanged(double)));
  connect(max_, SIGNAL(valueChanged(double)), this, SLOT(maxChanged(double)));
}

bool BoundedParameterWidget::changed() const
{
  if (current_.to_fit != original_.to_fit)
    return true;
  if (current_.min() != original_.min())
    return true;
  if (current_.max() != original_.max())
    return true;
  return  (current_.val() != original_.val());
}

void BoundedParameterWidget::update_values()
{
  modified_->setText(changed() ? "*" : "");

  to_fit_->setChecked(current_.to_fit);

  val_->setMinimum(current_.min());
  val_->setMaximum(current_.max());
  val_->setValue(current_.val());

  min_->setMinimum(current_.min() - std::abs(current_.min()));
  min_->setMaximum(current_.max());
  min_->setValue(current_.min());

  max_->setMinimum(current_.min());
  max_->setMaximum(current_.max() + std::abs(current_.max()));
  max_->setValue(current_.max());

  unc_val_->update(UncertainDouble(current_.val(), current_.uncert()));
}

void BoundedParameterWidget::to_fit_clicked()
{
  current_.to_fit = to_fit_->isChecked();
  update_values();
  emit updated();
}

void BoundedParameterWidget::valChanged(double v)
{
  current_.val(val_->value());
  update_values();
  emit updated();
}

void BoundedParameterWidget::minChanged(double v)
{
  current_.bound(v, max_->value());
  current_.val(val_->value());
  update_values();
  emit updated();
}

void BoundedParameterWidget::maxChanged(double v)
{
  current_.bound(min_->value(), v);
  current_.val(val_->value());
  update_values();
  emit updated();
}
