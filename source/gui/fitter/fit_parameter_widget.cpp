#include <gui/fitter/fit_parameter_widget.h>

#include <gui/fitter/uncertain_double_widget.h>

#include <QBoxLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QCheckBox>

FitParameterWidget::FitParameterWidget(const DAQuiri::SineBoundedParam& val,
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
