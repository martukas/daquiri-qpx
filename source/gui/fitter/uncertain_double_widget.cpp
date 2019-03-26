#include <gui/fitter/uncertain_double_widget.h>
#include <QBoxLayout>
#include <QLabel>

//#include <gui/widgets/qt_util.h>

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
