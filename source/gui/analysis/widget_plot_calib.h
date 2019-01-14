#pragma once

#include <QWidget>
#include <QPlot/QPlot1D.h>
#include <QPlot/QPlotAppearance.h>
#include <set>

class WidgetPlotCalib : public QPlot::Multi1D
{
  Q_OBJECT
private:
  struct PointSet
  {
    QPlot::Appearance appearance;
    QVector<double> x, y, x_sigma, y_sigma;
  };

public:
  explicit WidgetPlotCalib(QWidget *parent = 0);

  void clearPrimary() override;
  void replotExtras() override;
  void replotPrimary() override;

  std::set<double> get_selected_pts();
  void set_selected_pts(std::set<double>);

  void addPoints(QPlot::Appearance style,
                 const QVector<double>& x, const QVector<double>& y,
                 const QVector<double>& x_sigma, const QVector<double>& y_sigma);
  void setFit(const QVector<double>& x, const QVector<double>& y, QPlot::Appearance style);

  void set_log_x(bool);
  bool log_x() const;

private slots:
  void apply_scale_type_x();

protected:
  void executeButton(QPlot::Button *) override;

  PointSet fit_;
  QVector<PointSet> points_;

//  std::set<double> selection_;

//  QString scale_type_x_;
  bool scale_log_x_;

  void plotExtraButtons();
  void plotFit();
  void plotPoints();
};
