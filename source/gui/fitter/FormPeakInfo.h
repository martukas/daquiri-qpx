#pragma once

#include <QDialog>
#include <QCloseEvent>
#include <core/fitting/hypermet/Peak.h>

#include <core/fitting/fit_settings.h>

#include <QWidget>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QCheckBox>
class UncertainDoubleWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit UncertainDoubleWidget(double width, QWidget* parent = 0);
  void update(const UncertainDouble& val);

 private:
  QLabel* val_ {nullptr};
  QLabel* sigma_{nullptr};
  QLabel* percent_{nullptr};

  static QLabel* make_framed_label(double width);
};

class FitParameterWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit FitParameterWidget(const DAQuiri::Value& val,
      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  DAQuiri::Value parameter() const { return parameter_; }

 signals:
  void updated();

 private slots:
  void to_fit_clicked();
  void valChanged(double);
  void minChanged(double);
  void maxChanged(double);

 private:
  DAQuiri::Value original_;
  DAQuiri::Value parameter_;

  QLabel* modified_ {nullptr};
  QCheckBox* to_fit_{nullptr};
  QDoubleSpinBox* val_{nullptr};
  QDoubleSpinBox* min_{nullptr};
  QDoubleSpinBox* max_{nullptr};
  UncertainDoubleWidget* unc_val_ {nullptr};

  void update();
};

class TailWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit TailWidget(QString name, const DAQuiri::Tail& t,
                      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  DAQuiri::Tail tail() const { return tail_; }

 signals:
  void updated();

 private slots:
  void enabled_clicked();
  void override_clicked();
  void param_changed();

 private:
  DAQuiri::Tail original_;
  DAQuiri::Tail tail_;

  QLabel* modified_ {nullptr};
  QCheckBox* enabled_{nullptr};
  QCheckBox* override_{nullptr};
  FitParameterWidget* amplitude_ {nullptr};
  FitParameterWidget* slope_ {nullptr};

  void update();
};


namespace Ui
{
class FormPeakInfo;
}

class FormPeakInfo : public QDialog
{
 Q_OBJECT

 public:
  explicit FormPeakInfo(DAQuiri::Peak& hm, const DAQuiri::FCalibration& cal, QWidget* parent = 0);
  ~FormPeakInfo();

 protected:
  void closeEvent(QCloseEvent*);

 private slots:
  void update();

  void on_buttonBox_accepted();

  void on_buttonBox_rejected();

 private:
  Ui::FormPeakInfo* ui;

  DAQuiri::FCalibration calib_;

  DAQuiri::Peak& peak_;
  DAQuiri::Peak peak_backup_;

  FitParameterWidget* position_;
  UncertainDoubleWidget* energy_;

  UncertainDoubleWidget* area_;

  FitParameterWidget* width_;
  UncertainDoubleWidget* fwhm_;
  UncertainDoubleWidget* fwhm_energy_;

  FitParameterWidget* step_amp_;

  TailWidget* left_skew_;
  TailWidget* right_skew_;

  TailWidget* tail_;

};
