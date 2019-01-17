#pragma once

#include <QWidget>
#include <QAbstractTableModel>
#include <QItemSelectionModel>
#include <QSortFilterProxyModel>
#include <gui/analysis/isotope.h>

#include <gui/widgets/SettingDelegate.h>

namespace Ui {
class WidgetIsotopes;
}

class TableGammas : public QAbstractTableModel
{
  Q_OBJECT

private:
  std::vector<RadTypes::Radiation> gammas_;

public:
  void set_gammas(const Container<RadTypes::Radiation> &);
  Container<RadTypes::Radiation> get_gammas();
  void clear();

  explicit TableGammas(QObject *parent = 0);
  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;
  QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  Qt::ItemFlags flags(const QModelIndex & index) const;
  bool setData(const QModelIndex & index, const QVariant & value, int role);

  void update();

signals:
   void energiesChanged();

public slots:

};


class WidgetIsotopes : public QWidget
{
  Q_OBJECT

public:
  explicit WidgetIsotopes(QWidget *parent = 0);
  ~WidgetIsotopes();
  std::vector<double> current_gammas() const;
  std::list<RadTypes::Radiation> current_isotope_gammas() const;
  QString current_isotope() const;
  void set_current_isotope(QString);

  void push_energies(std::vector<double>);
  void select_energies(std::set<double>);

  bool save_close();
  void select_next_energy();

  void set_editable(bool);

signals:
  void energiesSelected();
  void isotopeSelected();


private slots:
  void isotopeChosen(QString);
  void on_pushSum_clicked();
  void on_pushRemove_clicked();
  void on_pushAddGamma_clicked();

  void on_pushRemoveIsotope_clicked();

  void on_pushAddIsotope_clicked();

  void selection_changed(QItemSelection, QItemSelection);
  void energies_changed();

private:
  Ui::WidgetIsotopes *ui;

  TableGammas table_gammas_;
  SettingDelegate special_delegate_;
  QSortFilterProxyModel sort_model_;

  Container<RadTypes::Isotope> isotopes_;

  std::vector<double> current_gammas_;

  bool modified_ {false};
};
