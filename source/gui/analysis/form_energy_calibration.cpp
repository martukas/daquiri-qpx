#include <gui/analysis/form_energy_calibration.h>
#include "ui_form_energy_calibration.h"

#include <QSettings>
//#include "widget_detectors.h"

#include <core/calibration/coef_function_factory.h>
#include <gui/widgets/qt_util.h>

FormEnergyCalibration::FormEnergyCalibration(XMLableDB<DAQuiri::Detector>& dets, DAQuiri::Fitter &fit, QWidget *parent) :
  QWidget(parent),
  ui(new Ui::FormEnergyCalibration),
  detectors_(dets),
  fit_data_(fit)
{
  ui->setupUi(this);

  loadSettings();

  QColor point_color;
  point_color.setHsv(180, 215, 150, 120);
  style_pts.default_pen = QPen(point_color, 9);
  QColor selected_color;
  selected_color.setHsv(225, 255, 230, 210);
  style_pts.themes["selected"] = QPen(selected_color, 9);
  style_fit.default_pen = QPen(Qt::darkCyan, 2);

  ui->PlotCalib->setAxisLabels("channel", "energy");


  ui->tablePeaks->verticalHeader()->hide();
  ui->tablePeaks->setColumnCount(3);
  ui->tablePeaks->setHorizontalHeaderLabels({"chan", "err(chan)", "energy"});
  ui->tablePeaks->setSelectionBehavior(QAbstractItemView::SelectRows);
  ui->tablePeaks->setSelectionMode(QAbstractItemView::ExtendedSelection);

  //all columns?
  ui->tablePeaks->setEditTriggers(QTableView::NoEditTriggers);

  ui->tablePeaks->horizontalHeader()->setStretchLastSection(true);
  ui->tablePeaks->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  ui->tablePeaks->show();
  connect(ui->tablePeaks, SIGNAL(itemSelectionChanged()), this, SLOT(selection_changed_in_table()));

  connect(ui->PlotCalib, SIGNAL(selection_changed()), this, SLOT(selection_changed_in_plot()));

  ui->isotopes->show();
  connect(ui->isotopes, SIGNAL(energiesSelected()), this, SLOT(isotope_energies_chosen()));

}

FormEnergyCalibration::~FormEnergyCalibration()
{
  delete ui;
}

bool FormEnergyCalibration::save_close() {
  if (ui->isotopes->save_close()) {
    saveSettings();
    return true;
  }
  else
    return false;
}

void FormEnergyCalibration::loadSettings()
{
  QSettings settings;

  settings.beginGroup("Program");
  ui->isotopes->setDir(settings.value("settingsdirectory", QDir::homePath() + "/qpx/settings").toString());
  data_directory_ = settings.value("save_directory", QDir::homePath() + "/qpx/data").toString();
  settings.endGroup();

  settings.beginGroup("Energy_calibration");
  ui->spinTerms->setValue(settings.value("fit_function_terms", 2).toInt());
  ui->isotopes->set_current_isotope(settings.value("current_isotope", "Co-60").toString());
  settings.endGroup();
}

void FormEnergyCalibration::saveSettings()
{
  QSettings settings;
  settings.beginGroup("Energy_calibration");
  settings.setValue("fit_function_terms", ui->spinTerms->value());
  settings.setValue("current_isotope", ui->isotopes->current_isotope());
  settings.endGroup();
}

void FormEnergyCalibration::clear()
{
  new_calibration_ = DAQuiri::Calibration();
  ui->tablePeaks->clearContents();
  ui->tablePeaks->setRowCount(0);
  toggle_push();
  ui->PlotCalib->clearAll();
  ui->PlotCalib->replot();
  ui->pushApplyCalib->setEnabled(false);
  ui->pushFromDB->setEnabled(false);
}

void FormEnergyCalibration::newSpectrum()
{
  new_calibration_ = fit_data_.settings().cali_nrg_;
  update_data();
}

void FormEnergyCalibration::update_data()
{
  rebuild_table();
  replot_calib();

  if (fit_data_.peaks().empty())
    selected_peaks_.clear();

  select_in_table();
  select_in_plot();
  toggle_push();
}

void FormEnergyCalibration::update_selection(std::set<double> selected_peaks)
{
  bool changed = (selected_peaks_ != selected_peaks);
  selected_peaks_ = selected_peaks;

  if (changed)
  {
    select_in_table();
    select_in_plot();
  }
}

void FormEnergyCalibration::select_in_table()
{
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);
  ui->tablePeaks->clearSelection();

  QItemSelectionModel *selectionModel = ui->tablePeaks->selectionModel();
  QItemSelection itemSelection = selectionModel->selection();

  for (int i=0; i < ui->tablePeaks->rowCount(); ++i)
    if (selected_peaks_.count(ui->tablePeaks->item(i, 0)->data(Qt::UserRole).toDouble())) {
      ui->tablePeaks->selectRow(i);
      itemSelection.merge(selectionModel->selection(), QItemSelectionModel::Select);
    }

  selectionModel->clearSelection();
  selectionModel->select(itemSelection,QItemSelectionModel::Select);

  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);
}

void FormEnergyCalibration::replot_calib() {
  ui->PlotCalib->clearAll();

  QVector<double> xx, yy;

  double xmin = std::numeric_limits<double>::max();
  double xmax = - std::numeric_limits<double>::max();

  for (auto &q : fit_data_.peaks()) {
    double x = q.first;
    double y = q.second.energy();

    xx.push_back(x);
    yy.push_back(y);

    if (x < xmin)
      xmin = x;
    if (x > xmax)
      xmax = x;
  }

  double x_margin = (xmax - xmin) / 10;
  xmax += x_margin;
  xmin -= x_margin;

  if (xx.size()) {
    QVector<double> xx_sigma(yy.size(), 0);
    QVector<double> yy_sigma(yy.size(), 0);

    ui->PlotCalib->addPoints(style_pts, xx, yy, xx_sigma, yy_sigma);
    ui->PlotCalib->set_selected_pts(selected_peaks_);
    if (new_calibration_.valid())
    {
      double step = (xmax-xmin) / 50.0;
      xx.clear(); yy.clear();

      for (double i=xmin; i < xmax; i+=step)
      {
        xx.push_back(i);
        // \todo bit shift first
        yy.push_back(new_calibration_.transform(i));
      }
      ui->PlotCalib->setFit(xx, yy, style_fit);
      ui->PlotCalib->setTitle("E = " + QString::fromStdString(new_calibration_.function()->to_UTF8(6, true)));
    }
  }

  ui->PlotCalib->replotAll();
}

void FormEnergyCalibration::rebuild_table() {
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);

  std::set<double> flagged;
  ui->tablePeaks->clearContents();
  ui->tablePeaks->setRowCount(fit_data_.peaks().size());
  int i=0;
  for (auto &p : fit_data_.peaks()) {
    bool close = false;
    for (auto &e : ui->isotopes->current_isotope_gammas())
      if (std::abs(p.second.energy() - e.energy) < 2.0) {//hardcoded
        flagged.insert(e.energy);
        close = true;
      }
    add_peak_to_table(p.second, i, close);
    ++i;
  }

//  ui->isotopes->select_energies(flagged);

  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);
}

void FormEnergyCalibration::selection_changed_in_plot()
{
  selected_peaks_ = ui->PlotCalib->get_selected_pts();
  select_in_table();
  if (isVisible())
    emit selection_changed(selected_peaks_);
  toggle_push();
}

void FormEnergyCalibration::selection_changed_in_table() {
  selected_peaks_.clear();
  foreach (QModelIndex i, ui->tablePeaks->selectionModel()->selectedRows())
    selected_peaks_.insert(ui->tablePeaks->item(i.row(), 0)->data(Qt::UserRole).toDouble());

  select_in_plot();
  if (isVisible())
    emit selection_changed(selected_peaks_);
  toggle_push();
}

void FormEnergyCalibration::toggle_push() {
  size_t sel = selected_peaks_.size();

  ui->pushEnergiesToPeaks->setEnabled((sel > 0) && (sel == ui->isotopes->current_gammas().size()));
  ui->pushPeaksToNuclide->setEnabled((sel > 0) && (ui->isotopes->current_gammas().empty()));

  if (static_cast<int>(fit_data_.peaks().size()) > 1) {
    ui->pushFit->setEnabled(true);
    ui->spinTerms->setEnabled(true);
  } else {
    ui->pushFit->setEnabled(false);
    ui->spinTerms->setEnabled(false);
  }

  if (detectors_.get(fit_data_.detector_).get_calibration({"energy", fit_data_.detector_.id()}, {"energy"}).valid())
    ui->pushFromDB->setEnabled(true);
  else
    ui->pushFromDB->setEnabled(false);

  ui->pushApplyCalib->setEnabled(fit_data_.settings().cali_nrg_ != new_calibration_);
}

void FormEnergyCalibration::on_pushFit_clicked()
{
  auto optimizer = DAQuiri::OptimizerFactory::getInstance().create_any();
  if (!optimizer)
    return;

  std::vector<double> x, y;
  x.resize(fit_data_.peaks().size());
  y.resize(fit_data_.peaks().size());
  int i = 0;
  for (auto &q : fit_data_.peaks())
  {
    x[i] = q.first;
    y[i] = q.second.energy();
    i++;
  }

//  std::vector<double> sigmas(y.size(), 1);

  auto p = DAQuiri::CoefFunctionFactory::singleton().create_type("Polynomial");
  p->set_coeff(0, {-50, 50, 0});
  p->set_coeff(1,   {0, 50, 0.5});
  for (int i=2; i <= ui->spinTerms->value(); ++i)
    p->set_coeff(i, {-5, 5, 0.5});

  optimizer->fit(p, x, y, std::vector<double>(), std::vector<double>());

  if (p->coeffs().size())
  {
    new_calibration_ = DAQuiri::Calibration({"energy", fit_data_.detector_.id()},{"energy", fit_data_.detector_.id(), "keV"});
    new_calibration_.function(p);
  }
  else
    WARN("<Energy calibration> DAQuiri::Calibration failed");

  replot_calib();
  select_in_plot();
  toggle_push();
  emit new_fit();
}

void FormEnergyCalibration::isotope_energies_chosen() {
  update_data();
}

void FormEnergyCalibration::on_pushApplyCalib_clicked()
{
  emit update_detector();
}

void FormEnergyCalibration::on_pushFromDB_clicked()
{
  new_calibration_ = detectors_.get(fit_data_.detector_).get_calibration({"energy", fit_data_.detector_.id()}, {"energy"});
  replot_calib();
  select_in_plot();
  toggle_push();
  emit new_fit();
}

void FormEnergyCalibration::on_pushDetDB_clicked()
{
  // \todo reenable this
//  WidgetDetectors *det_widget = new WidgetDetectors(this);
//  det_widget->setData(detectors_);
//  connect(det_widget, SIGNAL(detectorsUpdated()), this, SLOT(detectorsUpdated()));
//  det_widget->exec();
}

void FormEnergyCalibration::on_pushPeaksToNuclide_clicked()
{
  std::vector<double> gammas;
  for (auto &q : fit_data_.peaks())
    if (selected_peaks_.count(q.first))
      gammas.push_back(q.second.energy());
  ui->isotopes->push_energies(gammas);
}

void FormEnergyCalibration::on_pushEnergiesToPeaks_clicked()
{
  std::vector<double> gammas = ui->isotopes->current_gammas();
  std::sort(gammas.begin(), gammas.end());

  std::vector<double> peakIDs;
  double last_sel = -1;
  for (auto &q : fit_data_.peaks())
    if (selected_peaks_.count(q.first)) {
      peakIDs.push_back(q.first);
      last_sel = q.first;
    }

  if (gammas.size() != peakIDs.size())
    return;

  for (size_t i=0; i<gammas.size(); ++i)
    fit_data_.override_energy(peakIDs.at(i), gammas.at(i));

  selected_peaks_.clear();
  for (auto &q : fit_data_.peaks())
    if (q.first > last_sel) {
      selected_peaks_.insert(q.first);
      break;
    }

  ui->isotopes->select_next_energy();

  update_data();
  emit change_peaks();
  emit selection_changed(selected_peaks_);
}

void FormEnergyCalibration::add_peak_to_table(const DAQuiri::Peak &p, int row, bool gray) {
  QBrush background(gray ? Qt::lightGray : Qt::white);

  add_to_table(ui->tablePeaks, row, 0, QString::number(p.center()),
               QVariant::fromValue(p.center()), background);
  // \todo uncertainty
  //add_to_table(ui->tablePeaks, row, 1, p.center().error_percent(), QVariant(), background);
  add_to_table(ui->tablePeaks, row, 2, QString::number(p.energy()), QVariant(), background);

}

void FormEnergyCalibration::select_in_plot()
{
  ui->PlotCalib->set_selected_pts(selected_peaks_);
  ui->PlotCalib->replotAll();
}
