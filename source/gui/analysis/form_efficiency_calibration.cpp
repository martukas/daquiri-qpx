#include <gui/analysis/form_efficiency_calibration.h>
#include "ui_form_efficiency_calibration.h"

#include <core/consumer_factory.h>
//#include "widget_detectors.h"
#include <core/fitting/fitter.h>
#include <QInputDialog>
#include <gui/widgets/qt_util.h>
#include <gui/widgets/QFileExtensions.h>
#include <gui/daq/ConsumerDialog.h>

#include <core/calibration/coef_function_factory.h>

using namespace DAQuiri;

FormEfficiencyCalibration::FormEfficiencyCalibration(QWidget* parent)
    : QWidget(parent)
      , ui(new Ui::FormEfficiencyCalibration)
      , current_spectrum_(0)
{
//  optimizer_ = DAQuiri::OptimizerFactory::singleton().create_any();

  ui->setupUi(this);
  this->setWindowTitle("Efficiency calib");

  loadSettings();

  style_fit.default_pen = QPen(Qt::blue, 0);
  style_pts.themes["selected"] = QPen(Qt::black, 7);

  ui->PlotCalib->setAxisLabels("channel", "energy");

  ui->tablePeaks->verticalHeader()->hide();
  ui->tablePeaks->setColumnCount(3);
  ui->tablePeaks->setHorizontalHeaderLabels({"energy", "cps", QString(QChar(0x03B5)) + "-rel"});
  ui->tablePeaks->setSelectionBehavior(QAbstractItemView::SelectRows);
  ui->tablePeaks->setSelectionMode(QAbstractItemView::ExtendedSelection);
  ui->tablePeaks->setEditTriggers(QAbstractItemView::NoEditTriggers);
  ui->tablePeaks->horizontalHeader()->setStretchLastSection(true);
  ui->tablePeaks->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  ui->tablePeaks->show();
  connect(ui->tablePeaks, SIGNAL(itemSelectionChanged()), this, SLOT(selection_changed_in_table()));

  connect(ui->PlotCalib, SIGNAL(selection_changed()), this, SLOT(selection_changed_in_calib_plot()));
  ui->PlotCalib->set_log_x(true);
  ui->PlotCalib->setScaleType("Logarithmic");

  ui->isotopes->show();
  connect(ui->isotopes, SIGNAL(isotopeSelected()), this, SLOT(isotope_chosen()));

  QShortcut* shortcut = new QShortcut(QKeySequence(QKeySequence::Delete), ui->tablePeaks);
  connect(shortcut, SIGNAL(activated()), this, SLOT(on_pushMarkerRemove_clicked()));

  ui->isotopes->set_editable(false);



  //file formats for opening mca spectra
  std::vector<std::string> spectypes = DAQuiri::ConsumerFactory::singleton().types();
  QStringList filetypes;
  for (auto& q : spectypes)
  {
    ConsumerMetadata type_template = DAQuiri::ConsumerFactory::singleton().create_prototype("1D");
//    if (!type_template.input_types().empty())
//      filetypes.push_back(
//          "Spectrum " + QString::fromStdString(q) + "(" + catExtensions(type_template.input_types()) + ")");
  }
  mca_load_formats_ = catFileTypes(filetypes);

  ui->plotSpectrum->setFit(&fit_data_);

  connect(ui->plotSpectrum,
          SIGNAL(selection_changed(std::set<double>)),
          this,
          SLOT(update_selection(std::set<double>)));
  connect(ui->plotSpectrum, SIGNAL(data_changed()), this, SLOT(update_data()));

  connect(ui->spectrumSelector, SIGNAL(itemSelected(SelectorItem)), this, SLOT(spectrumDetails(SelectorItem)));
  connect(ui->spectrumSelector, SIGNAL(itemToggled(SelectorItem)), this, SLOT(spectrumLooksChanged(SelectorItem)));
}

FormEfficiencyCalibration::~FormEfficiencyCalibration()
{
  delete ui;
}

void FormEfficiencyCalibration::closeEvent(QCloseEvent* event)
{
  //if (!ui->isotopes->save_close()) {
  //  event->ignore();
  //  return;
  //}
  saveSettings();
  event->accept();
}

void FormEfficiencyCalibration::loadSettings()
{
  QSettings settings;
  ui->plotSpectrum->loadSettings(settings);

  settings.beginGroup("Program");
//  ui->isotopes->setDir(settings.value("settingsdirectory", QDir::homePath() + "/DAQuiri/settings").toString());
  data_directory_ = settings.value("save_directory", QDir::currentPath()).toString();
  settings.endGroup();

  settings.beginGroup("Efficiency_calibration");
  ui->spinTerms->setValue(settings.value("fit_function_terms", 2).toInt());
  ui->isotopes->set_current_isotope(settings.value("current_isotope", "Co-60").toString());
  ui->doubleEpsilonE->setValue(settings.value("epsilon_e", 2.0).toDouble());
  settings.endGroup();
}

void FormEfficiencyCalibration::saveSettings()
{
  QSettings settings;
  ui->plotSpectrum->saveSettings(settings);
  settings.beginGroup("Efficiency_calibration");
  settings.setValue("fit_function_terms", ui->spinTerms->value());
  settings.setValue("current_isotope", ui->isotopes->current_isotope());
  settings.setValue("epsilon_e", ui->doubleEpsilonE->value());
  settings.endGroup();
}

void FormEfficiencyCalibration::setDetector(DAQuiri::ProjectPtr newset, QString detector)
{
//  my_peak_fitter_->clear();

  current_detector_ = detector.toStdString();
  ui->labelCalibForDet->setText("Efficiency calibration for " + QString::fromStdString(current_detector_));

  project_ = newset;

  update_spectra();
}

void FormEfficiencyCalibration::setSpectrum(int64_t idx)
{
  if (!fit_data_.empty())  //should be ==Setting()
    peak_sets_[current_spectrum_] = fit_data_;

  auto spectrum = project_->get_consumer(idx);

  if (spectrum)
  {
    fit_data_ = Fitter();
    current_spectrum_ = idx;

    if (peak_sets_.count(idx))
    {
      fit_data_ = peak_sets_.at(idx);
      ui->isotopes->set_current_isotope(QString::fromStdString(fit_data_.sample_name_));
    }
    else
    {
      ConsumerMetadata md = spectrum->metadata();
      Setting descr = md.get_attribute("description");
      if (!descr.get_text().empty())
      {
        //find among data
        ui->isotopes->set_current_isotope(QString::fromStdString(descr.get_text()));
      }

      fit_data_.clear();
      fit_data_.setData(spectrum);
      peak_sets_[idx] = fit_data_;
    }
    ui->doubleScaleFactor->setValue(fit_data_.activity_scale_factor_);
    ui->doubleScaleFactor->setEnabled(true);
  }
  else
  {
    current_spectrum_ = 0;
    fit_data_ = Fitter();
    ui->plotSpectrum->update_spectrum();
    ui->doubleScaleFactor->setEnabled(false);
  }

  ui->plotSpectrum->update_spectrum();

  update_data();
}

void FormEfficiencyCalibration::update_spectrum()
{
  if (this->isVisible())
    ui->plotSpectrum->update_spectrum();
}

void FormEfficiencyCalibration::update_data()
{
  auto calib = fit_data_.settings().calib;
  auto live_time = fit_data_.metadata_.get_attribute("live_time").duration().count() * 0.001;

  UncertainDouble max {0.0, 0.0};
  std::set<UncertainDouble> flagged;
  fit_data_.sample_name_ = ui->isotopes->current_isotope().toStdString();
  for (auto& q : fit_data_.peaks())
  {
//    q.second.flagged = false;
    for (auto& p : ui->isotopes->current_isotope_gammas())
    {
      double diff = std::abs(q.second.peak_energy(calib.cali_nrg_).value() - p.energy.value());
      if (diff < ui->doubleEpsilonE->value())
      {
        Peak pk = q.second; //BROKEN
//        pk.flagged = true;
        pk.theoretical_intensity = p.abundance;
        pk.relative_efficiency = (pk.rate(live_time) / pk.theoretical_intensity);
        if (pk.relative_efficiency > max)
          max = pk.relative_efficiency;
        flagged.insert(pk.peak_energy(calib.cali_nrg_));
      }
    }
  }

  ui->isotopes->select_energies(flagged);

  //BROKEN
  if (max.value() > 0)
  {
    for (auto& q : fit_data_.peaks())
    {
//      if (q.second.theoretical_intensity > 0)
//        q.second.relative_efficiency = q.second.relative_efficiency / max;
    }
  }

  if (!fit_data_.empty()) //should be ==Setting()
    peak_sets_[current_spectrum_] = fit_data_;

//  ui->plotSpectrum->update_fit(contents_changed);
  replot_calib();
  rebuild_table(true);
  toggle_push();
}

void FormEfficiencyCalibration::update_detector_calibs()
{
//  std::string msg_text("Propagating calibration ");
//  msg_text += "<nobr>" + new_calibration_.debug() + "</nobr><br/>";
//
//  std::string question_text("Do you want to save this calibration to ");
//  question_text += current_detector_ + " in detector database?";
//
//  QMessageBox msgBox;
//  msgBox.setText(QString::fromStdString(msg_text));
//  msgBox.setInformativeText(QString::fromStdString(question_text));
//  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
//  msgBox.setDefaultButton(QMessageBox::No);
//  msgBox.setIcon(QMessageBox::Question);
//  int ret = msgBox.exec();
//
//  Detector modified;
//
//  if (ret == QMessageBox::Yes)
//  {
//    if (!detectors_.has_a(Detector(current_detector_)))
//    {
//      bool ok;
//      QString text = QInputDialog::getText(this, "New Detector",
//                                           "Detector name:", QLineEdit::Normal,
//                                           QString::fromStdString(current_detector_),
//                                           &ok);
//
//      if (!ok)
//        return;
//
//      if (!text.isEmpty())
//      {
//        modified = Detector(current_detector_);
//        if (detectors_.has_a(modified))
//        {
//          QMessageBox::warning(this,
//                               "Already exists",
//                               "Detector " + text + " already exists. Will not save to database.",
//                               QMessageBox::Ok);
//          modified = Detector();
//        }
//      }
//    }
//    else
//      modified = detectors_.get(Detector(current_detector_));
//
//    if (modified != Detector())
//    {
//      LINFO << "   applying new calibrations for " << modified.name()
//            << " in detector database";
//      modified.set_efficiency_calibration(new_calibration_);
//      detectors_.replace(modified);
//      emit detectorsChanged();
//    }
//  }
}

void FormEfficiencyCalibration::update_spectra()
{
  QVector<SelectorItem> items;

  for (auto& q : project_->get_consumers())
  {
    ConsumerMetadata md;
    if (q)
      md = q->metadata();

    if (!md.detectors.empty() && (md.detectors.front().id() == current_detector_))
    {
      SelectorItem new_spectrum;
      new_spectrum.text = QString::fromStdString(md.get_attribute("name").get_text());
      new_spectrum.color = QColor(QString::fromStdString(md.get_attribute("appearance").get_text()));
      new_spectrum.visible = md.get_attribute("visible").get_bool();
      items.push_back(new_spectrum);
    }
  }

  ui->spectrumSelector->setItems(items);

  setSpectrum(ui->spectrumSelector->selected().data.toLongLong());
}

void FormEfficiencyCalibration::spectrumLooksChanged(SelectorItem /*item*/)
{
  replot_calib();
}

void FormEfficiencyCalibration::spectrumDetails(SelectorItem item)
{
  setSpectrum(item.data.toLongLong());
}

void FormEfficiencyCalibration::add_peak_to_table(const Peak& p, int row, QColor bckg)
{
  QBrush background(bckg);

  auto calib = fit_data_.settings().calib;
  auto live_time = fit_data_.metadata_.get_attribute("live_time").duration().count() * 0.001;

  add_to_table(ui->tablePeaks, row, 0,
      QS(p.peak_energy(calib.cali_nrg_).to_string()),
      QVariant::fromValue(p.id()), background);
  add_to_table(ui->tablePeaks, row, 1, QS(p.rate(live_time).to_string()), {}, background);
  add_to_table(ui->tablePeaks, row, 2, QS(p.relative_efficiency.to_string(false)), {}, background);
}

void FormEfficiencyCalibration::replot_calib()
{
//  ui->PlotCalib->setFloatingText("");
  ui->PlotCalib->clearAll();

  QVector<double> xx, yy;

  double xmin = std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();

  QVector<SelectorItem> items = ui->spectrumSelector->items();

  bool have_data = false;

  for (auto& fit : peak_sets_)
  {
    bool visible = false;
    QColor color;
    for (auto& i : items)
    {
      if (i.visible && (fit.second.metadata_.get_attribute("name").get_text() == i.text.toStdString()))
      {
        visible = true;
        color = i.color;
      }
    }

    auto calib = fit.second.settings().calib;
    auto live_time = fit.second.metadata_.get_attribute("live_time").duration().count() * 0.001;

    if (visible)
    {
      have_data = true;
      xx.clear();
      yy.clear();

      for (auto& q : fit.second.peaks())
      {
        if (q.second.theoretical_intensity.value() == 0)
          continue;

        double x = q.second.peak_energy(calib.cali_nrg_).value();
        double y = q.second.relative_efficiency.value() * fit.second.activity_scale_factor_;

        xx.push_back(x);
        yy.push_back(y);

        if (x < xmin)
          xmin = x;
        if (x > xmax)
          xmax = x;
      }

      QVector<double> yy_sigma(yy.size(), 0);

      style_pts.default_pen = QPen(color, 7);
      ui->PlotCalib->addPoints(style_pts, xx, yy, QVector<double>(), yy_sigma);
    }
  }

  if (have_data)
  {
    ui->PlotCalib->set_selected_pts(selected_peaks_);
    if (new_calibration_.valid())
    {
      xx.clear();
      yy.clear();
      double step = (xmax - xmin) / 300.0;
      xmin -= step;
      if (xmin < 0)
        xmin = 0;

      xmax += (step * 100);

      for (double i = xmin; i < xmax; i += step)
      {
        double y = new_calibration_.transform(i);
        xx.push_back(i);
        yy.push_back(y);
      }
      ui->PlotCalib->setFit(xx, yy, style_fit);
      auto eqn = new_calibration_.function();
      QString eqn_str;
      if (eqn)
        eqn_str = QS(eqn->to_UTF8(3, true));
      ui->PlotCalib->setTitle("\u03B5 = " + eqn_str);
    }
  }

}

void FormEfficiencyCalibration::rebuild_table(bool contents_changed)
{
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);

  if (contents_changed)
  {
    ui->tablePeaks->clearContents();
    ui->tablePeaks->setRowCount(fit_data_.peaks().size());

    bool gray = false;
    QColor background_col;
    int i = 0;
    for (auto& q : fit_data_.peaks())
    {
/*      if (q.second.flagged && gray)
        background_col = Qt::darkGreen;
      else */if (gray)
        background_col = Qt::lightGray;
//      else if (q.second.flagged)
//        background_col = Qt::green;
      else
        background_col = Qt::white;
      add_peak_to_table(q.second, i, background_col);
      ++i;
//      if (!q.second.intersects_R)
//        gray = !gray;
    }
  }

  ui->tablePeaks->clearSelection();
  int i = 0;
  for (auto& q : fit_data_.peaks())
  {
    if (selected_peaks_.count(q.second.position.val()) > 0)
    {
      ui->tablePeaks->selectRow(i);
    }
    ++i;
  }
  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);
}

void FormEfficiencyCalibration::update_selection(std::set<double> selected_peaks)
{
  bool changed = (selected_peaks_ != selected_peaks);
  selected_peaks_ = selected_peaks;

  if (changed)
    select_in_table();
}

void FormEfficiencyCalibration::select_in_table()
{
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);
  ui->tablePeaks->clearSelection();
  int i = 0;
  for (auto& q : fit_data_.peaks())
  {
    if (selected_peaks_.count(q.second.position.val()) > 0)
    {
      ui->tablePeaks->selectRow(i);
    }
    ++i;
  }
  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);
}

void FormEfficiencyCalibration::selection_changed_in_calib_plot()
{
  selected_peaks_ = ui->PlotCalib->get_selected_pts();
  ui->plotSpectrum->set_selected_peaks(ui->PlotCalib->get_selected_pts());
  rebuild_table(false);
  toggle_push();
}

void FormEfficiencyCalibration::selection_changed_in_table()
{
  selected_peaks_.clear();
      foreach (QModelIndex i,
               ui->tablePeaks->selectionModel()->selectedRows())selected_peaks_.insert(ui->tablePeaks->item(i.row(),
                                                                                                            0)->data(Qt::EditRole).toDouble());
  ui->plotSpectrum->set_selected_peaks(selected_peaks_);
  replot_calib();
  toggle_push();
}

void FormEfficiencyCalibration::toggle_push()
{
  int unflagged = 0; //BROKEN

  ui->pushCullPeaks->setEnabled(unflagged > 0);

  int points_for_calib = 0;
  QVector<SelectorItem> items = ui->spectrumSelector->items();
  for (auto& fit : peak_sets_)
  {
    bool visible = false;

    for (auto& i : items)
      if (i.visible && (fit.second.metadata_.get_attribute("name").get_text() == i.text.toStdString()))
        visible = true;

    if (visible)
    { //BROKEN
//      for (auto &q : fit.second.peaks())
//        if (q.second.flagged)
//          points_for_calib++;
    }
  }

  ui->pushFit->setEnabled(points_for_calib > 1);
  ui->pushFit_2->setEnabled(points_for_calib > 1);
  ui->pushFitEffit->setEnabled(points_for_calib > 1);
  ui->spinTerms->setEnabled(points_for_calib > 1);

  if (new_calibration_ != Calibration())
    ui->pushApplyCalib->setEnabled(true);
  else
    ui->pushApplyCalib->setEnabled(false);
}

void FormEfficiencyCalibration::on_pushFit_clicked()
{
//  if (!optimizer_)
//    return;
//
//  QVector<SelectorItem> items = ui->spectrumSelector->items();
//  std::vector<double> xx, yy;
//
//  for (auto& fit : peak_sets_)
//  {
//    bool visible = false;
//
//    for (auto& i : items)
//      if (i.visible && (fit.second.metadata_.get_attribute("name").get_text() == i.text.toStdString()))
//        visible = true;
//
//    if (visible)
//    {
//      for (auto& q : fit.second.peaks())
//      { //bROKEN
////        if (!q.second.flagged) //BROKEN
////          continue;
//
//        double x = q.second.energy().value();
//        double y = q.second.relative_efficiency * fit.second.activity_scale_factor_;
//
//        xx.push_back(x);
//        yy.push_back(y);
//      }
//
//    }
//  }
//
////  PolyLog p = PolyLog(xx, yy, ui->spinTerms->value());
//
////  std::vector<double> sigmas(yy.size(), 1);
//
//  auto p = std::make_shared<PolyLog>();
//  p->add_coeff(0, -50, 50, 1);
//  for (int i = 1; i <= ui->spinTerms->value(); ++i)
//  {
//    if (i == 1)
//      p->add_coeff(i, -50, 50, 1);
//    else
//      p->add_coeff(i, -50, 50, 0);
//  }
//
//  optimizer_->fit(p, xx, yy, std::vector<double>(), std::vector<double>());
//
//  if (p->coeff_count())
//  {
////    new_calibration_.type_ = "Efficiency";
//    new_calibration_ = DAQuiri::Calibration(fit_data_.settings().bits_);
//    new_calibration_.set_units("ratio");
//    new_calibration_.set_function(p);
//    DBG << "<Efficiency calibration> new calibration fit " << new_calibration_.debug();
//  }
//  else
//    LINFO << "<Efficiency calibration> Calibration failed";
//
//  replot_calib();
//  toggle_push();
}

void FormEfficiencyCalibration::isotope_chosen()
{
  update_data();
}

void FormEfficiencyCalibration::on_pushApplyCalib_clicked()
{
  update_detector_calibs();
}

void FormEfficiencyCalibration::on_pushDetDB_clicked()
{
//  WidgetDetectors* det_widget = new WidgetDetectors(this);
//  det_widget->setData(detectors_);
//  connect(det_widget, SIGNAL(detectorsUpdated()), this, SLOT(detectorsUpdated()));
//  det_widget->exec();
}

void FormEfficiencyCalibration::on_pushCullPeaks_clicked()
{ //BROKEN
//  for (auto &q : fit_data_.peaks())
//    if (!q.second.flagged)
//      q.second.selected = true;
//  on_pushMarkerRemove_clicked();
}

void FormEfficiencyCalibration::on_doubleEpsilonE_editingFinished()
{
  update_data();
}

void FormEfficiencyCalibration::on_doubleScaleFactor_editingFinished()
{
  fit_data_.activity_scale_factor_ = ui->doubleScaleFactor->value();
  if (!fit_data_.empty()) //should be ==Setting()
    peak_sets_[current_spectrum_] = fit_data_;
  replot_calib();
}

void FormEfficiencyCalibration::on_doubleScaleFactor_valueChanged(double)
{
  fit_data_.activity_scale_factor_ = ui->doubleScaleFactor->value();
  if (!fit_data_.empty()) //should be ==Setting()
    peak_sets_[current_spectrum_] = fit_data_;
  replot_calib();
}

void FormEfficiencyCalibration::on_pushFit_2_clicked()
{
//  if (!optimizer_)
//    return;
//
//  QVector<SelectorItem> items = ui->spectrumSelector->items();
//  std::vector<double> xx, yy;
//
//  for (auto& fit : peak_sets_)
//  {
//    bool visible = false;
//
//    for (auto& i : items)
//      if (i.visible && (fit.second.metadata_.get_attribute("name").get_text() == i.text.toStdString()))
//        visible = true;
//
//    if (visible)
//    {
//      for (auto& q : fit.second.peaks())
//      {
//        //BROKEN
//
////        if (!q.second.flagged)
////          continue;
//
//        double x = q.second.energy().value();
//        double y = q.second.relative_efficiency * fit.second.activity_scale_factor_;
//
//        xx.push_back(x);
//        yy.push_back(y);
//      }
//    }
//  }
//
////  LogInverse p = LogInverse(xx, yy, ui->spinTerms->value());
//
////  std::vector<double> sigmas(yy.size(), 1);
//
//  auto p = std::make_shared<LogInverse>();
//  p->add_coeff(0, -50, 50, 1);
//  for (int i = 1; i <= ui->spinTerms->value(); ++i)
//  {
//    if (i == 1)
//      p->add_coeff(i, -50, 50, 1);
//    else
//      p->add_coeff(i, -50, 50, 0);
//  }
//
//  optimizer_->fit(p, xx, yy, std::vector<double>(), std::vector<double>());
//
//  if (p->coeff_count())
//  {
//    new_calibration_ = Calibration(fit_data_.settings().bits_);
//    new_calibration_.set_units("ratio");
//    new_calibration_.set_function(p);
//    DBG << "<Efficiency calibration> new calibration fit " << new_calibration_.debug();
//  }
//  else
//    LINFO << "<Efficiency calibration> Calibration failed";
//
//  replot_calib();
//  toggle_push();
}

void FormEfficiencyCalibration::on_pushFitEffit_clicked()
{
//  if (!optimizer_)
//    return;
//
//  QVector<SelectorItem> items = ui->spectrumSelector->items();
//  std::vector<double> xx, yy;
//
//  for (auto& fit : peak_sets_)
//  {
//    bool visible = false;
//
//    for (auto& i : items)
//      if (i.visible && (fit.second.metadata_.get_attribute("name").get_text() == i.text.toStdString()))
//        visible = true;
//
//    if (visible)
//    {
//      for (auto& q : fit.second.peaks())
//      {
//
//        // BROKEN
//        //        if (!q.second.flagged)
//        //          continue;
//
//        double x = q.second.energy().value();
//        double y = q.second.relative_efficiency * fit.second.activity_scale_factor_;
//
//        xx.push_back(x);
//        yy.push_back(y);
//      }
//
//    }
//  }
//
//  auto p = std::make_shared<Effit>();
//  optimizer_->fit(p, xx, yy, std::vector<double>(), std::vector<double>());
//
////  new_calibration_.type_ = "Efficiency";
//  new_calibration_ = Calibration(fit_data_.settings().bits_);
//  new_calibration_.set_units("ratio");
//  new_calibration_.set_function(p);
//  DBG << "<Efficiency calibration> new calibration fit " << new_calibration_.debug();
//
//  replot_calib();
//  toggle_push();
}
