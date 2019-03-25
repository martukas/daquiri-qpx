#include <gui/analysis/form_analysis_1d.h>
//#include "widget_detectors.h"
#include "ui_form_analysis_1d.h"
#include <core/fitting/fitter.h>
#include <QInputDialog>
#include <QSettings>

#include <gui/widgets/QFileExtensions.h>

using namespace DAQuiri;

FormAnalysis1D::FormAnalysis1D(QWidget* parent) :
    QWidget(parent), ui(new Ui::FormAnalysis1D)
{
  ui->setupUi(this);

  loadSettings();

  ui->plotSpectrum->setFit(&fit_data_);

  //connect(ui->widgetDetectors, SIGNAL(detectorsUpdated()), this, SLOT(detectorsUpdated()));

  form_energy_calibration_ = new FormEnergyCalibration(detector_, fit_data_);
  ui->tabs->addTab(form_energy_calibration_, "Energy calibration");
  connect(form_energy_calibration_, SIGNAL(detectorsChanged()), this, SLOT(detectorsUpdated()));
  connect(form_energy_calibration_, SIGNAL(update_detector()), this, SLOT(update_detector_calibs()));

  form_fwhm_calibration_ = new FormFwhmCalibration(detector_, fit_data_);
  ui->tabs->addTab(form_fwhm_calibration_, "FWHM calibration");
  connect(form_fwhm_calibration_, SIGNAL(detectorsChanged()), this, SLOT(detectorsUpdated()));
  connect(form_fwhm_calibration_, SIGNAL(update_detector()), this, SLOT(update_detector_calibs()));

  form_fit_results_ = new FormFitResults(fit_data_);
  ui->tabs->addTab(form_fit_results_, "Peak integration");
  connect(form_fit_results_, SIGNAL(save_peaks_request()), this, SLOT(save_report()));

  connect(form_energy_calibration_,
          SIGNAL(selection_changed(std::set<double>)),
          form_fwhm_calibration_,
          SLOT(update_selection(std::set<double>)));
  connect(form_energy_calibration_,
          SIGNAL(selection_changed(std::set<double>)),
          form_fit_results_,
          SLOT(update_selection(std::set<double>)));
  connect(form_energy_calibration_,
          SIGNAL(selection_changed(std::set<double>)),
          ui->plotSpectrum,
          SLOT(set_selected_peaks(std::set<double>)));

  connect(form_fwhm_calibration_,
          SIGNAL(selection_changed(std::set<double>)),
          form_energy_calibration_,
          SLOT(update_selection(std::set<double>)));
  connect(form_fwhm_calibration_,
          SIGNAL(selection_changed(std::set<double>)),
          form_fit_results_,
          SLOT(update_selection(std::set<double>)));
  connect(form_fwhm_calibration_,
          SIGNAL(selection_changed(std::set<double>)),
          ui->plotSpectrum,
          SLOT(set_selected_peaks(std::set<double>)));

  connect(form_fit_results_,
          SIGNAL(selection_changed(std::set<double>)),
          form_fwhm_calibration_,
          SLOT(update_selection(std::set<double>)));
  connect(form_fit_results_,
          SIGNAL(selection_changed(std::set<double>)),
          form_energy_calibration_,
          SLOT(update_selection(std::set<double>)));
  connect(form_fit_results_,
          SIGNAL(selection_changed(std::set<double>)),
          ui->plotSpectrum,
          SLOT(set_selected_peaks(std::set<double>)));

  connect(ui->plotSpectrum,
          SIGNAL(peak_selection_changed(std::set<double>)),
          form_fwhm_calibration_,
          SLOT(update_selection(std::set<double>)));
  connect(ui->plotSpectrum,
          SIGNAL(peak_selection_changed(std::set<double>)),
          form_energy_calibration_,
          SLOT(update_selection(std::set<double>)));
  connect(ui->plotSpectrum,
          SIGNAL(peak_selection_changed(std::set<double>)),
          form_fit_results_,
          SLOT(update_selection(std::set<double>)));

  connect(ui->plotSpectrum, SIGNAL(data_changed()), form_energy_calibration_, SLOT(update_data()));
  connect(ui->plotSpectrum, SIGNAL(data_changed()), form_fwhm_calibration_, SLOT(update_data()));
  connect(ui->plotSpectrum, SIGNAL(data_changed()), form_fit_results_, SLOT(update_data()));
  connect(ui->plotSpectrum, SIGNAL(fitting_done()), this, SLOT(update_fit()));

  connect(form_energy_calibration_, SIGNAL(change_peaks()), form_fwhm_calibration_, SLOT(update_data()));
  connect(form_energy_calibration_, SIGNAL(change_peaks()), form_fit_results_, SLOT(update_data()));
  connect(form_energy_calibration_, SIGNAL(change_peaks()), ui->plotSpectrum, SLOT(updateData()));

  connect(form_energy_calibration_, SIGNAL(new_fit()), this, SLOT(update_fits()));
  connect(form_fwhm_calibration_, SIGNAL(new_fit()), this, SLOT(update_fits()));

  ui->tabs->setCurrentWidget(form_energy_calibration_);
}

FormAnalysis1D::~FormAnalysis1D()
{
  delete ui;
}

void FormAnalysis1D::closeEvent(QCloseEvent* event)
{
  if (!form_energy_calibration_->save_close())
  {
    event->ignore();
    return;
  }

  if (!form_fwhm_calibration_->save_close())
  {
    event->ignore();
    return;
  }

  QSettings settings_;
  ui->plotSpectrum->saveSettings(settings_);

//  DBG << "closing analysis1d";
  saveSettings();
  event->accept();
}

void FormAnalysis1D::loadSettings()
{
  QSettings settings_;

  settings_.beginGroup("Program");
  data_directory_ = settings_.value("save_directory", QDir::homePath() + "/qpx/data").toString();
  settings_.endGroup();

  ui->plotSpectrum->loadSettings(settings_);
}

void FormAnalysis1D::saveSettings()
{
  QSettings settings_;

//  settings_.beginGroup("Program");
//  settings_.setValue("save_directory", data_directory_);
//  settings_.endGroup();

  ui->plotSpectrum->saveSettings(settings_);
}

void FormAnalysis1D::clear()
{
  current_spectrum_ = 0;
  fit_data_.clear();
  ui->plotSpectrum->update_spectrum();
  form_energy_calibration_->clear();
  form_fwhm_calibration_->clear();
  form_fit_results_->clear();
  new_energy_calibration_ = Calibration();
  new_fwhm_calibration_ = Calibration();
}

void FormAnalysis1D::setSpectrum(ProjectPtr newset, int64_t idx)
{
  form_energy_calibration_->clear();
  form_fwhm_calibration_->clear();
  form_fit_results_->clear();
  spectra_ = newset;
  new_energy_calibration_ = Calibration();
  new_fwhm_calibration_ = Calibration();

  if (!spectra_)
  {
    fit_data_.clear();
    ui->plotSpectrum->update_spectrum();
    current_spectrum_ = 0;
    return;
  }

  current_spectrum_ = idx;
  ConsumerPtr spectrum = spectra_->get_consumer(idx);

  if (spectrum)
  {
    fit_data_.clear();
    fit_data_.setData(spectrum);

    // \todo reenable saving of analysis
//    if (spectra_->has_fitter(idx))
//      fit_data_ = spectra_->get_fitter(idx);
//    else
    fit_data_.setData(spectrum);

    form_energy_calibration_->newSpectrum();
    form_fwhm_calibration_->newSpectrum();

    form_energy_calibration_->update_data();
    form_fwhm_calibration_->update_data();
    form_fit_results_->update_data();

    new_energy_calibration_ = form_energy_calibration_->get_new_calibration();
    new_fwhm_calibration_ = form_fwhm_calibration_->get_new_calibration();
  }

  ui->plotSpectrum->update_spectrum();
}

void FormAnalysis1D::update_spectrum()
{
  if (this->isVisible())
  {
    ConsumerPtr spectrum = spectra_->get_consumer(current_spectrum_);
    if (spectrum)
      fit_data_.setData(spectrum);
    ui->plotSpectrum->update_spectrum();
  }
}

void FormAnalysis1D::update_fits()
{
  DBG("calib fits updated locally");
  new_energy_calibration_ = form_energy_calibration_->get_new_calibration();
  new_fwhm_calibration_ = form_fwhm_calibration_->get_new_calibration();
}

void FormAnalysis1D::update_fit()
{
  if (!spectra_)
    return;
  if (!spectra_->get_consumer(current_spectrum_))
    return;
  // \todo reenable saving of analysis
  //spectra_->update_fitter(current_spectrum_, fit_data_);
  // \todo emit something...
}

void FormAnalysis1D::update_detector_calibs()
{
  std::string msg_text("Propagating calibrations ");
  msg_text += "<nobr>Energy:" + new_energy_calibration_.debug() + "</nobr><br/>"
                                                                  "<nobr>FWHM:" + new_fwhm_calibration_.debug()
      + "</nobr><br/>"
        "<nobr>  to all spectra in originating project </nobr><br/>";

  std::string question_text("Do you also want to save this calibration to ");
  question_text += fit_data_.detector_.id() + " in detector database?";

  QMessageBox msgBox;
  msgBox.setText(QString::fromStdString(msg_text));
  msgBox.setInformativeText(QString::fromStdString(question_text));
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::Abort);
  msgBox.setDefaultButton(QMessageBox::No);
  msgBox.setIcon(QMessageBox::Question);
  int ret = msgBox.exec();

  Detector modified;

  // \todo reenabled this
//  if (ret == QMessageBox::Yes) {
//    if (!detector_.has_a(fit_data_.detector_)) {
//      bool ok;
//      QString text = QInputDialog::getText(this, "New Detector",
//                                           "Detector name:", QLineEdit::Normal,
//                                           QString::fromStdString(fit_data_.detector_.id()),
//                                           &ok);
//
//      if (!ok)
//        return;
//
//      if (!text.isEmpty()) {
//        modified = fit_data_.detector_;
//        modified.set_id(text.toStdString());
//        if (detector_.has_a(modified)) {
//          QMessageBox::warning(this, "Already exists", "Detector " + text + " already exists. Will not save to database.", QMessageBox::Ok);
//          modified = Detector();
//        }
//      }
//    } else
//      modified = detector_.get(fit_data_.detector_);
//
//    if (modified != Detector())
//    {
//      INFO("   applying new calibrations for {} in detector database", modified.id());
//      modified.set_calibration(new_energy_calibration_);
//      modified.set_calibration(new_fwhm_calibration_);
//      detector_.replace(modified);
//      emit detectorsChanged();
//    }
//  }

  if (ret != QMessageBox::Abort)
  {
    for (auto& q : spectra_->get_consumers())
    {
      ConsumerMetadata md = q->metadata();
      for (auto& p : md.detectors)
      {
        if (p.shallow_equals(fit_data_.detector_))
        {
          INFO("   applying new calibrations for {} on {}",
               fit_data_.detector_.id(),
               q->metadata().get_attribute("name").get_text());
          p.set_calibration(new_energy_calibration_);
          p.set_calibration(new_fwhm_calibration_);
        }
      }
      q->set_detectors(md.detectors);
    }

//    std::set<Spill> spills = spectra_->spills();
//    Spill sp;
//    if (spills.size())
//      sp = *spills.rbegin();
//    for (auto &p : sp.detectors) {
//      if (p.shallow_equals(fit_data_.detector_)) {
//        LINFO << "   applying new calibrations for " << fit_data_.detector_.name_ << " in current project " << spectra_->identity();
//        p.energy_calibrations_.replace(new_energy_calibration_);
//        p.fwhm_calibration_ = new_fwhm_calibration_;
//      }
//    }
//    spectra_->add_spill(&sp);
    update_spectrum();
  }
}

void FormAnalysis1D::save_report()
{
  QString fileName = CustomSaveFileDialog(this, "Save analysis report",
                                          data_directory_, "Plain text (*.txt)");
  if (validateFile(this, fileName, true))
  {
    INFO("Writing report to {}", fileName.toStdString());
    fit_data_.save_report(fileName.toStdString());
  }
}

