#include "ImporterSPC.h"
#include <core/consumer_factory.h>
#include <date/date.h>
#include <iostream>
#include <core/calibration/polynomial.h>

#include <core/util/custom_logger.h>

struct __attribute__ ((packed)) Header
{
  int16_t INFTYP;        //Must be 1
  int16_t FILTYP;        //Must be 1 or 5
  int16_t ContentsFlag;  //Bit 0 = 1 for long filename, Bit 1 = 1 for ZDT spectrum and ROI in file
  int16_t Reserved;
  int16_t ACQIRP;        //Acquisition information record pointer
  int16_t SAMDRP;        //Sample description record pointer
  int16_t DETDRP;        //Detector description record pointer
  int16_t EBRDESC;       //EBR description record pointer
  int16_t ANARP1;        //First analysis parameters record pointer
  int16_t ANARP2;        //Second analysis parameters record pointer         [10]
  int16_t ANARP3;        //Third analysis parameters record pointer
  int16_t ANARP4;        //Fourth analysis parameters record pointer
  int16_t SRPDES;        //Absorption correction description record pointer
  int16_t IEQDESC;       //IEQ description record pointer
  int16_t GEODES;        //Geometry correction description record pointer
  int16_t MPCDESC;       //MPC description record pointer
  int16_t CALDES;        //Calibration description record pointer
  int16_t CALRP1;        //First calibration data record pointer
  int16_t CALRP2;        //Second calibration data record pointer
  int16_t EFFPRP;        //Efficiency pairs record pointer (first record)    [20]
  int16_t ROIRP1;        //Record number of the first of the two ROI

  int16_t Enprp;         //Energy pairs record pointer
  int16_t Nenpr;         //Number of energy pair records
  int16_t Reserved2;     //Reserved
  int16_t UnkPkDecDis;   //Disable deconvolution of unknown peaks
  int16_t ActivUnit;     //True = microcuries, false = becquerels
  int16_t PERPTR;        //Laboratory and operator name record pointer
  int16_t MAXRCS;        //Maximum record number ever used
  int16_t LSTREC;        //Maximum record number in use
  int16_t EFFPNM;        //Number of efficiency pairs records (See Word 20)  [30]

  int16_t SPCTRP;        //Spectrum record pointer (pointer to first record)
  int16_t SPCRCN;        //Number of records in the spectrum
  int16_t SPCCHN;        //Number of channels in spectrum
  int16_t ABSTCH;        //Physical start channel for data
  float ACQTIM;          //Date and time of acquisition start in DECDAY format
  double ACQTI8;         //Date and time as double precision DECDAY
  int16_t SEQNUM;        //Sequence number                                 [41]
  int16_t
      MCANU;         //MCA number as two ASCII characters (old) or Detector number as integer for systems with Connections
  int16_t
      SEGNUM;        //Segment number as two ASCII characters (old) or as integer value 1 for systems with Connections
  int16_t MCADVT;        //MCA device type
  int16_t CHNSRT;        //Start channel number
  float RLTMDT;          //Real Time in seconds
  float LVTMDT;          //Live Time in seconds

  int16_t MUCPtr;        //Pointer to MGA or U235 or CZTU records          [50]
  int16_t FramPtr;       //Pointer to FRAM records
  int16_t UNKNOWN;       //missing from documentation.
  int16_t TrifidPtr;     //Pointer to TRIFID records
  int16_t NaiPtr;        //Pointer to NaI records

  int16_t Reserved3[8];  //[55-62]

  float RRSFCT;          //Total random summing factor
};

struct __attribute__ ((packed)) AcqInfoRecord
{
  char DefSpeFileName[16];
  char SpeDate[12];
  char SpeTime[10];
  char LiveTime[10];
  char RealTime[10];

  char Reserved[34];

  char SampleCollStartDate[10];
  char SampleCollStartTime[8];
  char SampleCollStopDate[10];
  char SampleCollStopTime[8];
};

struct __attribute__ ((packed)) CalibrationRecord
{
  int16_t AFIT;    //Above knee efficiency calibration fit type
  int16_t BFIT;    //Below knee efficiency calibration fit type
  int16_t EFFPRS;  //Number of efficiency pairs
  int16_t NCH;     //Number of channels in spectrum
  float KNEE;      //Detector knee (keV)
  float ASIG;      //2-sigma uncertainty above knee
  float BSIG;      //2-sigma uncertainty below knee
  float EC_1;      //Energy vs. channel coefficient A
  float EC_2;      //Energy vs. channel coefficient B
  float EC_3;      //Energy vs. channel coefficient C
  float FC_1;      //FWHM vs. channel coefficient A
  float FC_2;      //FWHM vs. channel coefficient B
  float FC_3;      //FWHM vs. channel coefficient C
  float PE_1;      //Above knee efficiency vs. energy coefficient A or polynomial coefficient _1
  float PE_2;      //Above knee efficiency vs. energy coefficient B or polynomial coefficient _2
  float PE_3;      //Above knee efficiency vs. energy coefficient C or polynomial coefficient _3
  float SE_1;      //Below knee efficiency vs. energy coefficient A or polynomial coefficient _4
  float SE_2;      //Below knee efficiency vs. energy coefficient B or polynomial coefficient _5
  float SE_3;      //Below knee efficiency vs. energy coefficient C or polynomial coefficient _6
  int16_t FWHTYP;  //FWHM type
  bool PETYPE;     // True for p-type
  int16_t MAESTRO; //peak-search sensitivity
  int16_t ENGPRS;  //Number of energy pairs
  int16_t DETNUM;  //Detector number
  int16_t NBKNEE;  //Number of calibration points below knee
  float ENA2;      //Temp energy calibration
  float ENB2;      //Temp energy calibration
  float ENC2;      //Temp energy calibration
  float CALUNC;    //Calibration source uncertainty
  float CALDIF;    //Energy calibration difference
  float R_7;       //Polynomial coefficient 7
  float R_8;       //Polynomial coefficient 8
  float R_9;       //Polynomial coefficient 9
  float R_10;      //Polynomial coefficient 10
  float lcfe;      //Low channel FWHM error
  float hcfe;      //High channel FWHM error
  char lclfc;      //Low channel limit for calibrating
  bool STYPEFLAG;  //True = next record has TCC data
  //??? String   UNKNOWN[7] //???
};

struct __attribute__ ((packed)) OrtecHardwareRecord2
{
  int32_t ValidityFlag_1_32; //use bits to validate entries below
  int32_t ValidityFlag_33_64;
  int16_t INPUTPOL; //+1 positive, -1 negative
  float THERMISTOR; //in Ohms
  int16_t PUR_IN_USE; //0=off, 1=on
  int32_t PUR_WIDTH; //in ns
  int16_t OCTETE_CURRENT; //in nA
  int16_t OCTETE_VACUUM; //in mT
  char FIRMWARE_REVISION[16]; //16 characters
  int16_t HV_ENABLED; //1=Yes, 0=No
  int16_t SHAPING_TIME_INDEX;
  int16_t BATTERY_STATUS; //0=External, 1=Battery1, 2=Battery2
  float BATTERY1_VOLTAGE;
  float BATTERY2_VOLTAGE;
  int16_t DESIRED_HV;
  int16_t HV_SHUTDOWN_MODE; //0=TTL, 1=ORTEC, 2=OFF
  int16_t ADC_TYPE; //0=CI34, 1=CI36, 2=ORTEC, 3=SILENA
  int16_t MAX_VACUUM;
  int16_t AUTO_THRESHOLD_FLAG; //0=OFF, 1= ON (SBS-60)
  float MDA_PRESET_1;
  float MDA_PRESET_2;
  float MDA_PRESET_3;
  float MDA_PRESET_VALUE;
  float MDA_PRESET_USER;
  int16_t LOW_MDA_PRESET_ROI;
  int16_t HIGH_MDA_PRESET_ROI;
  int16_t MDA_UNIT;  //0=Bq, 1=uCi, -1=NONE
  char MDA_NUCL_NAME[8];
  int16_t ZDT_ENABLED_FLAG; //1=On, 0=Off
  int16_t ZDT_REFRESH_RATE; //in ms
  int16_t ZDT_VIEW; //1=normal, 2=corrected
  float TOTAL_COUNT_RATE; //in cps
  int16_t MINI_MCA166_ANAL_THRSH;
  int16_t MINI_MCA166_ANAL_ROUTING;
  int16_t MINI_MCA166_ANAL_FW_REV;
  int16_t MINI_MCA166_ANAL_HW_REV;
  int16_t MINI_MCA166_ANAL_LOW_PZ;
  int16_t MINI_MCA166_ANAL_HIGH_PZ;
  int16_t MINI_MCA166_ANAL_SLOW_DISCR;
  int16_t MINI_MCA166_ANAL_FAST_DISCR;
  int16_t MINI_MCA166_ANAL_POWER_SW;
  int16_t ZDT_MODE; //0=OFF, 1=NORM_CORR, 2=CORR_ERR
};

static constexpr int16_t ContentsFlagMaskHasZDT = 2;

bool ImporterSPC::validate(const boost::filesystem::path& path) const
{
  (void) path;
  return true;
}

void ImporterSPC::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  std::ifstream file(path.string(), std::ios::binary);

  Header header;
  file.read(reinterpret_cast<char*>(&header), sizeof(Header));

  if ((header.INFTYP != 1) || ((header.FILTYP != 1) && (header.FILTYP != 5)))
    throw std::runtime_error("<ImporterSPC> invalid SPC file, bad type markers");

  INFO("contents flag: {}", header.ContentsFlag);
  bool ZDT_included = ((header.ContentsFlag & ContentsFlagMaskHasZDT) == ContentsFlagMaskHasZDT);

  file.seekg((header.ACQIRP - 1) * 128);
  AcqInfoRecord acq_header;
  file.read(reinterpret_cast<char*>(&acq_header), sizeof(AcqInfoRecord));

  file.seekg((header.ANARP4 - 1) * 128 + 128);
  OrtecHardwareRecord2 hw2header;
  file.read(reinterpret_cast<char*>(&hw2header), sizeof(OrtecHardwareRecord2));

  std::vector<int32_t> spectrum(static_cast<size_t>(header.SPCCHN));
  std::vector<int32_t> spectrum2(static_cast<size_t>(header.SPCCHN));

  file.seekg((header.SPCTRP - 1) * 128);
  file.read(reinterpret_cast<char*>(spectrum.data()), spectrum.size() * sizeof(int32_t));
  if (ZDT_included)
  {
    file.seekg((header.SPCTRP - 1) * 128 + spectrum.size() * sizeof(int32_t));
    file.read(reinterpret_cast<char*>(spectrum2.data()), spectrum2.size() * sizeof(int32_t));
  }

  file.seekg((header.CALRP1 - 1) * 128);
  CalibrationRecord calheader;
  file.read(reinterpret_cast<char*>(&calheader), sizeof(CalibrationRecord));

//  INFO("mca number: {}", header.MCANU);
//  INFO("segment number: {}", header.SEGNUM);
//  INFO("tag numbrr: {}", header.MCADVT);
//  INFO("zdt_mode: {}", hw2header.ZDT_MODE);
//
//  INFO("fc_cal c0: {}", calheader.FC_1);
//  INFO("fc_cal c1: {}", calheader.FC_2);
//  INFO("fc_cal c2: {}", calheader.FC_3);

  DAQuiri::CalibID from("energy","unknown","");
  DAQuiri::CalibID to("energy","unknown","keV");
  DAQuiri::Calibration new_calib(from, to);
  std::vector<double> coefs {calheader.EC_1, calheader.EC_2, calheader.EC_3};
  new_calib.function(std::make_shared<DAQuiri::Polynomial>(coefs));

  DAQuiri::Detector det;
  det.set_calibration(new_calib);
//  INFO("det = {}", det.debug(""));

  std::stringstream stream;
  std::string start_day(&acq_header.SpeDate[0], 2);
  std::string start_month(&acq_header.SpeDate[3], 3);
  std::string start_year =
      ((acq_header.SpeDate[9] == '1') ? "20" : "19")
          + std::string(&acq_header.SpeDate[7], 2);
  std::string start_time(&acq_header.SpeTime[0], 10);
  stream << start_year << "-" << start_month << "-" << start_day
         << " " << start_time;

  date::sys_time<std::chrono::seconds> t;
  stream >> date::parse("%Y-%b-%d %H:%M:%S", t);
  if (stream.fail())
    throw std::runtime_error("<ImporterSPC> failed to parse date");

  auto lt_ms = static_cast<int64_t>(header.LVTMDT * 1000.0);
  auto rt_ms = static_cast<int64_t>(header.RLTMDT * 1000.0);

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterSPC could not get a valid Histogram 1D from factory");
  hist->set_attribute(DAQuiri::Setting::text("name", path.stem().string()));
  hist->set_attribute(DAQuiri::Setting::boolean("visible", true));
  hist->set_attribute(DAQuiri::Setting("start_time", t));
  hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));
  hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));
  for (size_t i = 0; i < spectrum.size(); ++i)
  {
    DAQuiri::Entry new_entry;
    new_entry.first.resize(1);
    new_entry.first[0] = i;
    new_entry.second = PreciseFloat(spectrum[i]);
    entry_list.push_back(new_entry);
  }
  hist->set_detectors({det});
  hist->set_attribute(DAQuiri::Setting::text("value_latch", "energy"));
  hist->import(*this);
  project->add_consumer(hist);

  if (ZDT_included) {
    auto hist2 = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
    if (!hist2)
      throw std::runtime_error("ImporterSPC could not get a valid Histogram 1D from factory");
    hist2->set_attribute(DAQuiri::Setting::text("name", path.stem().string() + " (ZDT)"));
    hist2->set_attribute(DAQuiri::Setting::boolean("visible", true));
    hist2->set_attribute(DAQuiri::Setting("start_time", t));
    hist2->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(rt_ms)));
    hist2->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));

    entry_list.clear();
    for (size_t i = 0; i < spectrum2.size(); ++i)
    {
      DAQuiri::Entry new_entry;
      new_entry.first.resize(1);
      new_entry.first[0] = i;
      new_entry.second = PreciseFloat(spectrum2[i]);
      entry_list.push_back(new_entry);
    }
    hist2->set_detectors({det});
    hist2->set_attribute(DAQuiri::Setting::text("value_latch", "energy"));
    hist2->import(*this);
    project->add_consumer(hist2);
  }
}
