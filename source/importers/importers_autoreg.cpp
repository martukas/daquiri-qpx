#include <core/importer_factory.h>
#include <importers/importers_autoreg.h>
#include <importers/ImporterAVA.h>
#include <importers/ImporterCHN.h>
#include <importers/ImporterCNF.h>
#include <importers/ImporterMCA.h>
#include <importers/ImporterN42.h>
#include <importers/ImporterSPC.h>
#include <importers/ImporterTKA.h>

using namespace DAQuiri;

void importers_autoreg()
{
  DAQUIRI_REGISTER_IMPORTER(ImporterAVA)
  DAQUIRI_REGISTER_IMPORTER(ImporterCHN)
  DAQUIRI_REGISTER_IMPORTER(ImporterCNF)
  DAQUIRI_REGISTER_IMPORTER(ImporterMCA)
  DAQUIRI_REGISTER_IMPORTER(ImporterN42)
  DAQUIRI_REGISTER_IMPORTER(ImporterSPC)
  DAQUIRI_REGISTER_IMPORTER(ImporterTKA)
}
