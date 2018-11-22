#include <core/importer_factory.h>

#include "importers_autoreg.h"
#include "ImporterAVA.h"
#include "ImporterCNF.h"
#include "ImporterMCA.h"
#include "ImporterN42.h"
#include "ImporterTKA.h"


using namespace DAQuiri;

void importers_autoreg()
{
  DAQUIRI_REGISTER_IMPORTER(ImporterAVA)
  DAQUIRI_REGISTER_IMPORTER(ImporterCNF)
  DAQUIRI_REGISTER_IMPORTER(ImporterMCA)
  DAQUIRI_REGISTER_IMPORTER(ImporterN42)
  DAQUIRI_REGISTER_IMPORTER(ImporterTKA)
}
