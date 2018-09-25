#include <core/importer_factory.h>

#include "importers_autoreg.h"
#include "ImporterTKA.h"


using namespace DAQuiri;

void importers_autoreg()
{
  DAQUIRI_REGISTER_IMPORTER(ImporterTKA)
}
