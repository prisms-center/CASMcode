#ifndef CASM_ConfigImport
#define CASM_ConfigImport

namespace CASM {

  struct ConfigImportData {

    fs::path structure_path;
    fs::path properties_calc_path;

    ConfigMappingResults map_result;

  };

}

#endif
