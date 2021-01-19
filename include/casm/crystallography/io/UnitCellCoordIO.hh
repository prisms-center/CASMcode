#ifndef UNITCELLCOORDIO_HH
#define UNITCELLCOORDIO_HH

#include "casm/casm_io/json/jsonParser.hh"
namespace CASM {
namespace xtal {
class UnitCellCoord;
}

/// \brief Print to json as [b, i, j, k]
jsonParser &to_json(const xtal::UnitCellCoord &ucc_val, jsonParser &fill_json);

/// \brief Read from json [b, i, j, k]
void from_json(xtal::UnitCellCoord &fill_value, const jsonParser &read_json);
}  // namespace CASM

#endif
