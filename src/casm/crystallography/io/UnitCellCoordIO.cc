#include "casm/crystallography/io/UnitCellCoordIO.hh"

#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
/// \brief Print to json as [b, i, j, k]
jsonParser &to_json(const xtal::UnitCellCoord &ucc_val, jsonParser &fill_json) {
  fill_json.put_array();
  fill_json.push_back(ucc_val.sublattice());
  fill_json.push_back(ucc_val.unitcell()(0));
  fill_json.push_back(ucc_val.unitcell()(1));
  fill_json.push_back(ucc_val.unitcell()(2));

  return fill_json;
}

/// \brief Read from json [b, i, j, k], assuming fill_value.unit() is already
/// set
void from_json(xtal::UnitCellCoord &fill_value, const jsonParser &read_json) {
  auto b = read_json[0].get<Index>();
  auto i = read_json[1].get<Index>();
  auto j = read_json[2].get<Index>();
  auto k = read_json[3].get<Index>();

  fill_value = xtal::UnitCellCoord(b, i, j, k);

  return;
}
}  // namespace CASM
