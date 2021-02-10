#include "casm/symmetry/io/json/SymInfo_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"

namespace CASM {

ENUM_JSON_IO_DEF(symmetry_type)

jsonParser &to_json(const SymInfoOptions &opt, jsonParser &json) {
  json.put_obj();
  json[traits<COORD_TYPE>::name] = opt.coord_type;
  json["tol"] = opt.tol;
  json["prec"] = opt.prec;
  json["print_matrix_tau"] = opt.print_matrix_tau;
  return json;
}

/// \brief Read from JSON
void from_json(SymInfoOptions &opt, const jsonParser &json) {
  json.get_if(opt.coord_type, traits<COORD_TYPE>::name);
  json.get_if(opt.tol, "tol");
  json.get_if(opt.prec, "prec");
  json.get_if(opt.print_matrix_tau, "print_matrix_tau");
}

SymInfoOptions jsonConstructor<SymInfoOptions>::from_json(
    const jsonParser &json) {
  SymInfoOptions res;
  CASM::from_json(res, json);
  return res;
}

/// \brief Adds to existing JSON object
void to_json(const SymInfo &info, jsonParser &j) {
  /// type of symmetry, given by one of the allowed values of symmetry_type
  j["type"] = info.op_type;

  if (info.op_type == symmetry_type::rotation_op ||
      info.op_type == symmetry_type::screw_op ||
      info.op_type == symmetry_type::rotoinversion_op) {
    to_json_array(info.axis.const_cart().normalized(),
                  j["rotation_axis"]["CART"]);
    to_json_array(info.axis.const_frac().normalized(),
                  j["rotation_axis"]["FRAC"]);
    j["rotation_angle"] = info.angle;
  } else if (info.op_type == symmetry_type::mirror_op ||
             info.op_type == symmetry_type::glide_op) {
    to_json_array(info.axis.const_cart().normalized(),
                  j["mirror_normal"]["CART"]);
    to_json_array(info.axis.const_frac().normalized(),
                  j["mirror_normal"]["FRAC"]);
  }

  if (info.op_type == symmetry_type::screw_op ||
      info.op_type == symmetry_type::glide_op) {
    to_json_array(info.screw_glide_shift.const_cart(), j["shift"]["CART"]);
    to_json_array(info.screw_glide_shift.const_frac(), j["shift"]["FRAC"]);
  }

  if (info.op_type != symmetry_type::identity_op) {
    to_json_array(info.location.const_cart(), j["invariant_point"]["CART"]);
    to_json_array(info.location.const_frac(), j["invariant_point"]["FRAC"]);
  }
}

}  // namespace CASM
