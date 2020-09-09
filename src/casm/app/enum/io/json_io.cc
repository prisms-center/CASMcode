#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/enum/io/json_io.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/enumerator/Enumerator.hh"
#include "casm/global/enum/json_io.hh"

namespace CASM {

  /// Parse `casm enum` options
  ///
  /// \code
  /// {
  ///    "filter": <string, optional>,
  ///    "dry_run": <bool, defaul=false>,
  ///    "coordinate_mode": <string, one of: "FRAC", "Direct", "direct", "Fractional",
  ///                                        "fractional", "CART", "Cartesian", "cartesian",
  ///                                        "INTEGRAL", "Integral", "integral">,
  ///    "primitive_only": <bool, default=true>,
  ///    "verbosity": <int or string, default=10, range=[0,100], "none" = 0, "quiet"=5,
  ///                  "standard"=10, "verbose"=20, "debug"=100>
  /// }
  /// \endcode
  ///
  void parse(InputParser<InsertAlgorithmsOptions> &parser) {
    parser.value = notstd::make_unique<InsertAlgorithmsOptions>();
    auto &insert_options = *parser.value;

    if(parser.self.contains("filter")) {
      std::string filter_expr;
      parser.require(filter_expr, "filter");
      insert_options.filter_expr.push_back(filter_expr);
    }
    parser.optional_else(insert_options.dry_run, "dry_run", false);
    parser.optional_else(insert_options.coord_type, traits<COORD_TYPE>::name, COORD_TYPE::FRAC);
    parser.optional_else(insert_options.primitive_only, "primitive_only", true);
    parse_verbosity(parser);
    parser.optional_else(insert_options.verbosity, "verbosity", Log::standard);
  }

}
