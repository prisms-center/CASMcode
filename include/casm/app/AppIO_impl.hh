#ifndef CASM_AppIO_impl
#define CASM_AppIO_impl

#include "casm/app/AppIO.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/PrimClex.hh"

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/clex.hh"
#include "casm/casm_io/SafeOfstream.hh"

#include "casm/database/Selection_impl.hh"


namespace CASM {

  // --------- CompositionAxes Definitions --------------------------------------------------

  /// \brief Read standard axes from JSON, and output to std::map<std::string, CompositionConverter>
  ///
  /// - json Format: {"standard_axes": {"key0" : {CompositionConverter}, ... more axes ...},
  ///                 "custom_axes":   {"key0" : {CompositionConverter}, ... more axes ...},
  ///                 "current_axes": "key"}
  ///   - "axes_type" is either "standard_axes" or "custom_axes"
  ///   - "key" is a string indicating which composition axes in the "standard_axes" or "custom_axes" JSON object
  ///
  /// - This read json["standard_axes"] or json["custom_axes"]
  ///
  template<typename OutputIterator>
  OutputIterator read_composition_axes(OutputIterator result, const jsonParser &json) {

    CompositionConverter conv;
    for(auto it = json.cbegin(); it != json.cend(); ++it) {
      from_json(conv, *it);
      *result++ = std::make_pair(it.name(), conv);
    }
    return result;
  }


  // ---------- casm bspecs --orbits -----------------------------------------------------------



  /// \brief Print IntegralCluster orbits
  ///
  /// \param begin,end Range of Orbit<IntegralCluster, SymCompareType>
  /// \param out output stream
  /// \param mode Coordinate output mode
  /// \param printer A functor to control printing for the orbit
  ///
  /// Printer is expected to have:
  /// - \code std::string Printer::indent(); \endcode
  /// - \code void Printer::coord_mode(std::ostream& out); \endcode
  /// - \code void Printer::operator()(const Orbit<IntegralCluster, SymCompareType>& orbit, std::ostream& out, Index orbit_index, Index Norbits); \endcode
  ///
  template<typename ClusterOrbitIterator, typename OrbitPrinter>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    std::ostream &out,
    OrbitPrinter printer) {

    printer.coord_mode(out);

    out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.precision(5);

    Index branch = -1;
    Index orbit_index = 0;
    Index Norbits = std::distance(begin, end);
    std::string indent = printer.indent();

    for(auto it = begin; it != end; ++it) {
      if(it->prototype().size() != branch) {
        branch = it->prototype().size();
        out << "** Branch " << branch << " ** " << std::endl;
      }
      out << printer.indent() << "** " << orbit_index << " of " << Norbits << " Orbits **"
          << "  Points: " << it->prototype().size()
          << "  Mult: " << it->size()
          << "  MinLength: " << it->prototype().min_length()
          << "  MaxLength: " << it->prototype().max_length() << std::endl;
      printer(*it, out, orbit_index, Norbits);
      out << std::endl;
      ++orbit_index;
    }

  }

  // ---------- clust.json IO ------------------------------------------------------------------

  /// \brief Read JSON containing Orbit<IntegralCluster, SymCompareType> prototypes
  ///
  /// - Uses 'prim', 'generating_grp', and 'sym_compare' to generate orbits from
  ///   prototypes read from the JSON
  /// - Ignores "prim" and "bspecs" info in the JSON
  ///
  template<typename ClusterOutputIterator, typename SymCompareType>
  ClusterOutputIterator read_clust(
    ClusterOutputIterator result,
    const jsonParser &json,
    const Structure &prim,
    const SymGroup &generating_grp,
    const SymCompareType &sym_compare,
    double xtal_tol) {

    typedef Orbit<IntegralCluster, SymCompareType> orbit_type;

    for(const auto &j : json["orbits"]) {
      *result++ = orbit_type(j["prototype"].get<IntegralCluster>(prim, xtal_tol), generating_grp, sym_compare);
    }
    return result;
  }

  /// \brief Write Orbit<IntegralCluster, SymCompareType> to JSON, including 'bspecs'
  ///
  /// Format:
  /// \code
  /// {
  ///   "orbits": [
  ///     {
  ///       "prototype" : {
  ///         "min_length" : number,
  ///         "max_length" : number,
  ///         "sites" : (JSON array of UnitCellCoord)
  ///       }
  ///     },
  ///     ... for each orbit ...
  ///   ],
  ///   "prim" : (JSON object, the contents of prim.json)
  /// }
  /// \endcode
  ///
  template<typename ClusterOrbitIterator, typename Printer>
  jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, jsonParser &json, Printer printer) {
    json = jsonParser::object();
    Index Norbits = std::distance(begin, end);
    json["orbits"] = jsonParser::array(Norbits, jsonParser::object());
    Index orbit_index = 0;
    for(auto it = begin; it != end; ++it, ++orbit_index) {
      printer.to_json(*it, json["orbits"][orbit_index], orbit_index, Norbits);
    }
    write_prim(begin->prototype().prim(), json["prim"], FRAC);
    return json;
  }

  /// \brief Write Orbit<IntegralCluster, SymCompareType> to JSON, including 'bspecs'
  ///
  /// Format:
  /// \code
  /// {
  ///   "orbits": [
  ///     {
  ///       "prototype" : {
  ///         "min_length" : number,
  ///         "max_length" : number,
  ///         "sites" : (JSON array of UnitCellCoord)
  ///       }
  ///     },
  ///     ... for each orbit ...
  ///   ],
  ///   "bspecs" : (JSON object, the contents of bspecs.json)
  ///   "prim" : (JSON object, the contents of prim.json)
  /// }
  /// \endcode
  ///
  template<typename ClusterOrbitIterator, typename Printer>
  jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, jsonParser &json, Printer printer, const jsonParser &bspecs) {
    write_clust(begin, end, json, printer);
    json["bspecs"] = bspecs;
    return json;
  }


  // --- Read Selections & Object names ----------------------------------------

  /// \brief Make a DB::Selection from JSON input
  ///
  /// \param primclex PrimClex
  /// \param kwargs jsonParser with JSON input
  /// \param name_key Read object names from kwargs[name_key]. Expects array of string.
  /// \param sel_key Read selection name for kwargs[sel_key]. Expect string.
  /// \param method_name Method name for error messages. Ex. "DiffTransConfigEnumOccPerturbations"
  /// \param on_error Indicates how to handle names that do not exist in database.
  ///
  /// Notes:
  /// - Reads selection from 'sel_key' first, using "NONE" if 'sel_key' does not
  ///   exist. Then reads names from 'name_key' array and adds them to the selection.
  /// - If 'name_key' and 'sel_key' do not exist, throw exception.
  ///
  template<typename DataObject>
  DB::Selection<DataObject> make_selection(
    const PrimClex &primclex,
    const jsonParser &kwargs,
    std::string name_key,
    std::string sel_key,
    std::string method_name,
    OnError on_error) {

    if(!kwargs.contains(name_key) && !kwargs.contains(sel_key)) {
      std::string msg = "Error in " + method_name + ": One of " + name_key +
                        " or " + sel_key + " must be given.";
      throw std::runtime_error(msg);
    }

    std::vector<std::string> obj_names;
    std::string sel_name;
    kwargs.get_else(obj_names, name_key, std::vector<std::string>());
    kwargs.get_else(sel_name, sel_key, std::string("NONE"));

    DB::Selection<DataObject> sel(primclex, sel_name);
    for(const auto &name : obj_names) {

      // validate obj_names & handle errors
      if(!sel.db().count(name)) {
        if(on_error == OnError::THROW) {
          std::string msg = "Error in " + method_name + ": \"" + name + "\" is "
                            "not the project database.";
          throw std::runtime_error(msg);
        }
        else if(on_error == OnError::WARN) {
          std::string msg = "Warning in " + method_name + ": \"" + name + "\" is "
                            "not the project database. Skipping...\n";
          primclex.err_log() << msg;
          continue;
        }
        else if(on_error == OnError::CONTINUE) {
          continue;
        }

      }

      sel.data()[name] = true;
    }
    return sel;
  }

}

#endif