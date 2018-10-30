#ifndef CASM_AppIO_impl
#define CASM_AppIO_impl

#include "casm/app/AppIO.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/basis_set/FunctionVisitor.hh"

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

  template<typename OrbitType>
  jsonParser &ProtoFuncsPrinter::to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) const {
    json.put_obj();
    json["prototype"] = orbit.prototype();
    json["linear_orbit_index"] = orbit_index;
    json["mult"] = orbit.size();

    jsonParser &orbitf = json["cluster_functions"];
    orbitf = jsonParser::array();



    // basis function info
    Index func_index = 0;
    for(Index i = 0; i < orbit_index; i++)
      func_index += clex_basis.clust_basis(i, 0).size();

    BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
    tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
    for(auto const &labeler : labelers) {
      tbasis.accept(labeler);
    }

    for(Index nf = 0; nf < tbasis.size(); ++nf) {
      orbitf.push_back(json_pair("\\Phi_" + std::to_string(func_index + nf), tbasis[nf]->tex_formula()));
    }

    return json;
  }

  template<typename OrbitType>
  void ProtoFuncsPrinter::operator()(const OrbitType &orbit, Log &_out, Index orbit_index, Index Norbits) const {
    std::ostream &out = _out.ostream();
    out << indent() << indent() << "Prototype" << " of " << orbit.size()
        << " Equivalent " << element_name << " in Orbit " << orbit_index << std::endl;

    out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.precision(5);

    COORD_TYPE _mode = mode;
    if(_mode == COORD_DEFAULT) {
      _mode = COORD_MODE::CHECK();
    }
    COORD_MODE printer_mode(_mode);

    auto const &clust = orbit.prototype();
    Index np = 0;
    for(const auto &coord : clust) {
      out << indent() << indent() << indent();
      if(_mode == INTEGRAL) {
        out << coord;
        out << " ";
        coord.site().site_occupant().print(out);
        out << std::flush;
      }
      else {
        out.setf(std::ios::showpoint, std::ios_base::fixed);
        out.precision(5);
        out.width(9);
        coord.site().print(out);
      }
      out << "  basis_index: " << coord.sublat() << "  clust_index: " << np++ << " ";
      if(delim)
        out << delim;
      out << std::flush;
    }

    Index func_index = 0;
    for(Index i = 0; i < orbit_index; i++)
      func_index += clex_basis.clust_basis(i, 0).size();

    // From clust:
    out << "            Basis Functions:\n";
    BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
    tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
    for(Index i = 0; i < tbasis.size(); i++) {
      out << "              \\Phi_" << func_index + i << " = " << tbasis[i]->tex_formula() << std::endl;
    }
    out << "\n\n" << std::flush;

  }

  /// \brief Print IntegralCluster orbits
  ///
  /// \param begin,end Range of Orbit<IntegralCluster, SymCompareType>
  /// \param out output stream
  /// \param mode Coordinate output mode
  /// \param printer A functor to control printing for the orbit
  ///
  /// Printer is expected to have:
  /// - \code std::string Printer::indent(); \endcode
  /// - \code void Printer::coord_mode(Log& out); \endcode
  /// - \code void Printer::operator()(const Orbit<IntegralCluster, SymCompareType>& orbit, Log& out, Index orbit_index, Index Norbits); \endcode
  ///
  template<typename ClusterOrbitIterator, typename OrbitPrinter>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    Log &out,
    OrbitPrinter printer) {

    printer.coord_mode(out);

    out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.ostream().precision(5);

    Index branch = -1;
    Index orbit_index = 0;
    Index Norbits = std::distance(begin, end);
    std::string indent = printer.indent();

    for(auto it = begin; it != end; ++it) {
      if(it->prototype().size() != branch) {
        branch = it->prototype().size();
        out << out.indent_str() << "** Branch " << branch << " ** " << std::endl;
      }
      printer.indent_level++;
      out << out.indent_str() << printer.indent() << "** " << orbit_index << " of " << Norbits << " Orbits **"
          << "  Points: " << it->prototype().size()
          << "  Mult: " << it->size()
          << "  MinLength: " << it->prototype().min_length()
          << "  MaxLength: " << it->prototype().max_length() << std::endl;
      printer.indent_level++;
      printer(*it, out, orbit_index, Norbits);
      out << std::endl;
      printer.indent_level--;
      printer.indent_level--;
      ++orbit_index;
    }

  }

  /// \brief Print IntegralCluster orbits
  template<typename ClusterOrbitIterator>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    Log &out,
    ORBIT_PRINT_MODE _orbit_print_mode,
    COORD_TYPE _coord_mode,
    int _indent_space,
    char _delim) {

    //typedef typename ClusterOrbitIterator::container_type container_type;
    //typedef typename container_type::value_type orbit_type;
    typedef typename std::iterator_traits<ClusterOrbitIterator>::value_type orbit_type;
    typedef typename orbit_type::Element Element;

    if(_orbit_print_mode == ORBIT_PRINT_MODE::PROTO) {
      OrbitPrinter<Element, ORBIT_PRINT_MODE::PROTO> orbit_printer(_indent_space, _delim, _coord_mode);
      print_clust(begin, end, out, orbit_printer);
    }
    else if(_orbit_print_mode == ORBIT_PRINT_MODE::FULL) {
      OrbitPrinter<Element, ORBIT_PRINT_MODE::FULL> orbit_printer(_indent_space, _delim, _coord_mode);
      print_clust(begin, end, out, orbit_printer);
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
