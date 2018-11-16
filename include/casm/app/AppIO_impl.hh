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
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"


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

  template<typename OrbitPrinter, typename Element>
  void print_coordinates(OrbitPrinter &printer, const Element &element, Log &out) {
    out << out.indent_str() << "Coordinates:" << std::endl;
    printer.increase_indent(out);
    printer.print(element, out);
    printer.decrease_indent(out);
  }

  template<typename OrbitType>
  void PrinterBase::print_equivalence_map(const OrbitType &orbit, Index equiv_index, Log &out) const {
    out << out.indent_str() << "Equivalence map: " << std::endl;
    this->increase_indent(out);
    int j = 0;
    for(const auto &op : orbit.equivalence_map()[equiv_index]) {
      out << out.indent_str() << j << ": (" << op.index() << ") "
          << brief_description(op, orbit.prototype().prim().lattice(), this->opt.coord_type)
          << std::endl;
      ++j;
    }
    this->decrease_indent(out);
  }

  template<typename OrbitType>
  void PrinterBase::print_equivalence_map(const OrbitType &orbit, Log &out) const {
    out << out.indent_str() << "Orbit equivalence map: " << std::endl;
    this->increase_indent(out);
    for(int i = 0; i < orbit.size(); ++i) {
      out << out.indent_str() << "Element: " << i << std::endl;
      this->increase_indent(out);
      int j = 0;
      for(const auto &op : orbit.equivalence_map()[i]) {
        out << out.indent_str() << j << ": (" << op.index() << ") "
            << brief_description(op, orbit.prototype().prim().lattice(), this->opt.coord_type)
            << std::endl;
        ++j;
      }
      this->decrease_indent(out);
    }
    this->decrease_indent(out);
  }

  template<typename OrbitType, typename Element>
  void PrinterBase::print_invariant_group(const OrbitType &orbit, const Element &element, Log &out) const {
    out << out.indent_str() << "Invariant group:" << std::endl;
    SymGroup invariant_group = make_invariant_subgroup(element, orbit.generating_group(), orbit.sym_compare());
    this->increase_indent(out);
    SymInfoOptions topt = this->opt.sym_info_opt;
    brief_description(out, invariant_group, orbit.prototype().prim().lattice(), topt);
    this->decrease_indent(out);
  }

  // --- OrbitPrinter templates ---


  template<typename _Element>
  template<typename OrbitType>
  void OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO>::operator()(
    const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const {

    out << out.indent_str() << "Prototype" << " of " << orbit.size()
        << " Equivalent " << element_name << " in Orbit " << orbit_index << std::endl;
    this->increase_indent(out);

    if(this->opt.print_coordinates) {
      print_coordinates(*this, orbit.prototype(), out);
    }
    if(this->opt.print_invariant_grp) {
      this->print_invariant_group(orbit, orbit.prototype(), out);
    }
    this->decrease_indent(out);

  }

  /// \brief Print to JSON
  ///
  /// Note: for 'read_clust' to work, "prototype" must be written
  template<typename _Element>
  template<typename OrbitType>
  jsonParser &OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO>::to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) const {
    json.put_obj();
    json["prototype"] = orbit.prototype();
    json["linear_orbit_index"] = orbit_index;
    return json;
  }

  template<typename _Element>
  template<typename OrbitType>
  void OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL>::operator()(const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const {

    for(Index equiv_index = 0; equiv_index != orbit.size(); ++equiv_index) {
      out << out.indent_str() << equiv_index << " of " << orbit.size()
          << " Equivalent " << element_name << " in Orbit " << orbit_index << std::endl;
      this->increase_indent(out);

      if(this->opt.print_coordinates) {
        print_coordinates(*this, orbit[equiv_index], out);
      }
      if(this->opt.print_invariant_grp) {
        this->print_invariant_group(orbit, orbit[equiv_index], out);
      }
      if(this->opt.print_equivalence_map) {
        this->print_equivalence_map(orbit, equiv_index, out);
      }
      this->decrease_indent(out);
    }
  }

  /// \brief Print to JSON
  ///
  /// Note: for 'read_clust' to work, "prototype" must be written
  template<typename _Element>
  template<typename OrbitType>
  jsonParser &OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL>::to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) const {
    json.put_obj();
    json["prototype"] = orbit.prototype();
    json["elements"].put_array(orbit.begin(), orbit.end());
    json["linear_orbit_index"] = orbit_index;
    return json;
  }

  template<typename OrbitType>
  void ProtoFuncsPrinter::operator()(const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const {

    out << out.indent_str() << "Prototype" << " of " << orbit.size()
        << " Equivalent " << element_name << " in Orbit " << orbit_index << std::endl;

    // out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    // out.precision(5);

    COORD_MODE printer_mode(opt.coord_type);

    auto const &clust = orbit.prototype();
    Index np = 0;
    this->increase_indent(out);
    for(const auto &coord : clust) {
      out.indent();
      if(opt.coord_type == INTEGRAL) {
        out << coord;
        out << " ";
        coord.site().site_occupant().print(out);
        out << std::flush;
      }
      else {
        coord.site().print(out);
      }
      out << "  basis_index: " << coord.sublat() << "  clust_index: " << np++ << " ";
      if(opt.delim)
        out << opt.delim;
      out << std::flush;
    }
    this->decrease_indent(out);

    Index func_index = 0;
    for(Index i = 0; i < orbit_index; i++)
      func_index += clex_basis.clust_basis(i, 0).size();

    // From clust:
    out.indent() << "Basis Functions:\n";
    BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
    tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
    this->increase_indent(out);
    for(Index i = 0; i < tbasis.size(); i++) {
      out.indent() << "\\Phi_" << func_index + i << " = " << tbasis[i]->tex_formula() << std::endl;
    }
    this->decrease_indent(out);
    out << "\n\n" << std::flush;

  }


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

  /// \brief Print IntegralCluster orbits
  ///
  /// \param begin,end Range of Orbit<IntegralCluster, SymCompareType>
  /// \param out output stream
  /// \param mode Coordinate output mode
  /// \param printer A functor to control printing for the orbit
  ///
  /// Printer is expected to have:
  /// - \code std::string Printer::indent(); \endcode
  /// - \code void Printer::coord_type(Log& out); \endcode
  /// - \code void Printer::operator()(const Orbit<IntegralCluster, SymCompareType>& orbit, Log& out, Index orbit_index, Index Norbits); \endcode
  ///
  template<typename ClusterOrbitIterator, typename OrbitPrinter>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    Log &out,
    OrbitPrinter printer) {

    printer.coord_type(out);

    out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.ostream().precision(5);

    Index branch = -1;
    Index orbit_index = 0;
    Index Norbits = std::distance(begin, end);

    for(auto it = begin; it != end; ++it) {
      if(it->prototype().size() != branch) {
        branch = it->prototype().size();
        out << out.indent_str() << "** Branch " << branch << " ** " << std::endl;
      }
      printer.increase_indent(out);
      out << out.indent_str() << "** " << orbit_index << " of " << Norbits << " Orbits **"
          << "  Points: " << it->prototype().size()
          << "  Mult: " << it->size()
          << "  MinLength: " << it->prototype().min_length()
          << "  MaxLength: " << it->prototype().max_length() << std::endl;
      printer.increase_indent(out);
      printer(*it, out, orbit_index, Norbits);
      out << std::endl;
      printer.decrease_indent(out);
      printer.decrease_indent(out);
      ++orbit_index;
    }

  }

  /// \brief Print IntegralCluster orbits
  template<typename ClusterOrbitIterator>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    Log &out,
    const OrbitPrinterOptions &opt) {

    //typedef typename ClusterOrbitIterator::container_type container_type;
    //typedef typename container_type::value_type orbit_type;
    typedef typename std::iterator_traits<ClusterOrbitIterator>::value_type orbit_type;
    typedef typename orbit_type::Element Element;

    if(opt.orbit_print_mode == ORBIT_PRINT_MODE::PROTO) {
      OrbitPrinter<Element, ORBIT_PRINT_MODE::PROTO> orbit_printer(opt);
      print_clust(begin, end, out, orbit_printer);
    }
    else if(opt.orbit_print_mode == ORBIT_PRINT_MODE::FULL) {
      OrbitPrinter<Element, ORBIT_PRINT_MODE::FULL> orbit_printer(opt);
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
