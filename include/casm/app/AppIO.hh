#ifndef CASM_AppIO
#define CASM_AppIO

#include <string>
#include <map>
#include <boost/filesystem/path.hpp>
#include "casm/CASM_global_definitions.hh"
#include "casm/CASM_global_enum.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/clusterography/ClusterDecl.hh"

namespace CASM {

  // --- These functions are for casm I/O -----------

  template<typename CoordType> class BasicStructure;
  class Site;
  class jsonParser;
  class ClexBasis;
  class ChemicalReference;
  class UnitCellCoord;
  class SymGroup;

  /** \defgroup ProjectIO
   *
   *  \ingroup Project
   *  \ingroup casmIO
   *
   *  \brief Relates to CASM project input/output
   *
   *  @{
  */


  // --------- PrimIO Declarations --------------------------------------------------

  BasicStructure<Site> read_prim(fs::path filename);

  BasicStructure<Site> read_prim(const jsonParser &json);

  /// \brief Write prim.json to file
  void write_prim(const BasicStructure<Site> &prim, fs::path filename, COORD_TYPE mode);

  /// \brief Write prim.json as JSON
  void write_prim(const BasicStructure<Site> &prim, jsonParser &json, COORD_TYPE mode);



  // --------- SymmetryIO Declarations --------------------------------------------------

  void write_symop(const SymGroup &grp, Index i, jsonParser &j);

  void write_symgroup(const SymGroup &grp, jsonParser &json);


  // --------- ChemicalReference IO Declarations --------------------------------------------------

  ChemicalReference read_chemical_reference(fs::path filename, const Structure &prim, double tol);

  ChemicalReference read_chemical_reference(const jsonParser &json, const Structure &prim, double tol);

  void write_chemical_reference(const ChemicalReference &chem_ref, fs::path filename);

  void write_chemical_reference(const ChemicalReference &chem_ref, jsonParser &json);

  void write_chemical_reference(const ChemicalReference &chem_ref, jsonParser &json);


  // --------- CompositionAxes Declarations --------------------------------------------------

  struct CompositionAxes {

    CompositionAxes() {}

    /// \brief Read CompositionAxes from file
    CompositionAxes(fs::path _filename);

    /// \brief Read CompositionAxes from JSON
    CompositionAxes(const jsonParser &json);

    /// \brief Read CompositionAxes from file
    void read(fs::path _filename);

    /// \brief Read CompositionAxes from JSON
    void read(const jsonParser &json, fs::path filename = fs::path());

    /// \brief Write CompositionAxes to file
    void write(fs::path _filename);

    /// \brief Write CompositionAxes to JSON
    void write(jsonParser &json) const;

    /// \brief Set this->curr using key
    void select(std::string key);


    std::map<std::string, CompositionConverter> standard;
    std::map<std::string, CompositionConverter> custom;
    bool has_current_axes = false;
    std::string curr_key;
    CompositionConverter curr;
    fs::path filename;

    int err_code = 0;
    std::string err_message;
  };


  /// \brief Read standard axes from JSON, and output to std::map<std::string, CompositionConverter>
  template<typename OutputIterator>
  OutputIterator read_composition_axes(OutputIterator result, const jsonParser &json);


  // ---------- casm bset --orbits, --clusters, --functions --------------------------------------

  struct PrinterBase {

    int indent_space;
    char delim;
    COORD_TYPE mode;


    PrinterBase(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC);

    std::string indent() const;

    void coord_mode(std::ostream &out);

  };

  template<typename _Element>
  struct Printer : public PrinterBase {

    typedef _Element Element;

    Printer(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      PrinterBase(_indent_space, _delim, _mode) {}

    void print(const Element &element, std::ostream &out) {
      COORD_MODE printer_mode(mode);
      out << element;
    }
  };

  template<>
  struct Printer<IntegralCluster> : public PrinterBase {

    typedef IntegralCluster Element;
    static const std::string element_name;

    Printer(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      PrinterBase(_indent_space, _delim, _mode) {}

    void print(const Element &element, std::ostream &out);
  };

  typedef Printer<IntegralCluster> SitesPrinter;

  /// \brief Print Orbit<IntegralCluster, SymCompareType>, including only prototypes
  template<typename _Element>
  struct PrototypePrinter : public Printer<_Element> {

    using Printer<_Element>::indent;
    using Printer<_Element>::element_name;
    using Printer<_Element>::print;

    PrototypePrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      Printer<_Element>(_indent_space, _delim, _mode) {}


    template<typename OrbitType>
    void operator()(const OrbitType &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      out << indent() << indent() << "Prototype" << " of " << orbit.size()
          << " Equivalent " << element_name << " in Orbit " << orbit_index << std::endl;
      print(orbit.prototype(), out);
    }

    /// \brief Print to JSON
    ///
    /// Note: for 'read_clust' to work, "prototype" must be written
    template<typename OrbitType>
    jsonParser &to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) {
      json.put_obj();
      json["prototype"] = orbit.prototype();
      json["linear_orbit_index"] = orbit_index;
      return json;
    }
  };

  typedef PrototypePrinter<IntegralCluster> ProtoSitesPrinter;

  /// \brief Print Orbit<IntegralCluster, SymCompareType>, including all equivalents
  template<typename _Element>
  struct FullOrbitPrinter : public Printer<_Element> {

    using Printer<_Element>::indent;
    using Printer<_Element>::element_name;
    using Printer<_Element>::print;

    FullOrbitPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      Printer<_Element>(_indent_space, _delim, _mode) {}


    template<typename OrbitType>
    void operator()(const OrbitType &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      for(Index equiv_index = 0; equiv_index != orbit.size(); ++equiv_index) {
        out << indent() << indent() << equiv_index << " of " << orbit.size()
            << " Equivalent " << element_name << " in Orbit " << orbit_index << std::endl;
        print(orbit.prototype(), out);
      }
    }

    /// \brief Print to JSON
    ///
    /// Note: for 'read_clust' to work, "prototype" must be written
    template<typename OrbitType>
    jsonParser &to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) {
      json.put_obj();
      json["prototype"] = orbit.prototype();
      json["elements"].put_array(orbit.begin(), orbit.end());
      json["linear_orbit_index"] = orbit_index;
      return json;
    }
  };

  typedef FullOrbitPrinter<IntegralCluster> FullSitesPrinter;

  /// \brief Print Orbit<IntegralCluster, SymCompareType> & ClexBasis, including prototypes and prototype basis functions
  struct ProtoFuncsPrinter : public SitesPrinter {

    const ClexBasis &clex_basis;

    ProtoFuncsPrinter(const ClexBasis &_clex_basis, int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      SitesPrinter(_indent_space, _delim, _mode),
      clex_basis(_clex_basis) {}

    /// \brief Print to JSON
    ///
    /// Note: for 'read_clust' to work, "prototype" must be written
    template<typename OrbitType>
    void operator()(const OrbitType &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      out << indent() << indent() << "Prototype" << " of " << orbit.size()
          << " Equivalent Clusters in Orbit " << orbit_index << std::endl;
      print(orbit.prototype(), out);

      throw std::runtime_error("Error printing basis functions: ProtoFuncsPrinter not implemented");
      //print_clust_basis(out, nf, 8, '\n');
      //nf += prototype(i, j).clust_basis.size();
      //out << "\n\n" << std::flush;
    }

    template<typename OrbitType>
    jsonParser &to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) {
      json.put_obj();
      json["prototype"] = orbit.prototype();
      json["linear_orbit_index"] = orbit_index;

      throw std::runtime_error("Error printing basis functions: ProtoFuncsPrinter not implemented");

      return json;
    }

  };



  /// \brief Print IntegralCluster orbits
  template<typename ClusterOrbitIterator, typename OrbitPrinter>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    std::ostream &out,
    OrbitPrinter printer);

  /// \brief Print site basis functions, as for 'casm bset --functions'
  template<typename ClusterOrbitIterator>
  void print_site_basis_funcs(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    const ClexBasis &clex_basis,
    std::ostream &out,
    COORD_TYPE mode);


  // ---------- clust.json IO ------------------------------------------------------------------

  /// \brief Read JSON containing Orbit<IntegralCluster, SymCompareType> prototypes
  template<typename ClusterOutputIterator, typename SymCompareType>
  ClusterOutputIterator read_clust(
    ClusterOutputIterator result,
    const jsonParser &json,
    const Structure &prim,
    const SymGroup &generating_grp,
    const SymCompareType &sym_compare,
    double xtal_tol);

  /// \brief Write Orbit<IntegralCluster, SymCompareType> to JSON, including 'bspecs'
  template<typename ClusterOrbitIterator, typename Printer>
  jsonParser &write_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    jsonParser &json,
    Printer printer);

  /// \brief Write Orbit<IntegralCluster, SymCompareType> to JSON, including 'bspecs'
  template<typename ClusterOrbitIterator, typename Printer>
  jsonParser &write_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    jsonParser &json,
    Printer printer,
    const jsonParser &bspecs);

}

#endif
