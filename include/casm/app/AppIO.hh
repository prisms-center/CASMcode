#ifndef CASM_AppIO
#define CASM_AppIO

#include <string>
#include <map>

#include "casm/CASM_global_definitions.hh"
#include "casm/clex/CompositionConverter.hh"

namespace CASM {

  // --- These functions are for casm I/O -----------

  template<typename CoordType> class BasicStructure;
  class Site;
  class jsonParser;
  class ClexBasis;
  class ChemicalReference;
  template<typename CoordType> class CoordCluster;
  class UnitCellCoord;
  typedef CoordCluster<UnitCellCoord> IntegralCluster;
  template<typename Element, typename SymCompareType> class Orbit;


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

  struct SitesPrinter {

    int indent_space;
    char delim;
    COORD_TYPE mode;


    SitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC);

    std::string indent() const;

    void coord_mode(std::ostream &out);

    void print_sites(const IntegralCluster &clust, std::ostream &out);
  };

  /// \brief Print Orbit<IntegralCluster, SymCompareType>, including only prototypes
  struct ProtoSitesPrinter : public SitesPrinter {

    ProtoSitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      SitesPrinter(_indent_space, _delim, _mode) {}


    template<typename OrbitType>
    void operator()(const OrbitType &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      out << indent() << indent() << "Prototype" << " of " << orbit.size()
          << " Equivalent Clusters in Orbit " << orbit_index << std::endl;
      print_sites(orbit.prototype(), out);
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

  /// \brief Print Orbit<IntegralCluster, SymCompareType>, including all equivalents
  struct FullSitesPrinter : public SitesPrinter {

    FullSitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC) :
      SitesPrinter(_indent_space, _delim, _mode) {}


    template<typename OrbitType>
    void operator()(const OrbitType &orbit, std::ostream &out, Index orbit_index, Index Norbits) {
      for(Index equiv_index = 0; equiv_index != orbit.size(); ++equiv_index) {
        out << indent() << indent() << equiv_index << " of " << orbit.size()
            << " Equivalent Clusters in Orbit " << orbit_index << std::endl;
        print_sites(orbit.prototype(), out);
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
      print_sites(orbit.prototype(), out);

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
    jsonParser &json,
    const Structure &prim,
    const SymGroup &generating_grp,
    const SymCompareType &sym_compare);

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
