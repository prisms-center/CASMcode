#ifndef CASM_AppIO
#define CASM_AppIO

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ChemicalReference.hh"

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/clex.hh"
#include "casm/casm_io/SafeOfstream.hh"

namespace CASM {

  // --- These functions are for casm I/O -----------

  class ClexBasis;

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

  /// \brief Print Orbit<IntegralCluster>, including only prototypes
  struct ProtoSitesPrinter : public SitesPrinter {

    ProtoSitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC);

    void operator()(const Orbit<IntegralCluster> &orbit, std::ostream &out, Index orbit_index, Index Norbits);
  };

  /// \brief Print Orbit<IntegralCluster>, including all equivalents
  struct FullSitesPrinter : public SitesPrinter {

    FullSitesPrinter(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC);

    void operator()(const Orbit<IntegralCluster> &orbit, std::ostream &out, Index orbit_index, Index Norbits);
  };

  /// \brief Print Orbit<IntegralCluster> & ClexBasis, including prototypes and prototype basis functions
  struct ProtoFuncsPrinter : public SitesPrinter {

    const ClexBasis &clex_basis;

    ProtoFuncsPrinter(const ClexBasis &_clex_basis, int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = FRAC);

    void operator()(const Orbit<IntegralCluster> &orbit, std::ostream &out, Index orbit_index, Index Norbits);
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

  /// \brief Print clust.json file
  template<typename ClusterOrbitIterator>
  jsonParser &write_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    jsonParser &bspecs,
    jsonParser &json);


  // ---------- basis.json IO ------------------------------------------------------------------


  /// \brief Write summary of cluster expansion basis
  template<typename ClusterOrbitIterator>
  void write_basis(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    const ClexBasis &clex_basis,
    jsonParser &json,
    double tol);
}

#endif
