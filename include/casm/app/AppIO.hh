#ifndef CASM_AppIO
#define CASM_AppIO

#include <string>
#include <map>
#include <boost/filesystem/path.hpp>
#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/symmetry/SymInfo.hh"

namespace CASM {
  namespace xtal {
    class Structure;
    class AtomPosition;
    class Molecule;
    class Site;
    class UnitCellCoord;
    template<typename CoordType> class BasicStructure;
  }
  using namespace xtal;

  // --- These functions are for casm I/O -----------

  class jsonParser;
  class ClexBasis;
  class HamiltonianModules;
  class ChemicalReference;
  class SymGroup;
  class PrimClex;

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

  /// \brief Print AtomPosition to json after applying affine transformation cart2frac*cart()+trans
  jsonParser &to_json(const AtomPosition &apos, jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &cart2frac);

  /// \brief Read AtomPosition from json and then apply affine transformation cart2frac*cart()
  void from_json(AtomPosition &apos, const jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &frac2cart);

  template<>
  struct jsonConstructor<AtomPosition> {

    /// \brief Read from json [b, i, j, k], using 'unit' for AtomPosition::unit()
    static AtomPosition from_json(const jsonParser &json, Eigen::Matrix3d const &f2c_mat);
  };

  jsonParser &to_json(const Molecule &mol, jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &c2f_mat);

  void from_json(Molecule &mol, const jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat);

  template<>
  struct jsonConstructor<Molecule> {
    static Molecule from_json(const jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat);
  };

  template<>
  struct jsonConstructor<Site> {

    static Site from_json(const jsonParser &json,
                          Lattice const &_home,
                          COORD_TYPE coordtype,
                          std::map<std::string, Molecule> const &mol_map,
                          HamiltonianModules const &_modules);
  };

  jsonParser &to_json(const Site &value,
                      jsonParser &json,
                      COORD_TYPE coordtype);

  void from_json(Site &value,
                 const jsonParser &json,
                 Lattice const &_home,
                 COORD_TYPE coordtype,
                 std::map<std::string, Molecule> const &mol_map,
                 HamiltonianModules const &_modules);



  BasicStructure<Site> read_prim(fs::path filename, double xtal_tol, HamiltonianModules const *_modules = nullptr);

  BasicStructure<Site> read_prim(jsonParser const &json, double xtal_tol, HamiltonianModules const *_modules = nullptr);

  /// \brief Write prim.json to file
  void write_prim(const BasicStructure<Site> &prim, fs::path filename, COORD_TYPE mode);

  /// \brief Write prim.json as JSON
  void write_prim(const BasicStructure<Site> &prim, jsonParser &json, COORD_TYPE mode);



  // --------- SymmetryIO Declarations --------------------------------------------------

  void write_symop(const SymGroup &grp, Index i, jsonParser &j);

  void write_symgroup(const SymGroup &grp, jsonParser &json);


  // --------- ChemicalReference IO Declarations --------------------------------------------------

  ChemicalReference read_chemical_reference(fs::path filename, const BasicStructure<Site> &prim, double tol);

  ChemicalReference read_chemical_reference(const jsonParser &json, const BasicStructure<Site> &prim, double tol);

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

  enum class ORBIT_PRINT_MODE {
    PROTO, FULL
  };

  ENUM_IO_DECL(ORBIT_PRINT_MODE)
  ENUM_JSON_IO_DECL(ORBIT_PRINT_MODE)
  ENUM_TRAITS(ORBIT_PRINT_MODE)

  struct OrbitPrinterOptions {
    int indent_space = 6;
    char delim = '\n';
    int prec = 7;
    COORD_TYPE coord_type = FRAC;
    ORBIT_PRINT_MODE orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
    SymInfoOptions sym_info_opt;
    bool print_coordinates = true;
    bool print_equivalence_map = false;
    bool print_invariant_grp = false;
  };

  jsonParser &to_json(const OrbitPrinterOptions &opt, jsonParser &json);

  /// \brief Read from JSON
  void from_json(OrbitPrinterOptions &opt, const jsonParser &json);

  template<>
  struct jsonConstructor<OrbitPrinterOptions> {
    static OrbitPrinterOptions from_json(const jsonParser &json);
  };


  struct PrinterBase {

    OrbitPrinterOptions opt;


    PrinterBase(const OrbitPrinterOptions &_opt = OrbitPrinterOptions());

    void coord_type(Log &out);

    void increase_indent(Log &out) const {
      out.increase_indent_spaces(opt.indent_space);
    }
    void decrease_indent(Log &out) const {
      out.decrease_indent_spaces(opt.indent_space);
    }

    template<typename OrbitType>
    void print_equivalence_map(const OrbitType &orbit, Index equiv_index, Log &out) const;

    template<typename OrbitType>
    void print_equivalence_map(const OrbitType &orbit, Log &out) const;

    template<typename OrbitType, typename Element>
    void print_invariant_group(const OrbitType &orbit, const Element &element, Log &out) const;

  };

  template<typename _Element>
  struct Printer : public PrinterBase {

    typedef _Element Element;

    Printer(const OrbitPrinterOptions &_opt = OrbitPrinterOptions()) :
      PrinterBase(_opt) {}

    void print(const Element &element, Log &out) const {
      COORD_MODE printer_mode(opt.coord_type);
      out << element;
    }
  };

  template<>
  struct Printer<IntegralCluster> : public PrinterBase {

    typedef IntegralCluster Element;
    static const std::string element_name;

    Printer(const OrbitPrinterOptions &_opt = OrbitPrinterOptions()) :
      PrinterBase(_opt) {}

    void print(const Element &element, Log &out) const;
  };

  typedef Printer<IntegralCluster> SitesPrinter;

  template<typename _Element, ORBIT_PRINT_MODE>
  struct OrbitPrinter {};

  /// \brief Print Orbit<SymCompareType>, including only prototypes
  template<typename _Element>
  struct OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO> : public Printer<_Element> {

    using Printer<_Element>::element_name;
    using Printer<_Element>::print;

    OrbitPrinter(const OrbitPrinterOptions &_opt = OrbitPrinterOptions()) :
      Printer<_Element>(_opt) {}


    template<typename OrbitType>
    void operator()(const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const;

    template<typename OrbitType>
    jsonParser &to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) const;
  };

  template<typename _Element>
  using PrototypePrinter = OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO>;

  typedef PrototypePrinter<IntegralCluster> ProtoSitesPrinter;

  /// \brief Print Orbit<SymCompareType>, including all equivalents
  template<typename _Element>
  struct OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL> : public Printer<_Element> {

    using Printer<_Element>::element_name;
    using Printer<_Element>::print;

    OrbitPrinter(const OrbitPrinterOptions &_opt = OrbitPrinterOptions()) :
      Printer<_Element>(_opt) {}


    template<typename OrbitType>
    void operator()(const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const;

    template<typename OrbitType>
    jsonParser &to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) const;
  };

  template<typename _Element>
  using FullOrbitPrinter = OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL>;

  typedef FullOrbitPrinter<IntegralCluster> FullSitesPrinter;

  /// \brief Print Orbit<SymCompareType> & ClexBasis, including prototypes and prototype basis functions
  struct ProtoFuncsPrinter : public SitesPrinter {

    ClexBasis const &clex_basis;

    std::vector<SubExpressionLabeler> labelers;

    ProtoFuncsPrinter(ClexBasis const &_clex_basis, OrbitPrinterOptions const &_opt = OrbitPrinterOptions());

    /// \brief Print to JSON
    ///
    /// Note: for 'read_clust' to work, "prototype" must be written
    template<typename OrbitType>
    void operator()(const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const;

    template<typename OrbitType>
    jsonParser &to_json(const OrbitType &orbit, jsonParser &json, Index orbit_index, Index Norbits) const;

  };



  /// \brief Print IntegralCluster orbits
  template<typename ClusterOrbitIterator, typename OrbitPrinter>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    Log &out,
    OrbitPrinter printer);

  /// \brief Print IntegralCluster orbits
  template<typename ClusterOrbitIterator>
  void print_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    Log &out,
    const OrbitPrinterOptions &opt = OrbitPrinterOptions());

  // /// \brief Print site basis functions, as for 'casm bset --functions'
  // template<typename ClusterOrbitIterator>
  // void print_site_basis_funcs(
  //   ClusterOrbitIterator begin,
  //   ClusterOrbitIterator end,
  //   const ClexBasis &clex_basis,
  //   Log &out,
  //   COORD_TYPE mode);

  /// \brief Print site basis functions, as for 'casm bset --functions'
  void print_site_basis_funcs(
    Structure const &prim,
    ClexBasis const &clex_basis,
    Log &out,
    Index indent_space = 6,
    COORD_TYPE mode = FRAC);

  void write_site_basis_funcs(
    Structure const &prim,
    ClexBasis const &clex_basis,
    jsonParser &json);


  // ---------- clust.json IO ------------------------------------------------------------------

  /// \brief Read JSON containing Orbit<SymCompareType> prototypes
  template<typename ClusterOutputIterator, typename SymCompareType>
  ClusterOutputIterator read_clust(
    ClusterOutputIterator result,
    const jsonParser &json,
    const Structure &prim,
    const SymGroup &generating_grp,
    const SymCompareType &sym_compare,
    double xtal_tol);

  /// \brief Write Orbit<SymCompareType> to JSON, including 'bspecs'
  template<typename ClusterOrbitIterator, typename Printer>
  jsonParser &write_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    jsonParser &json,
    Printer printer);

  /// \brief Write Orbit<SymCompareType> to JSON, including 'bspecs'
  template<typename ClusterOrbitIterator, typename Printer>
  jsonParser &write_clust(
    ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    jsonParser &json,
    Printer printer,
    const jsonParser &bspecs);


  // --- Read Selections & Object names ----------------------------------------

  namespace DB {
    template<typename T> class Selection;
  }

  /// \brief Make a DB::Selection from JSON input
  template<typename DataObject>
  DB::Selection<DataObject> make_selection(
    const PrimClex &primclex,
    const jsonParser &kwargs,
    std::string name_key,
    std::string sel_key,
    std::string method_name,
    OnError on_error);
}

#endif
