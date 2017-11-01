#ifndef CASM_DiffTransConfigMapping
#define CASM_DiffTransConfigMapping

#include <vector>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {
  class Supercell;
  class Lattice;
  class SymGroup;
  class PrimClex;
  class Configuration;
  class ConfigDoF;

  namespace Kinetics {
    class DiffTransConfiguration;
  }

  namespace Completer {
    class ImportOption;
  }

  /// Data structure holding results of ConfigMapper algorithm
  struct DiffTransConfigMapperResult {

    DiffTransConfigMapperResult() :
      success(false) {}

    /// Output structures, after applying lattice similarity and/or rotation to
    /// input structures.
    std::vector<BasicStructure<Site>> structures;

    /// The configuration the input structure was mapped onto
    std::unique_ptr<Kinetics::DiffTransConfiguration> config;

    /// relaxation_properties is populated by relaxation properties:
    ///
    /// - 'lattice_deformation': lattice mapping score
    /// - 'basis_deformation': atomic mapping score
    /// - 'volume_relaxation': V/V_ideal
    /// - 'relaxation_deformation': 3x3 tensor describing cell relaxation
    /// - 'relaxation_displacement': Nx3 matrix describing basis displacements
    /// - 'relaxed_energy': the energy of the relaxed configuration
    jsonParser relaxation_properties;

    /// best_assignment is populated by the permutation of sites in the imported
    /// structure that maps them onto sites of the ideal crystal (excluding vacancies)
    std::vector<Index> best_assignment;

    /// cart_op is populated by the cartesian isometry that rotates the imported
    /// structure onto the coordinate system of the ideal crystal
    Eigen::Matrix3d cart_op;

    /// True if could map to prim, false if not
    bool success;

    /// Failure message if could not map to prim
    std::string fail_msg;

  };


  /// A class for mapping an arbitrary list of crystal structures as a difftransconfiguration of a crystal template
  /// as described by a PrimClex.  DiffTransConfigMapper manages options for the mapping algorithm and mapping cost function
  /// It also caches some information about supercell lattices so that batch imports are more efficient
  ///
  /// \ingroup DiffTransConfiguration
  class DiffTransConfigMapper {
  public:
    enum NullInitializer {null_initializer};
    enum Options {none = 0,
                  rotate = (1u << 0),
                  strict = (1u << 1),
                  robust = (1u << 2)
                 };

    ///\brief Default construction not allowed -- this constructor provides an override
    DiffTransConfigMapper(NullInitializer) :
      m_pclex(nullptr),
      m_lattice_weight(0.5),
      m_max_volume_change(0.5),
      m_min_va_frac(0.),
      m_max_va_frac(1.) {
    }

    ///\brief Construct and initialize a DiffTransConfigMapper
    ///\param _pclex the PrimClex that describes the crystal template
    ///
    ///\param _lattice_weight
    ///\parblock
    ///          free parameter 'w' in the cost function: total_cost = w*lattice_deformation+(1-w)*basis_deformation
    ///          can vary between 0 (completely basis-focused) and 1 (completely lattice-focused)
    ///\endparblock
    ///
    ///\param _max_volume_change
    ///\parblock
    ///          constrains the search space by assuming a limit on allowed volume change
    ///          only taken into account when non-interstitial vacancies are allowed
    ///\endparblock
    ///
    ///\param _options
    ///\parblock
    ///          specify a combination of DiffTransConfigMapper::Options using bitwise OR: Ex. _options=DiffTransConfigMapper::rotate|ConfigMapper::strict
    ///          Options are:
    ///             'rotate': removes rigid rotation of the imported crystal, in a least-squares sense (i.e., yields a symmetric deformation tensor)
    ///             'robust': does not assume the imported structure might be ideal ('robust' is much slower for importing ideal structures, but if 'robust' is not
    ///                       set and a non-ideal structure is passed, this will be almost always be detected and robust methods will be used instead. Thus, 'robust'
    ///                       is slightly faster if imported Structures are *not* ideal)
    ///             'strict': prevents transformation into canonical form. Tries to preserve original orientation of imported structure if possible
    ///\endparblock
    ///
    ///\param _tol tolerance for mapping comparisons
    DiffTransConfigMapper(const PrimClex &_pclex,
                          double _lattice_weight,
                          double _max_volume_change = 0.5,
                          int _options = robust, // this should actually be a bitwise-OR of ConfigMapper::Options
                          double _tol = TOL);


    const PrimClex &primclex() const {
      return *m_pclex;
    }

    void set_primclex(PrimClex &_pclex) {
      m_pclex = &_pclex;
    }

    double lattice_weight() const {
      return m_lattice_weight;
    }

    void set_lattice_weight(double _lw) {
      m_lattice_weight = max(min(_lw, 1.0), 1e-9);
    }

    double min_va_frac() const {
      return m_min_va_frac;
    }

    void set_min_va_frac(double _min_va) {
      m_min_va_frac = max(_min_va, 0.);
    }

    double max_va_frac() const {
      return m_max_va_frac;
    }

    void set_max_va_frac(double _max_va) {
      m_max_va_frac = min(_max_va, 1.);
    }


    ///\brief imports structures specified by 'pos_path' into primclex() by finding optimal mapping
    ///       and then setting displacements and strain to zero (only the mapped occupation is preserved)
    ///
    DiffTransConfigMapperResult import_structure_occupation(const fs::path &pos_path) const;

    ///\brief imports structure specified by '_struc' into primclex()
    /// Not viable for diff_trans_configs
    ///ConfigMapperResult import_structure_occupation(const BasicStructure<Site> &_struc) const;

    ///\brief imports structure specified by 'pos_path' into primclex()
    ///\param hint_ptr[in]
    ///\parblock
    ///                provides a suggestion for which DiffTransConfiguration _struc should map onto
    ///                The hint is used to reduce search times, but may be used in the future
    ///                in combination with Option 'strict' to force mapping onto a particular configuration
    ///                or be used to provide user reports of the form "Suggested mapping: 0.372; Optimal mapping: 0.002"
    ///\endparblock
    ///
    DiffTransConfigMapperResult import_structure_occupation(const fs::path &pos_path,
                                                            const Kinetics::DiffTransConfiguration *hint_ptr) const;


    ///\brief imports structure specified by 'pos_path' into primclex() by finding optimal mapping
    ///       unlike import_structure_occupation, displacements and strain are preserved
    ///
    DiffTransConfigMapperResult import_structure(const fs::path &pos_path) const;

    //private:

    std::vector<BasicStructure<Site>> _get_structures(const fs::path &pos_path) const;

  private:

    const PrimClex *m_pclex;
    mutable std::map<Index, std::vector<Lattice> > m_superlat_map;
    double m_lattice_weight;
    double m_max_volume_change;
    double m_min_va_frac;
    double m_max_va_frac;
    bool m_robust_flag, m_strict_flag, m_rotate_flag;
    double m_tol;
    std::vector<std::pair<std::string, Index> > m_fixed_components;
    const std::vector<Lattice> &_lattices_of_vol(Index prim_vol) const;
  };



}

#endif
