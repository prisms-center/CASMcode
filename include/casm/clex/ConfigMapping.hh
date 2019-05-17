#ifndef CASM_ConfigMapping
#define CASM_ConfigMapping

#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/StrucMapping.hh"


namespace CASM {
  class Supercell;
  class Lattice;
  class SymGroup;
  class PermuteIterator;
  class PrimClex;
  class Configuration;
  class ConfigDoF;
  class Structure;

  namespace Completer {
    class ImportOption;
  }

  class PrimStrucMapCalculator : public StrucMapCalculatorInterface {
  public:
    PrimStrucMapCalculator(BasicStructure<Site> const &_prim,
                           StrucMapping::SpeciesMode _species_mode = StrucMapping::ATOM);


    SimpleStructure::Info const &info(SimpleStructure const &_struc) const override {
      return m_species_mode == StrucMapping::ATOM ?  _struc.atom_info : _struc.mol_info;
    }

    Index max_n_va() const override;

    std::vector<Eigen::Vector3d> get_translations(SimpleStructure const &child_struc,
                                                  LatticeNode const &lat_node) const override;

    bool populate_cost_mat(MappingNode &_node,
                           SimpleStructure const &child_struc) const override;

    void populate_displacement(MappingNode &_node,
                               SimpleStructure const &child_struc) const override;

    void populate_properties(MappingNode &_node,
                             const SimpleStructure &child_struc) const override;

    bool validate(MappingNode const &_node) const override {
      return true;
    }

  private:
    /// \brief Make an exact copy of the calculator (including any initialized members)
    StrucMapCalculatorInterface *_clone() const override {
      return new PrimStrucMapCalculator(*this);
    }

    BasicStructure<Site> m_prim;
    StrucMapping::SpeciesMode m_species_mode;
    std::vector<std::pair<std::string, Index> > m_fixed_components;

  };

  /// Data structure holding results of ConfigMapper algorithm
  struct ConfigMapperResult {

    ConfigMapperResult() {}

    bool success() const {
      return !StrucMapping::is_inf(result.cost);
    }

    MappingNode result;

    /// Output structure, after applying lattice similarity and/or rotation to
    /// input structure.
    SimpleStructure structure;

    /// The configuration the input structure was mapped onto
    std::unique_ptr<Configuration> config;

    /// relaxation_properties is populated by relaxation properties:
    ///
    /// - 'lattice_deformation': lattice mapping score
    /// - 'basis_deformation': atomic mapping score
    /// - 'volume_relaxation': V/V_ideal
    /// - 'relaxation_deformation': 3x3 tensor describing cell relaxation
    /// - 'relaxation_displacement': Nx3 matrix describing basis displacements
    /// - 'relaxed_energy': the energy of the relaxed configuration
    jsonParser relaxation_properties;

    /// Failure message if could not map to prim
    std::string fail_msg;

  };


  /// A class for mapping an arbitrary crystal structure as a configuration of a crystal template
  /// as described by a PrimClex.  ConfigMapper manages options for the mapping algorithm and mapping cost function
  /// It also caches some information about supercell lattices so that batch imports are more efficient
  ///
  /// \ingroup Configuration
  class ConfigMapper {
  public:

    enum NullInitializer {null_initializer};
    enum Options {none = 0,
                  rotate = (1u << 0),
                  strict = (1u << 1),
                  robust = (1u << 2)
                 };

    ///\brief Default construction not allowed -- this constructor provides an override
    ConfigMapper(NullInitializer) :
      m_pclex(nullptr),
      m_struc_mapper(StrucMapper::null_initializer) {
    }

    ///\brief Construct and initialize a ConfigMapper
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
    ///          specify a combination of ConfigMapper::Options using bitwise OR: Ex. _options=ConfigMapper::rotate|ConfigMapper::strict
    ///          Options are:
    ///             'rotate': removes rigid rotation of the imported crystal, in a least-squares sense (i.e., yields a symmetric deformation tensor)
    ///             'robust': does not assume the imported structure might be ideal ('robust' is much slower for importing ideal structures, but if 'robust' is not
    ///                       set and a non-ideal structure is passed, this will be almost always be detected and robust methods will be used instead. Thus, 'robust'
    ///                       is slightly faster if imported Structures are *not* ideal)
    ///             'strict': prevents transformation into canonical form. Tries to preserve original orientation of imported structure if possible
    ///\endparblock
    ///
    ///\param _tol tolerance for mapping comparisons (default is _pclex.crystallography_tol())
    ConfigMapper(const PrimClex &_pclex,
                 double _lattice_weight,
                 double _max_volume_change = 0.5,
                 int _options = robust, // this should actually be a bitwise-OR of ConfigMapper::Options
                 double _tol = -1.);


    const PrimClex &primclex() const {
      return *m_pclex;
    }

    void set_primclex(const PrimClex &_pclex) {
      m_pclex = &_pclex;
    }

    StrucMapper &struc_mapper() {
      return m_struc_mapper;
    }

    StrucMapper const &struc_mapper() const {
      return m_struc_mapper;
    }

    void add_allowed_lattices(std::vector<std::string> const &_lattice_names) {

      throw std::runtime_error("ConfigMapper::add_allowed_lattices() not implemented");
    }

    void clear_allowed_lattices() {
      struc_mapper().clear_allowed_lattices();
    }
    //STEPS:
    //    0) [If Hint] Do SimpleStructure -> SimpleStructure(Config) mapping => HintMapping (Default, HintMapping.cost=inf())
    //    1) If HintMapping.cost>tol  Do SimpleStructure -> PrimClex mapping => ClexMapping (Default, ClexMapping.cost=inf())
    //    2) If HintMapping.cost<ClexMapping.cost, use HintMapping, else use ClexMapping => BestMapping
    //    3) Convert BestMapping to ConfigDoF
    //       [a] - BestMapping attributes that define ConfigDoF are mapped 'DoF', all others mapped 'property'
    //       [b] - 'property' attributes are subsumed into 'relaxation_properties' object
    //    4) Construct Configuration as ConfigDoF + relation_properties

    ///\brief imports structure specified by 'pos_path' into primclex() by finding optimal mapping
    ///       unlike import_structure_occupation, displacements and strain are preserved
    ///
    ConfigMapperResult import_structure(const std::string &pos_path) const;

    ///\brief imports structure specified by '_struc' into primclex() by finding optimal mapping
    ///       unlike import_structure_occupation, displacements and strain are preserved
    ///
    ConfigMapperResult import_structure(const SimpleStructure &_struc) const;

    ///\brief imports structure specified by '_struc' into primclex()
    ///\param hint_ptr[in]
    ///\parblock
    ///                provides a suggestion for which Configuration _struc should map onto
    ///                The hint is used to reduce search times, but may be used in the future
    ///                in combination with Option 'strict' to force mapping onto a particular configuration
    ///                or be used to provide user reports of the form "Suggested mapping: 0.372; Optimal mapping: 0.002"
    ///\endparblock
    ///
    ConfigMapperResult import_structure(const SimpleStructure &_struc,
                                        const Configuration *hint_ptr) const;


    std::vector<Index> occupation(const SimpleStructure &_struc, MappingNode const &_node) const {
      throw std::runtime_error("ConfigMapper::occupation() is not implemented");
      return {};
    }
  private:

    const PrimClex *m_pclex;
    ///Maps the supercell volume to a vector of Lattices with that volume
    StrucMapper m_struc_mapper;
  };


}

#endif
