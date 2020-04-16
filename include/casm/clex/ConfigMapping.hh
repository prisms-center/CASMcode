#ifndef CASM_ConfigMapping
#define CASM_ConfigMapping

#include "casm/global/definitions.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
  namespace xtal {
    class Lattice;
    class Site;
    class BasicStructure;
    struct MappingNode;
    class SimpleStructure;
    class SimpleStrucMapCalculator;
    class StrucMapper;
  }
  using xtal::Lattice;
  using xtal::Site;
  using xtal::BasicStructure;
  using xtal::MappingNode;
  using xtal::SimpleStructure;
  using xtal::SimpleStrucMapCalculator;
  using xtal::StrucMapper;

  class Supercell;
  class PermuteIterator;
  class PrimClex;
  class Configuration;
  class ConfigDoF;
  class SupercellSymInfo;

  namespace ConfigMapping {

    /// \brief Struct with optional parameters for Config Mapping
    /// Specifies default parameters for all values, in order to simplify
    /// parsing from JSON

    struct Settings {
      Settings(double _lattice_weight = 0.5,
               bool _ideal = false,
               bool _strict = false,
               bool _robust = false,
               bool _primitive_only = false,
               bool _fix_volume = false,
               bool _fix_lattice = false,
               Index _k_best = 1,
               std::vector<std::string> _forced_lattices = {},
               std::string _filter = "",
               double _cost_tol = CASM::TOL,
               double _min_va_frac = 0.,
               double _max_va_frac = 0.5,
               double _max_vol_change = 0.3) :
        lattice_weight(_lattice_weight),
        ideal(_ideal),
        strict(_strict),
        robust(_robust),
        primitive_only(_primitive_only),
        fix_volume(_fix_volume),
        fix_lattice(_fix_lattice),
        k_best(_k_best),
        forced_lattices(_forced_lattices),
        filter(_filter),
        cost_tol(_cost_tol),
        min_va_frac(_min_va_frac),
        max_va_frac(_max_va_frac),
        max_vol_change(_max_vol_change) {}

      int options()const {
        int opt = 0;
        if(robust)
          opt |= StrucMapper::robust;
        return opt;
      }

      void set_default() {
        *this = Settings();
      }

      /// lattice_weight specifies the cost function in terms of lattice deformation cost and
      /// atomic deformation cost (i.e., atomic displacement)
      /// cost = lattice_weight*lattice_cost + (1l-lattice_weight)*atomic_displacement_cost
      double lattice_weight;

      /// True if child structure's lattice should be assumed to be ideal integer supercell of parent structure
      bool ideal;

      /// True invokes post-processing step to find a symmetry operation of the parent structure that preserves
      /// setting of child structure as much as possible after mapping
      bool strict;

      /// True invokes additional checks which determine whether there are any other mappings that are distinct
      /// from the best mapping, but have the same cost
      bool robust;

      /// If true, non-primitive configurations are only inserted in the database in the form of their primitive form
      /// The primitive form is always inserted into the database, regardless of setting value
      bool primitive_only;

      /// If true, search for potential mappings will be constrained to the supercell volume of the starting config
      /// (update operations only)
      bool fix_volume;


      /// If true, search for potential mappings will be constrained to the exact supercell of the starting config
      /// (update operations only)
      bool fix_lattice;

      /// Specify the number, k, of k-best mappings to include in solution set (default is 1)
      Index k_best;

      /// List of superlattices of parent structure to consider when searching for mappings
      std::vector<std::string> forced_lattices;

      /// casm-query expression used to filter list of potential supercells of parent structure to search over
      std::string filter;

      /// Tolerance used to determine if two mappings have identical cost
      double cost_tol;

      /// minimum fraction of vacant sites, below this fraction a mapping will not be considered
      double min_va_frac;

      /// maximum fraction of vacant sites, above this fraction a mapping will not be considered
      double max_va_frac;

      /// constrains the search space by assuming a limit on allowed volume change
      /// only taken into account when non-interstitial vacancies are allowed in parent structure
      double max_vol_change;

    };
  }

  /// \brief Reorders the permutation and compounds the spatial isometry (rotation + translation) of _node with that of _it
  MappingNode copy_apply(PermuteIterator const &_it, MappingNode const &_node, bool transform_cost_mat = true);

  /// \brief Initializes configdof of Supercell '_scel' corresponding to an idealized child structure (encoded by _child_struc)
  /// _child_struc is assumed to have been idealized via structure-mapping or to be the result of converting a configuration to
  /// a SimpleStructure. result.second gives list of properties that were utilized in the course of building the configdof
  std::pair<ConfigDoF, std::set<std::string> > to_configdof(SimpleStructure const &_child_struc, Supercell const  &_scel);

  class PrimStrucMapCalculator : public SimpleStrucMapCalculator {
  public:
    PrimStrucMapCalculator(BasicStructure const &_prim,
                           std::vector<SymOp> const &symgroup = {},
                           SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM);

  private:
    /// \brief Make an exact copy of the calculator (including any initialized members)
    StrucMapCalculatorInterface *_clone() const override {
      return new PrimStrucMapCalculator(*this);
    }

    BasicStructure m_prim;

  };



  /// Data structure holding results of ConfigMapper algorithm
  struct ConfigMapperResult {
    /// Specify degree to which the hinted configuration matches the imported structure:
    ///    None : unspecified/unknown
    ///    Derivate : same occupation, but other DoFs are different
    ///    Equivalent : same occupation and DoFs, but related to Hint by a SymOp of the ideal crystal
    ///    Identical : Exact same occupation and DoFs
    ///    NewOcc : Occupation has no relation to Hint (no statement about other DoFs, but they should be presumed different)
    ///    NewScel : Mapped configuration corresponds to different supercell than Hint
    enum class HintStatus { None, Derivative, Equivalent, Identical, NewOcc, NewScel};

    struct Individual {
      Individual(Configuration _config,
                 SimpleStructure _resolved_struc,
                 std::set<std::string> _dof_managed_properties,
                 HintStatus _hint_status = HintStatus::None,
                 double _hint_cost = xtal::StrucMapping::big_inf()) :
        config(std::move(_config)),
        resolved_struc(std::move(_resolved_struc)),
        dof_managed_properties(std::move(_dof_managed_properties)),
        hint_status(_hint_status),
        hint_cost(_hint_cost) {}


      Configuration config;

      //Child structure, converted to setting of reference structure (i.e., prim)
      //It is equivalent to initial child structure, but the following transformations have been performed:
      // - Atoms converted to molecules and ordered as in reference structure
      // - Vacancies explicitly inserted, if any were inferred from mapping
      // - Sites are translated such that average displacement of child sites relative to parent sites is 0
      // - Lattice and coordinates rotated
      // - Integer combinations of child lattice vectors (Lc) to match parent lattice vectors (Lp) such that
      //      Lc = U * Lp
      //   where U is symmetric right stretch tensor
      // - Permutation and rotation also applied to all properties of resolved_struc
      SimpleStructure resolved_struc;

      // list of properties that are handled by DoFs and are thus not considered properties
      std::set<std::string> dof_managed_properties;

      HintStatus hint_status;
      double hint_cost;
    };

    ConfigMapperResult() {}

    bool success() const {
      return !maps.empty();
    }

    Index n_optimal(double tol = TOL) const;

    /// Mapped structure, before applying lattice similarity and/or rotation to
    /// input structure.
    SimpleStructure structure;

    /// The configurations that the input structure mapped onto
    std::map<MappingNode, Individual> maps;


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

    using HintStatus = ConfigMapperResult::HintStatus;

    ///\brief Construct and initialize a ConfigMapper
    ///\param _pclex the PrimClex that describes the crystal template
    ///
    ///\param _strain_weight
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
    ConfigMapper(PrimClex const &_pclex,
                 ConfigMapping::Settings const &_settings,
                 double _tol = -1.);


    const PrimClex &primclex() const {
      return *m_pclex;
    }

    ConfigMapping::Settings const &settings() const {
      return m_settings;
    }

    void set_primclex(const PrimClex &_pclex) {
      m_pclex = &_pclex;
    }

    StrucMapper const &struc_mapper() const {
      return m_struc_mapper;
    }

    void add_allowed_lattices(std::vector<std::string> const &_lattice_names);

    void clear_allowed_lattices();

    //STEPS:
    //    0) [If Hint] Do SimpleStructure -> SimpleStructure(Config) mapping => HintMapping (Default, HintMapping.cost=inf())
    //    1) If HintMapping.cost>tol  Do SimpleStructure -> PrimClex mapping => ClexMapping (Default, ClexMapping.cost=inf())
    //    2) If HintMapping.cost<ClexMapping.cost, use HintMapping, else use ClexMapping => BestMapping
    //    3) Convert BestMapping to ConfigDoF
    //       [a] - BestMapping attributes that define ConfigDoF are mapped 'DoF', all others mapped 'property'
    //       [b] - 'property' attributes are subsumed into 'relaxation_properties' object
    //    4) Construct Configuration as ConfigDoF + relation_properties

    ///\brief imports structure specified by '_struc' into primclex()
    ///\param hint_ptr[in]
    ///\parblock
    ///                provides a suggestion for which Configuration _struc should map onto
    ///                The hint is used to reduce search times, but may be used in the future
    ///                in combination with Option 'strict' to force mapping onto a particular configuration
    ///                or be used to provide user reports of the form "Suggested mapping: 0.372; Optimal mapping: 0.002"
    ///\endparblock
    ///
    ConfigMapperResult import_structure(SimpleStructure const &_struc,
                                        Configuration const *hint_ptr = nullptr,
                                        std::vector<DoFKey> const &_hint_dofs = {"occ"}) const;

    ConfigMapperResult import_structure(SimpleStructure const &_struc,
                                        Index k,
                                        Configuration const *hint_ptr = nullptr,
                                        std::vector<DoFKey> const &_hint_dofs = {"occ"}) const;


  private:
    const PrimClex *m_pclex;
    ///Maps the supercell volume to a vector of Lattices with that volume
    StrucMapper m_struc_mapper;

    ConfigMapping::Settings m_settings;
  };


}

#endif
