#ifndef CASM_ConfigMapping
#define CASM_ConfigMapping

#include "casm/CASM_global_definitions.hh"
#include "casm/database/MappedProperties.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
  namespace xtal {
    class Lattice;
    class Site;
    template<typename T>
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


  namespace Completer {
    class ImportOption;
  }

  //TO IMPLEMENT
  /// \brief Find symop (as PermuteIterator) that gives the most 'faithful' equivalent mapping
  /// This means that
  ///   (1) the site permutation is as close to identity as possible (i.e., maximal character)
  ///   (2) ties at (1) are broken by ensuring _node.isometry is proper and close to zero rotation (i.e., maximal character*determinant)
  ///   (3) if (1) and (2) are ties, then we minimize _node.translation.norm()

  /// \brief Reorders the permutation and compounds the spatial isometry (rotation + translation) of _node with that of _it
  MappingNode copy_apply(PermuteIterator const &_it, MappingNode const &_node, bool transform_cost_mat = true);

  /// \brief Initializes configdof corresponding to a mapping (encoded by _node) of _child_struc onto _pclex
  ConfigDoF to_configdof(MappingNode const _node, SimpleStructure const &_child_struc, Supercell const  &_scel);
  //\TO IMPLEMENT

  class PrimStrucMapCalculator : public SimpleStrucMapCalculator {
  public:
    PrimStrucMapCalculator(BasicStructure<Site> const &_prim,
                           SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM);

    /// \brief Creates copy of _child_struc by applying isometry, lattice transformation, translation, and site permutation of _node
    SimpleStructure resolve_setting(MappingNode const &_node,
                                    SimpleStructure const &_child_struc) const override;

  private:
    /// \brief Make an exact copy of the calculator (including any initialized members)
    StrucMapCalculatorInterface *_clone() const override {
      return new PrimStrucMapCalculator(*this);
    }

    BasicStructure<Site> m_prim;

  };



  /// Data structure holding results of ConfigMapper algorithm
  struct ConfigMapperResult {
    enum class HintStatus { None, Derivative, Equivalent, Identical};

    struct MapData {
      MapData(std::string _name,
              HintStatus _hint_status = HintStatus::None) :
        name(std::move(_name)),
        hint_status(_hint_status) {}

      std::string name;
      MappedProperties props;
      HintStatus hint_status;
    };

    using Individual = std::pair<MapData, Configuration>;
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
                 double _strain_weight,
                 double _max_volume_change = 0.5,
                 int _options = StrucMapper::robust, // this should actually be a bitwise-OR of StrucMapper::Options
                 double _tol = -1.);


    const PrimClex &primclex() const {
      return *m_pclex;
    }

    bool strict() const {
      return struc_mapper().options() && StrucMapper::strict;
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

    void add_allowed_lattices(std::vector<std::string> const &_lattice_names);

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

    template <typename IterType, typename ScelType, typename OutputIter>
    void _to_configmaps(SimpleStructure const &child_struc,
                        IterType begin,
                        IterType end,
                        ScelType const &scel,
                        SupercellSymInfo const &sym_info,
                        OutputIter out) const;

    const PrimClex *m_pclex;
    ///Maps the supercell volume to a vector of Lattices with that volume
    StrucMapper m_struc_mapper;
  };


}

#endif
