#ifndef CASM_DiffTransConfigMapping
#define CASM_DiffTransConfigMapping

#include <vector>
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/clex/ConfigMapping.hh"
namespace CASM {
  namespace xtal {
    class Site;
    class UnitCellCoord;
    template<typename CoordType>
    class BasicStructure;
  }
  using xtal::Site;
  using xtal::UnitCellCoord;
  using xtal::BasicStructure;

  class Supercell;
  class PrimClex;

  namespace ConfigMapping {
    struct Settings;
  }

  namespace Kinetics {
    class DiffTransConfiguration;
    class DiffusionTransformation;
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

    double kra;

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
  namespace Kinetics {
    class DiffTransConfigMapper {
    public:
      enum NullInitializer {null_initializer};
      ///\brief Default construction not allowed -- this constructor provides an override
      DiffTransConfigMapper(NullInitializer) :
        m_pclex(nullptr),
        m_tol(CASM::TOL) {
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
      ///          specify a combination of StrucMapper::Options using bitwise OR: Ex. _options=StrucMapper::robust|StrucMapper::strict
      ///          Options are:
      ///             'strict': prevents transformation into canonical form. Tries to preserve original orientation of imported structure if possible
      ///             'robust': does not assume the imported structure might be ideal ('robust' is much slower for importing ideal structures, but if 'robust' is not
      ///                       set and a non-ideal structure is passed, this will be almost always be detected and robust methods will be used instead. Thus, 'robust'
      ///                       is slightly faster if imported Structures are *not* ideal)
      ///             'sym_strain': only calculates contribution to lattice mapping score due to symmetry-breaking strains
      ///             'sym_basis': only calculates contribution to basis mapping score due to symmetry-breaking displacements
      ///\endparblock
      ///
      ///\param _tol tolerance for mapping comparisons
      DiffTransConfigMapper(const PrimClex &_pclex,
                            ConfigMapping::Settings const &_settings,
                            double _tol = TOL);


      const PrimClex &primclex() const {
        return *m_pclex;
      }

      void set_primclex(PrimClex &_pclex) {
        m_pclex = &_pclex;
      }

      ConfigMapping::Settings const &settings()const {
        return m_settings;
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

    private:

      /// \brief loads structures from folder of poscars or compound properties.calc.json
      std::vector<BasicStructure<Site> > _get_structures(const fs::path &pos_path) const;

      /// \brief Helper function that creates a diff_trans from the information of which atoms are moving
      Kinetics::DiffusionTransformation _make_hop(BasicStructure<Site> const &from_struc,
                                                  std::vector<UnitCellCoord> const &from_coords,
                                                  std::vector<UnitCellCoord> const &to_coords,
                                                  std::set<UnitCellCoord> const &vacancy_from,
                                                  std::set<UnitCellCoord> const &vacancy_to,
                                                  std::vector<Index> const &moving_atoms) const;

      /// \brief Adjusts the diff_trans to have a direction that is the shortest distance across a periodic boundary
      Kinetics::DiffusionTransformation _shortest_hop(const Kinetics::DiffusionTransformation &diff_trans, const Supercell &scel) const;

      /// \brief transforms a site into a unit cell coord for ease of diff_trans construction
      UnitCellCoord _site_to_uccoord(const Site &site, const PrimClex &pclex, double tol) const;

      /// Rotates,Strains  and translates the structures being imported/updated into the frame of reference of the prim
      void _precondition_from_and_to(const Eigen::Matrix3d &cart_op, const Eigen::Matrix3d &strain, const Eigen::Vector3d &trans, BasicStructure<Site> &from, BasicStructure<Site> &to) const;

      /// Takes an aligned  set of initial and final position and determines which atoms are moving
      std::vector<Index> _analyze_atoms_ideal(BasicStructure<Site> const &from,
                                              BasicStructure<Site> const &to,
                                              const Supercell &scel,
                                              double uccoord_tol,
                                              std::vector<UnitCellCoord> &from_uccoords,
                                              std::vector<UnitCellCoord> &to_uccoords,
                                              std::set<UnitCellCoord> &vacancy_from,
                                              std::set<UnitCellCoord> &vacancy_to) const;

      const PrimClex *m_pclex;
      ConfigMapping::Settings m_settings;
      double m_tol;
    };
  }


}

#endif
