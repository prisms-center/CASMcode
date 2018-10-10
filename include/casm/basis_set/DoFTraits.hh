#ifndef CASM_DoFTraits
#define CASM_DoFTraits

#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/symmetry/OrbitDecl.hh"
#include "casm/clusterography/ClusterDecl.hh"

// remove these once implementation of derived classes gets moved out of this file
#include "casm/crystallography/Site.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {
  class jsonParser;
  class MasterSymGroup;
  class PrimNeighborList;
  class Structure;
  //namespace DoF_impl{
  //class OccupationDoFTraits;
  //}

  namespace DoFType {
    class Traits;

    Traits const &traits(std::string const &dof_key);

    BasicTraits const &basic_traits(std::string const &dof_key);

    //DoF_impl::OccupationDoFTraits occupation();

    /// \brief Collection of all the traits specific to a DoF type

    class Traits : public BasicTraits {
    public:
      Traits(std::string const &_type_name,
             std::vector<std::string> const &_std_var_names,
             DOF_DOMAIN _domain,
             DOF_MODE _mode) :
        BasicTraits(_type_name,
                    _std_var_names,
                    _domain,
                    _mode) {

      }

      /// \brief Allow destruction through base pointer
      virtual ~Traits() {}

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      virtual std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                         std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                         jsonParser const &_bspecs) const = 0;

      /// \brief Populate @param _in from JSON
      virtual void from_json(DoFSet &_in, jsonParser const &_json) const { }

      /// \brief Output @param _in to JSON
      virtual void to_json(DoFSet const &_out, jsonParser &_json) const { };

      /// \brief Generate a symmetry representation for the supporting vector space
      virtual Eigen::MatrixXd symop_to_matrix(SymOp const &op) const = 0;

      /*
      /// \brief Generate a symmetry representation for this DoF
      virtual SymGroupRepID generate_symrep(MasterSymGroup const &_group,
                                            Structure const &_prim,
                                            Index _nb) const = 0;

      */

      virtual std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors() const = 0;

      virtual std::vector<std::unique_ptr<FunctionVisitor> > clust_function_visitors() const = 0;

      virtual std::string site_basis_description(BasisSet site_bset, Site site) const = 0;

      virtual std::vector<std::pair<std::string, Index> > param_pack_allocation(std::vector<BasisSet> const &_bases) const = 0;

      virtual std::string clexulator_constructor_string(Structure const &_prim,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const = 0;

      virtual std::string clexulator_point_prepare_string(Structure const &_prim,
                                                          std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                          PrimNeighborList &_nlist,
                                                          std::vector<BasisSet> const &site_bases,
                                                          std::string const &indent) const = 0;

      virtual std::string clexulator_global_prepare_string(Structure const &_prim,
                                                           std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                           PrimNeighborList &_nlist,
                                                           std::vector<BasisSet> const &site_bases,
                                                           std::string const &indent) const = 0;

      virtual std::string clexulator_member_declarations_string(Structure const &_prim,
                                                                std::vector<BasisSet> const &site_bases,
                                                                std::string const &indent) const = 0;

      virtual std::string clexulator_private_method_declarations_string(Structure const &_prim,
                                                                        std::vector<BasisSet> const &site_bases,
                                                                        std::string const &indent) const = 0;

      virtual std::string clexulator_public_method_declarations_string(Structure const &_prim,
                                                                       std::vector<BasisSet> const &site_bases,
                                                                       std::string const &indent) const = 0;

      virtual std::string clexulator_private_method_implementations_string(Structure const &_prim,
                                                                           std::vector<BasisSet> const &site_bases,
                                                                           std::string const &indent) const = 0;

      virtual std::string clexulator_public_method_implementations_string(Structure const &_prim,
                                                                          std::vector<BasisSet> const &site_bases,
                                                                          std::string const &indent) const = 0;

      /// \brief non-virtual method to obtain copy through Traits pointer
      std::unique_ptr<Traits> clone() const {
        return std::unique_ptr<Traits>(static_cast<Traits *>(_clone()));
      }
    };



    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// \brief  Parsing dictionary for obtaining the correct BasicTraits given a name
    using TraitsDictionary = ParsingDictionary<BasicTraits>;

    /// This will eventually be managed by ProjectSettings
    //TraitsDictionary const &traits_dict();


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    inline
    Traits const &traits(std::string const &dof_key) {
      return static_cast<Traits const &>(DoF::traits(dof_key));
    }

    inline
    BasicTraits const &basic_traits(std::string const &dof_key) {
      return DoF::traits(dof_key);
    }

  }

  template<>
  DoFType::TraitsDictionary make_parsing_dictionary<DoF::BasicTraits>();


}
#endif
