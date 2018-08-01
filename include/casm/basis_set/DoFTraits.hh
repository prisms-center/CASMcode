#ifndef CASM_DoFTraits
#define CASM_DoFTraits

#include "casm/basis_set/DoF.hh"
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

  namespace DoF_impl {
    class Traits;
  }

  namespace DoFType {
    DoF_impl::Traits const &traits(std::string const &dof_key);

    DoF_impl::BasicTraits const &basic_traits(std::string const &dof_key);

    notstd::cloneable_ptr<typename DoF_impl::BasicTraits> occupation();
  }

  namespace DoF_impl {
    /// \brief Collection of all the traits specific to a DoF type

    class Traits : public BasicTraits {
    public:
      Traits(std::string const &_type_name,
             DOF_DOMAIN _domain,
             DOF_MODE _mode) :
        BasicTraits(_type_name,
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
      virtual void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const = 0;

      /// \brief Output @param _in to JSON
      virtual void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const = 0;

      /// \brief Generate a symmetry representation for this DoF
      virtual SymGroupRepID generate_symrep(MasterSymGroup const &_group,
                                            Structure const &_prim,
                                            Index _nb) const = 0;


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
    class OccupationDoFTraits : public Traits {
    public:
      OccupationDoFTraits():
        Traits("occ",
               DISCRETE,
               LOCAL) {
      }

      /// \brief Populate @param _in from JSON
      void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const override {
        throw std::runtime_error("OccupationDoFTraits::from_json not implemented!");
      }

      /// \brief Output @param _in to JSON
      void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const override {
        throw std::runtime_error("OccupationDoFTraits::to_json not implemented!");
      }

      /// \brief Generate a symmetry representation for this DoF
      SymGroupRepID generate_symrep(MasterSymGroup const &_group,
                                    Structure const &_prim,
                                    Index _nb) const override {
        throw std::runtime_error("OccupationDoFTraits::generate_symrep not implemented!");
      }

      std::string site_basis_description(BasisSet site_bset, Site site) const override;

      std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors() const override;

      std::vector<std::unique_ptr<FunctionVisitor> > clust_function_visitors() const override;

      std::vector<std::pair<std::string, Index> > param_pack_allocation(std::vector<BasisSet> const &_bases) const override;

      std::string clexulator_constructor_string(Structure const &_prim,
                                                std::vector<BasisSet> const &site_bases,
                                                std::string const &indent) const override;

      std::string clexulator_point_prepare_string(Structure const &_prim,
                                                  std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                  PrimNeighborList &_nlist,
                                                  std::vector<BasisSet> const &site_bases,
                                                  std::string const &indent) const override;

      std::string clexulator_global_prepare_string(Structure const &_prim,
                                                   std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                   PrimNeighborList &_nlist,
                                                   std::vector<BasisSet> const &site_bases,
                                                   std::string const &indent) const override;

      std::string clexulator_member_declarations_string(Structure const &_prim,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const override;

      std::string clexulator_private_method_declarations_string(Structure const &_prim,
                                                                std::vector<BasisSet> const &site_bases,
                                                                std::string const &indent) const override;

      std::string clexulator_public_method_declarations_string(Structure const &_prim,
                                                               std::vector<BasisSet> const &site_bases,
                                                               std::string const &indent) const override {
        //todo
        return std::string();
      }

      std::string clexulator_private_method_implementations_string(Structure const &_prim,
                                                                   std::vector<BasisSet> const &site_bases,
                                                                   std::string const &indent) const override {
        // todo
        return std::string();
      }

      std::string clexulator_public_method_implementations_string(Structure const &_prim,
                                                                  std::vector<BasisSet> const &site_bases,
                                                                  std::string const &indent) const override {
        // todo
        return std::string();
      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;
    protected:
      BasicTraits *_clone() const override;
    };


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class DisplacementDoFTraits : public Traits {
    public:
      DisplacementDoFTraits():
        Traits("disp",
               CONTINUOUS,
               LOCAL) {
      }

      /// \brief Populate @param _in from JSON
      void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const override {
        throw std::runtime_error("DisplacementDoFTraits::from_json not implemented!");
      }

      /// \brief Output @param _in to JSON
      void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const override {
        throw std::runtime_error("DisplacementDoFTraits::to_json not implemented!");
      }

      /// \brief Generate a symmetry representation for this DoF
      SymGroupRepID generate_symrep(MasterSymGroup const &_group,
                                    Structure const &_prim,
                                    Index _nb) const override {
        throw std::runtime_error("DisplacementDoFTraits::generate_symrep not implemented!");
      }

      std::string site_basis_description(BasisSet site_bset, Site site) const override;

      std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors() const override;

      std::vector<std::unique_ptr<FunctionVisitor> > clust_function_visitors() const override;

      std::vector<std::pair<std::string, Index> > param_pack_allocation(std::vector<BasisSet> const &_bases) const override;

      std::string clexulator_constructor_string(Structure const &_prim,
                                                std::vector<BasisSet> const &site_bases,
                                                std::string const &indent) const override;

      std::string clexulator_point_prepare_string(Structure const &_prim,
                                                  std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                  PrimNeighborList &_nlist,
                                                  std::vector<BasisSet> const &site_bases,
                                                  std::string const &indent) const override;

      std::string clexulator_global_prepare_string(Structure const &_prim,
                                                   std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                   PrimNeighborList &_nlist,
                                                   std::vector<BasisSet> const &site_bases,
                                                   std::string const &indent) const override;

      std::string clexulator_member_declarations_string(Structure const &_prim,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const override;

      std::string clexulator_private_method_declarations_string(Structure const &_prim,
                                                                std::vector<BasisSet> const &site_bases,
                                                                std::string const &indent) const override;

      std::string clexulator_public_method_declarations_string(Structure const &_prim,
                                                               std::vector<BasisSet> const &site_bases,
                                                               std::string const &indent) const override {
        //todo
        return std::string();
      }

      std::string clexulator_private_method_implementations_string(Structure const &_prim,
                                                                   std::vector<BasisSet> const &site_bases,
                                                                   std::string const &indent) const override {
        // todo
        return std::string();
      }

      std::string clexulator_public_method_implementations_string(Structure const &_prim,
                                                                  std::vector<BasisSet> const &site_bases,
                                                                  std::string const &indent) const override {
        // todo
        return std::string();
      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;
    protected:
      BasicTraits *_clone() const override;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    class StrainDoFTraits : public Traits {
    public:
      StrainDoFTraits():
        Traits("strain",
               CONTINUOUS,
               GLOBAL) {
      }

      /// \brief Populate @param _in from JSON
      void from_json(std::vector<ContinuousDoF> &_in, jsonParser const &_json) const override {
        throw std::runtime_error("StrainDoFTraits::from_json not implemented!");
      }

      /// \brief Output @param _in to JSON
      void to_json(std::vector<ContinuousDoF> const &_out, jsonParser &_json) const override {
        throw std::runtime_error("StrainDoFTraits::to_json not implemented!");
      }

      /// \brief Generate a symmetry representation for this DoF
      SymGroupRepID generate_symrep(MasterSymGroup const &_group,
                                    Structure const &_prim) const override {
        throw std::runtime_error("StrainDoFTraits::generate_symrep not implemented!");
      }

      std::string site_basis_description(BasisSet site_bset, Site site) const override;

      std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors() const override;

      std::vector<std::unique_ptr<FunctionVisitor> > clust_function_visitors() const override;

      std::vector<std::pair<std::string, Index> > param_pack_allocation(std::vector<BasisSet> const &_bases) const override;

      std::string clexulator_constructor_string(Structure const &_prim,
                                                std::vector<BasisSet> const &site_bases,
                                                std::string const &indent) const override;

      std::string clexulator_point_prepare_string(Structure const &_prim,
                                                  std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                  PrimNeighborList &_nlist,
                                                  std::vector<BasisSet> const &site_bases,
                                                  std::string const &indent) const override;

      std::string clexulator_global_prepare_string(Structure const &_prim,
                                                   std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                   PrimNeighborList &_nlist,
                                                   std::vector<BasisSet> const &site_bases,
                                                   std::string const &indent) const override;

      std::string clexulator_member_declarations_string(Structure const &_prim,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const override;

      std::string clexulator_private_method_declarations_string(Structure const &_prim,
                                                                std::vector<BasisSet> const &site_bases,
                                                                std::string const &indent) const override;

      std::string clexulator_public_method_declarations_string(Structure const &_prim,
                                                               std::vector<BasisSet> const &site_bases,
                                                               std::string const &indent) const override {
        //todo
        return std::string();
      }

      std::string clexulator_private_method_implementations_string(Structure const &_prim,
                                                                   std::vector<BasisSet> const &site_bases,
                                                                   std::string const &indent) const override {
        // todo
        return std::string();
      }

      std::string clexulator_public_method_implementations_string(Structure const &_prim,
                                                                  std::vector<BasisSet> const &site_bases,
                                                                  std::string const &indent) const override {
        // todo
        return std::string();
      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;
    protected:
      BasicTraits *_clone() const override;
    };


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// \brief  Parsing dictionary for obtaining the correct BasicTraits given a name
    using TraitsDictionary = ParsingDictionary<BasicTraits>;

    /// This will eventually be managed by ProjectSettings
    TraitsDictionary const &traits_dict();
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  namespace DoFType {
    inline
    DoF_impl::Traits const &traits(std::string const &dof_key) {
      return static_cast<DoF_impl::Traits const &>(DoF::traits(dof_key));
    }

    inline
    DoF_impl::BasicTraits const &basic_traits(std::string const &dof_key) {
      return DoF::traits(dof_key);
    }

    inline
    notstd::cloneable_ptr<typename DoF_impl::BasicTraits> occupation() {
      return DoF_impl::OccupationDoFTraits().clone();
    }

    inline
    notstd::cloneable_ptr<typename DoF_impl::BasicTraits> displacement() {
      return DoF_impl::DisplacementDoFTraits().clone();
    }

    inline
    notstd::cloneable_ptr<typename DoF_impl::BasicTraits> strain() {
      return DoF_impl::StrainDoFTraits().clone();
    }
  }

}
#endif
