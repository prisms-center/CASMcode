#ifndef CASM_ClexBasisWriter
#define CASM_ClexBasisWriter

#include "casm/global/definitions.hh"
#include "casm/clex/ClexBasis.hh"
namespace CASM {
  namespace xtal {
    class Structure;
  }
  class ClexBasis;
  class PrimNeighborList;
  class ParamPackMixIn;
  class OrbitFunctionTraits;

  using namespace xtal;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class ClexBasisWriter {
  public:

    /// \brief Construct ClexBasisWriter, collecting requisite DoF info from '_prim'
    ClexBasisWriter(Structure const &_prim, std::string const &parampack_type);

    /// \brief Construct ClexBasisWriter, collecting requisite DoF info from '_prim'
    ClexBasisWriter(Structure const &_prim, ParamPackMixIn const &parampack_mix_in);


    /// \brief Print clexulator
    template <typename OrbitType>
    void print_clexulator(std::string class_name,
                          ClexBasis const &clex,
                          std::vector<OrbitType > const &_tree,
                          PrimNeighborList &_nlist,
                          std::ostream &stream,
                          double xtal_tol);

  private:
    void _initialize(Structure const &_prim, ParamPackMixIn const &parampack_mix_in);

    std::vector<std::unique_ptr<FunctionVisitor> > const &_site_function_visitors() const {
      return m_site_visitors;
    }

    std::vector<std::unique_ptr<FunctionVisitor> > const &_clust_function_visitors() const {
      return m_clust_visitors;
    }

    std::vector<std::unique_ptr<OrbitFunctionTraits> > const &_orbit_func_traits() const {
      return m_orbit_func_traits;
    }

    /// \brief Print ClexParamPack specialization
    template <typename OrbitType>
    void print_param_pack(std::string class_name,
                          ClexBasis const &clex,
                          std::vector<OrbitType > const &_tree,
                          PrimNeighborList &_nlist,
                          std::ostream &stream,
                          std::string const &_indent)const;

    std::vector<std::unique_ptr<FunctionVisitor> > m_site_visitors;
    std::vector<std::unique_ptr<FunctionVisitor> > m_clust_visitors;
    std::vector<std::unique_ptr<OrbitFunctionTraits> > m_orbit_func_traits;
    notstd::cloneable_ptr<ParamPackMixIn> m_param_pack_mix_in;
  };



  namespace ClexBasisWriter_impl {


    std::string clexulator_member_declarations(std::string const &class_name,
                                               ClexBasis const &clex,
                                               ParamPackMixIn const &_param_pack_mix_in,
                                               std::vector<std::unique_ptr<OrbitFunctionTraits> > const &orbit_func_traits,
                                               std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                               std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_private_method_declarations(std::string const &class_name,
                                                       ClexBasis const &clex,
                                                       std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_public_method_declarations(std::string const &class_name,
                                                      ClexBasis const &clex,
                                                      std::string const &indent);
    //*******************************************************************************************

    template <typename OrbitType>
    std::tuple<std::string, std::string> clexulator_orbit_function_strings(std::string const &class_name,
                                                                           ClexBasis::BSetOrbit const &_bset_orbit,
                                                                           OrbitType const &_clust_orbit,
                                                                           std::function<std::string(Index, Index)> method_namer,
                                                                           PrimNeighborList &_nlist,
                                                                           std::vector<std::unique_ptr<FunctionVisitor> > const &visitors,
                                                                           std::string const &indent);
    //*******************************************************************************************

    /// \brief Print the flower function formulae for orbit @param _clust_orbit specified by BasisSet @param _bset_orbit
    /// The pivot of the flower is specified by @param _sublat_index

    /// The flower function of site @param _sublat_index and orbit @param _clust_orbit is obtained by summing the contributions of all
    /// cluster functions from @param _bset_orbit that 'touch' the site (b,i,j,k)=(sublat_index,0,0,0), including functions that are
    /// found by translations of equivalent clusters in @param _clust_orbit.
    /// Depending on the orbit periodicity (i.e., Orbit::sym_compare()), not all translations of the cluster that touch (sublat_index,0,0,0)
    /// are translationally equivalent. Thus, the result is the std::map that associates UnitCell (i.e, translation) to a set of formulae,
    /// (i.e., std::vector<std::string>), with one formula per function in _clust_orbit[i] (some or all formulae may evaluate to zero, if
    /// if @param _clust_orbit doesn't include site of type @param _sublat_index.

    /// @param _bset_transfrom is a function/functor that applies a transformation to each _bset_orbit[i].
    /// @param _nlist is the PrimNeighborList, used to index sites in the neighborhood
    /// @param _labelers is a set of FunctionVisitors that can be used to control formatting of the formulae

    template <typename OrbitType>
    std::tuple<std::string, std::string> clexulator_flower_function_strings(std::string const &class_name,
                                                                            ClexBasis::BSetOrbit const &_bset_orbit,
                                                                            OrbitType const &_clust_orbit,
                                                                            std::function<std::string(Index, Index)> method_namer,
                                                                            std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                                            PrimNeighborList &_nlist,
                                                                            std::vector<std::unique_ptr<FunctionVisitor> > const &visitor,
                                                                            std::string const &indent);

    //*******************************************************************************************

    template <typename OrbitType>
    std::tuple<std::string, std::string> clexulator_dflower_function_strings(std::string const &class_name,
                                                                             ClexBasis::BSetOrbit const &_bset_orbit,
                                                                             ClexBasis::BSetOrbit const &_site_bases,
                                                                             OrbitType const &_clust_orbit,
                                                                             std::function<std::string(Index, Index)> method_namer,
                                                                             std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                                             PrimNeighborList &_nlist,
                                                                             std::vector<std::unique_ptr<FunctionVisitor> > const &visitors,
                                                                             FunctionVisitor const &prefactor_labeler,
                                                                             std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_interface_declaration(std::string const &class_name,
                                                 ClexBasis const &clex,
                                                 ParamPackMixIn const &_param_pack_mix_in,
                                                 std::string const &indent);

    //*******************************************************************************************

    template <typename OrbitType>
    std::string clexulator_constructor_definition(std::string const &class_name,
                                                  ClexBasis const &clex,
                                                  std::vector<OrbitType> const &_tree,
                                                  std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                  PrimNeighborList &_nlist,
                                                  ParamPackMixIn const &_param_pack_mix_in,
                                                  std::vector<std::string> const &orbit_method_names,
                                                  std::vector< std::vector<std::string> > const &flower_method_names,
                                                  std::vector< std::vector<std::string> > const &dflower_method_names,
                                                  std::string const &indent);

    //*******************************************************************************************


    template <typename OrbitType>
    std::string clexulator_point_prepare_definition(std::string const &class_name,
                                                    ClexBasis const &clex,
                                                    std::vector<OrbitType> const &_tree,
                                                    std::vector<std::unique_ptr<OrbitFunctionTraits> > const &orbit_func_traits,
                                                    std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                    PrimNeighborList &_nlist,
                                                    std::string const &indent);
    //*******************************************************************************************

    template <typename OrbitType>
    std::string clexulator_global_prepare_definition(std::string const &class_name,
                                                     ClexBasis const &clex,
                                                     std::vector<OrbitType> const &_tree,
                                                     std::vector<std::unique_ptr<OrbitFunctionTraits> > const &orbit_func_traits,
                                                     std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                     PrimNeighborList &_nlist,
                                                     std::string const &indent);

    //*******************************************************************************************

    // Divide by multiplicity. Same result as evaluating correlations via orbitree.
    template<typename OrbitType>
    std::vector<std::string> orbit_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                        OrbitType const &_clust_orbit,
                                                        PrimNeighborList &_nlist,
                                                        std::vector<std::unique_ptr<FunctionVisitor> > const &visitors);
    //*******************************************************************************************
    /// nlist_index is the index of the basis site in the neighbor list
    template<typename OrbitType>
    std::vector<std::string>  flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                          std::function<BasisSet(BasisSet const &)> _bset_transform,
                                                          OrbitType const &_clust_orbit,
                                                          std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                          PrimNeighborList &_nlist,
                                                          std::vector<std::unique_ptr<FunctionVisitor> > const &_visitors,
                                                          UnitCellCoord const &_nbor);

    //*******************************************************************************************

    template<typename OrbitType>
    void print_proto_clust_funcs(ClexBasis const &clex,
                                 std::ostream &out,
                                 std::vector<OrbitType > const &_tree);

    //*******************************************************************************************

    template<typename OrbitIterType>
    std::map<UnitCellCoord, std::set<UnitCellCoord> > dependency_neighborhood(OrbitIterType begin,
                                                                              OrbitIterType end);

    //*******************************************************************************************

    template<typename UCCIterType, typename IntegralClusterSymCompareType>
    std::set<UnitCellCoord>  equiv_ucc(UCCIterType begin,
                                       UCCIterType end,
                                       UnitCellCoord const &_pivot,
                                       Structure const &_prim,
                                       IntegralClusterSymCompareType const &_sym_compare);
  }
}


#include "casm/clex/ClexBasisWriter_impl.hh"

#endif
