#ifndef CASM_ClexBasisWriter
#define CASM_ClexBasisWriter

#include "casm/CASM_global_definitions.hh"
#include "casm/clex/ClexBasis.hh"
namespace CASM {
  class ClexBasis;
  class Structure;
  class PrimNeighborList;
  class ParamPackMixIn;
  class OrbitFunctionTraits;



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class ClexBasisWriter {
  public:

    /// \brief Construct ClexBasisWriter, collecting requisite DoF info from '_prim'
    ClexBasisWriter(Structure const &_prim, ParamPackMixIn const *parampack_mix_in = null_ptr);

    /// \brief Print clexulator
    template <typename OrbitType>
    void print_clexulator(std::string class_name,
                          ClexBasis const &clex,
                          std::vector<OrbitType > const &_tree,
                          PrimNeighborList &_nlist,
                          std::vector<UnitCellCoord> const &_flower_pivots,
                          std::ostream &stream,
                          double xtal_tol);

  private:
    std::vector<std::unique_ptr<FunctionVisitor> > const &_function_label_visitors() const;
    std::vector<std::unique_ptr<OrbitFunctionTraits> > const &_orbit_func_traits() const;

    /// \brief Print ClexParamPack specialization
    template <typename OrbitType>
    void print_param_pack(std::string class_name,
                          ClexBasis const &clex,
                          std::vector<OrbitType > const &_tree,
                          PrimNeighborList &_nlist,
                          std::vector<UnitCellCoord> const &_flower_pivots,
                          std::ostream &stream,
                          double xtal_tol);

  };



  namespace ClexBasisWriter_impl {


    std::string clexulator_member_definitions(std::string const &class_name,
                                              ClexBasis const &clex,
                                              std::vector<std::unique_ptr<OrbitFunctionTraits> > const &orbit_func_traits,
                                              std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_private_method_definitions(std::string const &class_name,
                                                      ClexBasis const &clex,
                                                      std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_public_method_definitions(std::string const &class_name,
                                                     ClexBasis const &clex,
                                                     std::string const &indent);
    //*******************************************************************************************

    template <typename OrbitType>
    std::tuple<std::string, std::string> clexulator_orbit_function_strings(std::string const &class_name,
                                                                           ClexBasis::BSetOrbit const &_bset_orbit,
                                                                           OrbitType const &_clust_orbit,
                                                                           std::function<std::string(Index, Index)> method_namer,
                                                                           PrimNeighborList &_nlist,
                                                                           std::vector<std::unique_ptr<FunctionVisitor> > const &labelers,
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
                                                                            PrimNeighborList &_nlist,
                                                                            std::vector<std::unique_ptr<FunctionVisitor> > const &labelers,
                                                                            std::string const &indent);

    //*******************************************************************************************

    template <typename OrbitType>
    std::tuple<std::string, std::string> clexulator_dflower_function_strings(std::string const &class_name,
                                                                             ClexBasis::BSetOrbit const &_bset_orbit,
                                                                             ClexBasis::BSetOrbit const &_site_bases,
                                                                             OrbitType const &_clust_orbit,
                                                                             std::function<std::string(Index, Index)> method_namer,
                                                                             PrimNeighborList &_nlist,
                                                                             std::vector<std::unique_ptr<FunctionVisitor> > const &labelers,
                                                                             FunctionVisitor const &prefactor_labeler,
                                                                             std::string const &indent);
    //*******************************************************************************************

    std::string clexulator_interface_implementation(std::string const &class_name,
                                                    ClexBasis const &clex,
                                                    std::string const &indent);

    //*******************************************************************************************

    template <typename OrbitType>
    std::string clexulator_constructor_implementation(std::string const &class_name,
                                                      ClexBasis const &clex,
                                                      std::vector<OrbitType> const &_tree,
                                                      PrimNeighborList &_nlist,
                                                      std::vector<std::string> const &orbit_method_names,
                                                      std::vector< std::vector<std::string> > const &flower_method_names,
                                                      std::vector< std::vector<std::string> > const &dflower_method_names,
                                                      std::string const &indent);

    //*******************************************************************************************

    template <typename OrbitType>
    std::string clexulator_point_prepare_implementation(std::string const &class_name,
                                                        ClexBasis const &clex,
                                                        std::vector<OrbitType> const &_tree,
                                                        PrimNeighborList &_nlist,
                                                        std::vector<std::unique_ptr<OrbitFunctionTraits> > const &orbit_func_traits,
                                                        //Something that contains info about DoF requirements
                                                        std::string const &indent);
    //*******************************************************************************************

    template <typename OrbitType>
    std::string clexulator_global_prepare_implementation(std::string const &class_name,
                                                         ClexBasis const &clex,
                                                         std::vector<OrbitType> const &_tree,
                                                         PrimNeighborList &_nlist,
                                                         std::vector<std::unique_ptr<OrbitFunctionTraits> > const &orbit_func_traits,
                                                         //Something that contains info about DoF requirements
                                                         std::string const &indent);

    //*******************************************************************************************

    // Divide by multiplicity. Same result as evaluating correlations via orbitree.
    template<typename OrbitType>
    std::vector<std::string> orbit_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                        OrbitType const &_clust_orbit,
                                                        PrimNeighborList &_nlist,
                                                        std::vector<std::unique_ptr<FunctionVisitor> > const &labelers);
    //*******************************************************************************************
    /// nlist_index is the index of the basis site in the neighbor list
    template<typename OrbitType>
    std::map< UnitCell, std::vector< std::string > > flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                                                 std::function<BasisSet(BasisSet const &)> _bset_transform,
                                                                                 OrbitType const &_clust_orbit,
                                                                                 PrimNeighborList &_nlist,
                                                                                 std::vector<std::unique_ptr<FunctionVisitor> > const &labelers,
                                                                                 Index sublat_index);

    //*******************************************************************************************

    template<typename OrbitType>
    void print_proto_clust_funcs(ClexBasis const &clex,
                                 std::ostream &out,
                                 std::vector<OrbitType > const &_tree);

    //*******************************************************************************************

    template<typename UCCIterType, typename IntegralClusterSymCompareType>
    std::map<UnitCellCoord, std::set<UnitCellCoord> > unique_ucc(UCCIterType begin,
                                                                 UCCIterType end,
                                                                 IntegralClusterSymCompareType const &sym_compare);
  }
}

#endif
