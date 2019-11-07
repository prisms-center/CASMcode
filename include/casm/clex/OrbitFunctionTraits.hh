#ifndef CASM_OrbitFunctionTraits
#define CASM_OrbitFunctionTraits

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include "casm/CASM_global_definitions.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/basis_set/DoFDecl.hh"

namespace CASM {
  namespace xtal {
    class UnitCellCoord;
    class Structure;
  }
  using xtal::UnitCellCoord;
  using xtal::Structure;

  class BasisSet;
  class SymGroup;
  class ClexBasis;
  class PrimNeighborList;


  class ClexBasisBuilder {
  public:
    ClexBasisBuilder(std::string const &_name) : m_name(_name) {}

    virtual ~ClexBasisBuilder() {}

    std::string const &name() const {
      return m_name;
    }

    virtual void prepare(Structure const &_prim) {

    }

    virtual std::vector<DoFKey> filter_dof_types(std::vector<DoFKey> const &_dof_types) {
      return _dof_types;
    }

    virtual void pre_generate() {

    }

    virtual BasisSet build_proto(IntegralCluster const &_prototype,
                                 SymGroup const &_generating_group,
                                 std::vector<BasisSet const *> const &_arg_bases,
                                 Index max_poly_order,
                                 Index min_poly_order) const = 0;

    std::unique_ptr<ClexBasisBuilder> clone()const {
      return std::unique_ptr<ClexBasisBuilder>(_clone());
    }

  private:
    virtual ClexBasisBuilder *_clone()const = 0;
    std::string m_name;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// \brief virtual base class for printing orbit functions of type specified by implementation.
  class OrbitFunctionTraits {
  public:
    static std::string class_desc() {
      return "Orbit Function Traits";
    }


    OrbitFunctionTraits(std::string const &_name,
                        std::string const &_short_desc,
                        std::string const &_long_desc) :
      m_name(_name),
      m_short_desc(_short_desc),
      m_long_desc(_long_desc) {}


    std::string const &name() const {
      return m_name;
    }

    std::string const &short_desc() const {
      return m_short_desc;
    }

    std::string const &long_desc() const {
      return m_long_desc;
    }

    virtual std::unique_ptr<ClexBasisBuilder> basis_builder() const = 0;

    virtual void print_param_pack_initilialization() const {}

    virtual std::string clexulator_point_prepare_string(Structure const &_prim,
                                                        std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                        PrimNeighborList &_nlist,
                                                        std::string const &indent) const {
      return "";
    }

    virtual std::string clexulator_global_prepare_string(Structure const &_prim,
                                                         std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                         PrimNeighborList &_nlist,
                                                         std::string const &indent) const {
      return "";
    }
    virtual void print_typedefs(std::ostream &out,
                                std::string const &class_name,
                                std::string const &indent)const {}


    virtual void print_eval_table_declarations(std::ostream &out,
                                               std::string const &class_name,
                                               ClexBasis const &clex,
                                               std::string const &indent)const {}

  private:
    std::string m_name;
    //std::vector<std::string> m_signature;
    //std::vector<std::string> m_arg_names;
    std::string m_short_desc;
    std::string m_long_desc;

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class InvariantPolyBasisBuilder : public ClexBasisBuilder {
  public:
    InvariantPolyBasisBuilder(std::string const &_name) :
      ClexBasisBuilder(_name) {}

    BasisSet build_proto(IntegralCluster const &_prototype,
                         SymGroup const &_generating_group,
                         std::vector<BasisSet const *> const &_arg_bases,
                         Index max_poly_order,
                         Index min_poly_order) const override;

  private:
    ClexBasisBuilder *_clone()const override {
      return new InvariantPolyBasisBuilder(*this);
    }

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class InvariantPolyOrbitFunctionTraits : public OrbitFunctionTraits {
  public:

    InvariantPolyOrbitFunctionTraits() :
      OrbitFunctionTraits("invariant_poly",
                          "Invariant Polynomials",
                          "Orbit functions are constructed as invariant polynomials of site basis functions.") {}


    std::unique_ptr<ClexBasisBuilder> basis_builder() const override {
      return std::unique_ptr<ClexBasisBuilder>(new InvariantPolyBasisBuilder(name()));
    }


  };


}

#endif
