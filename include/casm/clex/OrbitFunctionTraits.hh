#ifndef CASM_OrbitFunctionTraits
#define CASM_OrbitFunctionTraits

#include <vector>
#include <string>

#include "casm/CASM_global_definitions.hh"
#include "casm/clusterography/ClusterDecl.hh"

namespace CASM {

  class Structure;
  class BasisSet;
  class SymGroup;
  class ClexBasis;

  class ClexBasisBuilder {
  public:
    virtual ~ClexBasisBuilder() {}

    virtual void prepare(Structure const &_prim) {

    }

    virtual std::vector<ClexBasis::DoFKey> filter_dof_types(std::vector<ClexBasis::DoFKey> const &_dof_types) {
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

    virtual void print_to_point_prepare() const {}

    virtual void print_to_global_prepare() const {}

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
      return std::unique_ptr<ClexBasisBuilder>(new InvariantPolyBasisBuilder());
    }


  };


}

#endif
