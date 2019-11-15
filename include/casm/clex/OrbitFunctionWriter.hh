#ifndef CASM_OrbitFunctionTraits
#define CASM_OrbitFunctionTraits

#include "casm/global/definitions.hh"

namespace CASM {

  /// \brief virtual base class for printing orbit functions of type specified by implementation.
  class OrbitFunctionTraits {
  public:
    static string class_desc() {
      return "Orbit Function Traits";
    }


    // Constructor? OrbitFunctionTraits()
    //

    virtual std::string name() const = 0;

    virtual void print_param_pack_initilialization() const {}

    virtual void print_to_point_prepare() const {}

    virtual void print_to_global_prepare() const {}

    virtual void print_typedefs(std::ostream &out,
                                std::string const &class_name,
                                std::string const &indent)const {}


    virtual void print_eval_table_definitions(std::ostream &out,
                                              std::string const &class_name,
                                              ClexBasis const &clex,
                                              std::string const &indent)const {}

  private:
    std::string m_name;
    std::vector<std::string> m_signature;
    std::vector<std::string> m_arg_names;
    std::string m_short_desc;
    std::string m_long_desc;

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}

#endif
