#ifndef FUNCTIONVISITOR_HH
#define FUNCTIONVISITOR_HH

#include <iostream>
#include <sstream>  //John G 010413
#include <map>

#include "casm/container/Array.hh"

namespace CASM {

  class Function;
  class MonomialFunction;
  class PolynomialFunction;

  class OccupantFunction;

  /// Defines visitor pattern for abstract Function objects.
  /// FunctionVisitor has virtual methods for each type of derived Function object,
  /// which by default do nothing. Classes that derive from FunctionVisitor implement
  /// specific visitation behavior by overloading the virtual visit() methods, as needed.
  /// An abstract FunctionVisitor object operates on a Function tree by being passed to the
  /// Function::accept(FunctionVisitor*) method of the root of the tree.
  class FunctionVisitor {

  public:
    virtual ~FunctionVisitor() {};

    virtual std::string type_name() const = 0;

    virtual bool visit(OccupantFunction &host)const {
      return false;
    };
    virtual bool visit(MonomialFunction &host)const {
      return false;
    };
    virtual bool visit(PolynomialFunction &host)const {
      return false;
    };
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Selectively relabel all OccupantFunctions in a Function tree, using their
  /// set_formula() method. OccFuncLabeler is constructed using a template string,
  /// which should be of the form (substr1 + "%n" + substr2 + "%f" + substr3 + "%b" + substr4),
  /// where substr1, etc, are user-defined substrings. "%n" indicates the placement of
  /// the neighbor_list index, "%f" indicates the placement of the function index, and
  /// "%b" indicates the placement of the basis site index. The relative order of "%n", "%f", and "%b"
  /// is determined by the user.


  class OccFuncLabeler : public FunctionVisitor {
    Array<std::string> m_sub_strings;
    mutable std::stringstream m_ss;
  public:
    OccFuncLabeler(const std::string &_template);

    std::string type_name() const {
      return "OccFuncLabeler";
    };

    bool visit(OccupantFunction &host)const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class OccFuncBasisIndexer : public FunctionVisitor {
    int m_new_index;
  public:
    OccFuncBasisIndexer(int _new_index) : m_new_index(_new_index) {};

    std::string type_name() const {
      return "OccFuncBasisIndexer";
    };

    bool visit(OccupantFunction &host)const;
  };

}
#endif
