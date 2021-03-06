#ifndef FUNCTIONVISITOR_HH
#define FUNCTIONVISITOR_HH

#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "casm/global/definitions.hh"

namespace CASM {
class BasisSet;

class Function;
class Variable;
class PolynomialFunction;
class OccupantFunction;

template <typename T>
class Array;

namespace FunctionLabel_impl {
std::pair<bool, std::string> make_label(Function const &host,
                                        BasisSet const *bset_ptr);
}

/// Defines visitor pattern for abstract Function objects.
/// FunctionVisitor has virtual methods for each type of derived Function
/// object, which by default do nothing. Classes that derive from
/// FunctionVisitor implement specific visitation behavior by overloading the
/// virtual visit() methods, as needed. An abstract FunctionVisitor object
/// operates on a Function tree by being passed to the
/// Function::accept(FunctionVisitor*) method of the root of the tree.
class FunctionVisitor {
 public:
  virtual ~FunctionVisitor() {}

  virtual std::string type_name() const = 0;

  virtual bool visit(Variable const &host, BasisSet const *bset_ptr) const;

  virtual bool visit(Variable &host, BasisSet const *bset_ptr) const;

  virtual bool visit(OccupantFunction const &host,
                     BasisSet const *bset_ptr) const;

  virtual bool visit(OccupantFunction &host, BasisSet const *bset_ptr) const;

  virtual bool visit(PolynomialFunction const &host,
                     BasisSet const *bset_ptr) const;

  virtual bool visit(PolynomialFunction &host, BasisSet const *bset_ptr) const;

 private:
  virtual bool _generic_visit(Function &host, BasisSet const *bset_ptr) const;

  virtual bool _generic_visit(Function const &host,
                              BasisSet const *bset_ptr) const;

  virtual bool _visit(const Array<Function *> &host_list,
                      BasisSet const *bset_ptr) const {
    return false;
  }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Selectively relabel all OccupantFunctions in a Function tree, using their
/// set_formula() method. OccFuncLabeler is constructed using a template string,
/// which should be of the form (substr1 + "%n" + substr2 + "%f" + substr3 +
/// "%b" + substr4), where substr1, etc, are user-defined substrings. "%n"
/// indicates the placement of the neighbor_list index, "%f" indicates the
/// placement of the function index, and
/// "%b" indicates the placement of the basis site index. The relative order of
/// "%n", "%f", and "%b" is determined by the user.

class OccFuncLabeler : public FunctionVisitor {
  std::vector<std::string> m_sub_strings;

 public:
  OccFuncLabeler(const std::string &_template);

  std::string type_name() const { return "OccFuncLabeler"; }

  bool visit(OccupantFunction &host, BasisSet const *bset_ptr) const;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class OccFuncBasisIndexer : public FunctionVisitor {
  int m_new_index;

 public:
  OccFuncBasisIndexer(int _new_index) : m_new_index(_new_index) {}

  std::string type_name() const { return "OccFuncBasisIndexer"; }

  bool visit(OccupantFunction &host, BasisSet const *bset_ptr) const;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// \brief Evaluates
class OccFuncEvaluator : public FunctionVisitor {
 public:
  OccFuncEvaluator(int state) : m_state(state) {}

  double value() { return m_value; }

  std::string type_name() const { return "OccFuncEvaluator"; }

  bool visit(OccupantFunction &host, BasisSet const *bset_ptr) const;

  bool visit(OccupantFunction const &host, BasisSet const *bset_ptr) const;

 private:
  bool _generic_visit(Function const &host, BasisSet const *bset_ptr) const {
    throw std::runtime_error(
        "OccFuncEvaluator can only be applied to OccupantFunction!");
    return false;
  }

  int m_state;
  mutable double m_value;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Selectively relabel all Variabless in a Function tree, using their
/// set_formula() method. VariableLabeler is constructed using a template
/// string, which should be of the form (substr1 + "%n" + substr2 + "%p" +
/// substr3 + "%s" + substr4), where substr1, etc, are user-defined substrings.
/// "%n" indicates the placement of the neighbor_list index, "%p" indicates
/// placement of the DoF-type prefix, and
/// "%s" indicates placement of the DoF-type suffix (all ContinuousDoFs are
/// assumed to report a type-name of the form "prefix_suffix", where 'prefix'
/// denotes the general class of DoF and 'suffix' denotes the particular
/// designation of the type (e.g., STRAIN_F, STRAIN_GL, disp_x, disp_y, etc)

class VariableLabeler : public FunctionVisitor {
 public:
  VariableLabeler(std::string const &m_type_name, std::string const &_template);

  std::string type_name() const { return "VariableLabeler"; }

  bool visit(Variable &host, BasisSet const *bset_ptr) const;

 private:
  std::string m_type_name;
  std::vector<std::string> m_sub_strings;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class SubExpressionLabeler : public FunctionVisitor {
 public:
  SubExpressionLabeler(const std::string &_bset_name,
                       const std::string &_template);

  SubExpressionLabeler(SubExpressionLabeler const &) = default;

  std::string type_name() const { return "SubExpressionLabeler"; }

 private:
  bool _generic_visit(Function &host, BasisSet const *bset_ptr) const;

  bool _generic_visit(Function const &host, BasisSet const *bset_ptr) const {
    throw std::runtime_error(
        "Application of a SubExpressionLabeler to a const Function is not "
        "allowed!!\n");
  }

  std::string m_bset_name;
  std::vector<std::string> m_sub_strings;
};

}  // namespace CASM
#endif
