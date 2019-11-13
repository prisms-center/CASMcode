#ifndef CASM_AnisoValTraits
#define CASM_AnisoValTraits

#include <set>
#include <string>
#include <memory>
#include <vector>
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/crystallography/SymRepBuilder.hh"

//Defines notstd::make_unique<>()
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  ///\brief Class for specifying essential traits of anisotropic values (e.g., degrees of freedom, calculated properties, and discretized species attributes)
  class AnisoValTraits {
  public:
    static const unsigned char LOCAL = 0;
    static const unsigned char GLOBAL = (1u << 0);
    static const unsigned char UNIT_LENGTH = (1u << 1);
    static const unsigned char DESCRIBES_ORIENTATION = (1u << 2);

    ///\brief Named constructor for uninitialized AnisoValTraits
    static AnisoValTraits null();

    ///\brief Named constructor for total energy AnisoValTraits
    static AnisoValTraits energy();

    ///\brief Named constructor for site displacement AnisoValTraits
    static AnisoValTraits disp();

    ///\brief Named constructor for site force AnisoValTraits
    static AnisoValTraits force();

    ///\brief Named constructor for global strain AnisoValTraits
    /// @param _metric specifies which strain metric. Choices are:
    ///  - "GL" : Green-Lagrange
    ///  - "AE" : Almansi-Euler
    ///  - "H"  : Hencky
    static AnisoValTraits strain(std::string const &_metric);

    ///\brief Named constructor for site magnetic spin AnisoValTraits
    static AnisoValTraits magspin();

    ///\brief Named constructor for magnetic moment AnisoValTraits
    /// Same as AnisoValTraits::magspin(), but requires unit length
    static AnisoValTraits magmom();

    /// \brief Given a string, returns string with all characters before the final @delim character deleted
    /// For example, if a trajector has properties with  key1="step1_force", key2="step2_force", etc, then
    /// AnisoValTraits::name_suffix(key1) will return "force" if delim='_'
    static std::string name_suffix(std::string const &_name, char delim = '_') {
      std::string result;
      for(char c : _name) {
        if(c == delim)
          result.clear();
        else
          result.push_back(c);
      }
      return result;
    }

    /// \brief Explicit constructor for AnisoValTraits must specify *all* attributes at construction
    /// There is no mutator for AnisoValTraits. Explicitly constructing two AnisoValTraits objects
    /// with the same name but non-identical data members will result in a thrown exception.
    AnisoValTraits(std::string const &_name,
                   std::vector<std::string> const &_std_var_names,
                   unsigned char _options,
                   SymRepBuilderInterface const &_symrep_builder = NullSymRepBuilder(),
                   std::set<std::string> const &_incompatible = {},
                   std::set<std::string> const &_must_apply_before = {},
                   std::set<std::string> const &_must_apply_after = {},
                   bool _default = false);

    /// \brief Returns previously explicitly initialized AnisoValTraits with name AnisoValTraits::name_suffix(_name)
    /// If no AnisoValTraits with matching name has been initialized, throws an exception
    AnisoValTraits(std::string const &_name);

    static std::string class_desc() {
      return "AnisoValTraits";
    }

    /// \brief Generate a symmetry representation for the supporting vector space
    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                    Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                    bool time_reversal) const {
      if(m_symrep_builder != nullptr) {
        return m_symrep_builder->symop_to_matrix(_matrix, _tau, time_reversal, dim());
      }
      //else
      return Eigen::MatrixXd::Identity(dim(), dim());
    }

    /// \brief const access of name
    std::string const &name() const {
      return m_name;
    }

    /// \brief return true if *this has 'default' designation, meaning it can be overridden
    bool is_default() const {
      return m_default;
    }

    /// \brief return 'options' bitflag
    unsigned char options() const {
      return m_opt;
    }

    // Please use !AnisoValTraits::global() instead of implementing AnisoValTraits::local()
    //bool local()const {
    //  return !global();
    //}

    /// \brief returns true if DoF is global
    bool global()const {
      return m_opt & GLOBAL;
    }

    /// \brief returns true if DoF must always have unit length
    bool unit_length()const {
      return m_opt & UNIT_LENGTH;
    }

    /// \brief returns true if time-reversal changes the DoF value
    bool time_reversal_active() const {
      return m_symrep_builder->time_reversal_active();
    }

    /// \brief returns true if value tracks the orientation of an occupying molecule (not typical)
    bool describes_occupant_orientation() const {
      return m_opt & DESCRIBES_ORIENTATION;
    }

    /// \brief conventional dimensionality of this DoF, returns -1 if always variable
    Index dim() const {
      return standard_var_names().size();
    }

    /// \brief equality comparison of name
    bool operator==(std::string const &other_name) const {
      return name() == other_name;
    }

    /// \brief lexicographic comparison of name
    bool operator<(std::string const &other_name) const {
      return name() < other_name;
    }

    /// \brief allow implicit conversion to std::string (name)
    operator std::string const &() const {
      return name();
    }

    /// \brief return standard coordinate axes for continuous variable space
    std::vector<std::string> const &standard_var_names() const {
      return m_standard_var_names;
    }

    std::set<std::string> const &incompatible() const {
      return m_incompatible;
    }

    /// \brief Return list of DoFs that *must* be applied before this DoF is applied
    std::set<std::string> const &must_apply_before() const {
      return m_apply_before;
    }

    /// \brief Return list of DoFs that *must* be applied after this DoF is applied
    std::set<std::string> const &must_apply_after() const {
      return m_apply_after;
    }

    std::unique_ptr<AnisoValTraits> clone() const {
      return notstd::make_unique<AnisoValTraits>(*this);
    }

    std::string symrep_builder_name() const {
      if(m_symrep_builder)
        return m_symrep_builder->name();
      else
        return "NULL";
    }

  protected:
    std::string m_name;
    bool m_default;
    std::vector<std::string> m_standard_var_names;
    unsigned char m_opt;
    notstd::cloneable_ptr<const SymRepBuilderInterface> m_symrep_builder;
    std::set<std::string> m_incompatible;
    std::set<std::string> m_apply_before;
    std::set<std::string> m_apply_after;
  };


  /// \brief comparison of name, domain (discrete/continuous) and mode (local/global)
  inline
  bool identical(AnisoValTraits const &A, AnisoValTraits const &B) {
    return (A.name() == B.name()
            && A.is_default() == B.is_default()
            && A.symrep_builder_name() == B.symrep_builder_name()
            && A.standard_var_names() == B.standard_var_names()
            && A.options() == B.options()
            && A.must_apply_before() == B.must_apply_before()
            && A.must_apply_after() == B.must_apply_after()
            && A.incompatible() == B.incompatible());
  }

  namespace AnisoVal_impl {
    /// \brief A class to manage dynamic evaluation of BasisFunctions
    //struct TraitsConverter {
    /*
    inline
    notstd::cloneable_ptr<DoFType::BasicTraits> traits2cloneable_ptr(const DoFType::BasicTraits &value) {
      return notstd::cloneable_ptr<DoFType::BasicTraits>(value.clone().release());
    }
    */
    //};
  }

}

#endif
