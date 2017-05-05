#ifndef CASM_ContainerTraits_HH
#define CASM_ContainerTraits_HH

#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {


  /// \brief Specialize container traits
  /// Expects:
  /// - typedef _Container Container;
  /// - typedef value_type;
  /// - typedef size_type;
  /// - typedef Access; // such as CASM_TMP::BracketAccess or CASM_TMP::ParenthesesAccess
  /// - size_type size(const Eigen::VectorXd& vec) const
  ///
  template<typename _Container>
  struct ContainerTraits {};

  /// \brief Eigen::VectorXd container traits
  template<>
  struct ContainerTraits<Eigen::VectorXd> {
    typedef Eigen::VectorXd Container;
    typedef Eigen::VectorXd::Index size_type;
    typedef Eigen::VectorXd::Scalar value_type;
    typedef CASM_TMP::ParenthesesAccess<Container, value_type, size_type> Access;

    /// \brief Return size of container
    static size_type size(const Eigen::VectorXd &vec) {
      return vec.size();
    }
  };

  /// \brief Eigen::VectorXd container traits
  template<>
  struct ContainerTraits<Eigen::VectorXi> {
    typedef Eigen::VectorXi Container;
    typedef Eigen::VectorXi::Index size_type;
    typedef Eigen::VectorXi::Scalar value_type;
    typedef CASM_TMP::ParenthesesAccess<Container, value_type, size_type> Access;

    /// \brief Return size of container
    static size_type size(const Eigen::VectorXi &vec) {
      return vec.size();
    }
  };

}

#endif
