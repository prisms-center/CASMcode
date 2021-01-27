#ifndef CASM_ContainerTraits_HH
#define CASM_ContainerTraits_HH

#include <vector>

#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {

/// \brief Specialize container traits
/// Expects:
/// - typedef _Container Container;
/// - typedef value_type;
/// - typedef size_type;
/// - typedef Access; // such as CASM_TMP::BracketAccess or
/// CASM_TMP::ParenthesesAccess
/// - size_type size(const Eigen::VectorXd& vec) const
///
template <typename _Container>
struct ContainerTraits {};

template <typename _value_type, typename Allocator>
struct ContainerTraits<std::vector<_value_type, Allocator> > {
  using Container = std::vector<_value_type, Allocator>;
  using value_type = typename Container::value_type;
  using size_type = typename Container::size_type;
  using Access = CASM_TMP::BracketAccess<Container, value_type, size_type>;

  static size_type size(const Container &vec) { return vec.size(); }

  static size_type rows(const Container &vec) { return vec.size(); }

  static size_type cols(const Container &vec) { return 1; }
};

/// \brief Eigen::MatrixXd container traits

template <typename _value_type, int a, int b, int c, int d, int e>
struct ContainerTraits<Eigen::Matrix<_value_type, a, b, c, d, e> > {
  typedef Eigen::Matrix<_value_type, a, b, c, d, e> Container;
  typedef typename Container::Index size_type;
  typedef typename Container::Scalar value_type;
  typedef typename Container::Scalar value_type2D;
  typedef CASM_TMP::ParenthesesAccess<Container, value_type, size_type> Access;

  /// \brief Return size of container
  static size_type size(const Container &mat) { return mat.size(); }

  /// \brief Return size of container
  static size_type rows(const Container &mat) { return mat.rows(); }

  /// \brief Return size of container
  static size_type cols(const Container &mat) { return mat.cols(); }
};
/*
/// \brief Eigen::VectorXd container traits
template<>
struct ContainerTraits<Eigen::VectorXd> {
  typedef Eigen::VectorXd Container;
  typedef Eigen::VectorXd::Index size_type;
  typedef Eigen::VectorXd::Scalar value_type;
  typedef Eigen::VectorXd::Scalar value_type2D;
  typedef CASM_TMP::ParenthesesAccess<Container, value_type, size_type> Access;

  /// \brief Return size of container
  static size_type size(const Eigen::VectorXd &vec) {
    return vec.size();
  }

  /// \brief Return size of container
  static size_type size(const Eigen::VectorXd &vec) {
    return vec.size();
  }

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
*/
}  // namespace CASM

#endif
