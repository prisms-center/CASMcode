#ifndef CASM_clexulator_ClexParamPack
#define CASM_clexulator_ClexParamPack
#include <cstddef>
#include <map>
#include <vector>

#include "casm/misc/cloneable_ptr.hh"

namespace CASM {
namespace clexulator {
/** \ingroup Clexulator
 * @{ */

// ClexParamPack is an abstract data structure for storing and retrieving 2D and
// 1D Arrays of real values The data stored in a ClexParamPack is used to
// perform calculations within a Clexulator. Values stores within a
// ClexParamPack represent independent variables, intermediate computations, and
// final evaluated functions. ClexParamPack allows the underlying Scalar value
// type to be obscured, in order to implement special functionality, such as
// automatic differentiation Clexulators may implement their own ClexParamPack,
// or utilitize a pre-defined ClexParamPack specialization
class ClexParamPack;

// ClexParamKey is used externally to Clexulator and ClexParamPack to access
// stored and derived values from ClexParamPack. ClexParamKey holds a pointer to
// ClexParamPack_impl::BaseKey, which implements read and write functionality.
// Keys are allocated and registered by ClexParamPack.
class ClexParamKey;

// Struct for templated read/write access of parameter values
template <typename Scalar>
struct ValAccess {};

namespace ClexParamPack_impl {
/// \brief BaseKey is base class from which all ClexParamKey implementation
/// classes inherit. Hides implementation-specific access attributes through an
/// abstract interface. A ClexParamKey object manages access for a particular
/// type of data valuess. Additional specification of referenced data can be
/// provided via linear or pair indices as necessary
class BaseKey {
 public:
  /// \brief Must be constructed with at least a name
  /// \param _name { unique name that describes managed data }
  /// \param _standalone { true if key refers directly to managed data, false if
  /// refers to value derived from managed data } \param _offset { list of
  /// constant offsets for each supplemental linear identifier index } \param
  /// _stride { list of strides to convert supplemental pair identifier indices
  /// to linear indices}
  BaseKey(std::string const &_name, bool _standalone,
          std::vector<Index> const &_offset = {},
          std::vector<Index> const &_stride = {})
      : m_offset(_offset),
        m_stride(_stride),
        m_identifiers(_offset.size(), 0),
        m_name(_name),
        m_standalone(_standalone) {}

  /// \brief Abstract class must define virtual destructor
  virtual ~BaseKey() {}

  /// \brief Name of this Key, should be descriptive of managed data. Ex.:
  /// "disp_var" of "magspin_site_basis"
  std::string const &name() const { return m_name; }

  /// \brief Specify supplemental identifiers, either via linear indices or pair
  /// indices. Internally, identifiers are always passed as linear indices. The
  /// number and specification of identifiers depends on derived ClexParamKey
  /// implementation. Linear identifiers should be passed as an integral value
  /// (Index, int, etc) pair index identifiers must be passed as
  /// std::pair<Index,Index>
  template <typename... Args>
  void set_identifiers(Args const &...args) {
    _set_identifiers(0, args...);
  }

  /// \brief Clone the ClexParamKey
  std::unique_ptr<BaseKey> clone() const {
    return std::unique_ptr<BaseKey>(_clone());
  }

  /// \brief Returns true if key refers directly to managed data, false if it
  /// refers to values derived from managed data
  bool standalone() const { return m_standalone; }

 protected:
  /// \brief Returns i'th linear index identifier
  Index _l(Index i) const { return m_identifiers[i]; }

  /// \brief Recursive unrolling of variadic template to set all identifiers
  template <typename T, typename... Args>
  void _set_identifiers(Index i, T const &_val, Args const &...args) {
    _set_identifiers(i, _val);
    _set_identifiers(++i, args...);
  }

  /// \brief Overload for linear indexed identifier
  void _set_identifiers(Index i, Index _ind) {
    m_identifiers[i] = m_offset[i] + _ind;
  }

  /// \brief Overload for pair indexed identifier
  void _set_identifiers(Index i, std::pair<Index, Index> _inds) {
    _set_identifiers(i, m_stride[i] * _inds.first + _inds.second);
  }

  /// \brief Clone the ClexParamKey
  virtual BaseKey *_clone() const = 0;

  std::vector<Index> m_offset;
  std::vector<Index> m_stride;
  std::vector<Index> m_identifiers;

 private:
  std::string m_name;
  bool m_standalone;
};
}  // namespace ClexParamPack_impl

/// \brief ParamPackMixIn is interface class to control ClexParamPack portion of
/// Clexulator printing Used primarily by ClexBasisWriter
class ParamPackMixIn {
 public:
  /// \brief static factory function for BasicClexParamPack
  static ParamPackMixIn basic_mix_in() {
    return ParamPackMixIn("BasicClexParamPack",
                          "casm/clexulator/BasicClexParamPack.hh",
                          {{"ParamPack::DEFAULT", "double"}});
  }

  /// \brief static factory function for DiffClexParamPack
  static ParamPackMixIn diff_mix_in() {
    return ParamPackMixIn("DiffClexParamPack",
                          "casm/clexulator/DiffClexParamPack.hh",
                          {{"ParamPack::DEFAULT", "double"},
                           {"ParamPack::DIFF", "ParamPack::DiffScalar"}});
  }

  /// \brief Constructor
  /// Constructed with
  /// \param _name Name of derived ClexParamPack specialization. Expected to be
  ///     a type that exists in the namespace `CASM::clexulator`. Ex:
  ///     "BasicClexParamPack".
  /// \param _filename File to include for ClexParamPack specialization. Ex:
  ///      "casm/clexulator/BasicClexParamPack.hh"
  /// \param _specializations See `basic_mix_in` and `diff_mix_in` examples.
  ParamPackMixIn(std::string const &_name, std::string const &_filename,
                 std::map<std::string, std::string> const &_specializations)
      : m_name(_name),
        m_filename(_filename),
        m_scalar_specializations(_specializations) {}

  /// \brief Abstract class must define virtual destructor
  virtual ~ParamPackMixIn() {}

  /// \brief typename of the corresponding ClexParamPack
  std::string const &name() const { return m_name; }

  /// \brief Dictionary of pairs ("EvalMode", "ScalarType")
  ///  These correspond to the underlying scalar type to be used for each
  ///  Evaluation mode
  std::map<std::string, std::string> const &scalar_specializations() const {
    return m_scalar_specializations;
  }

  /// \brief returns string with include directives for Clexulator.
  virtual std::string cpp_includes_string() const {
    return "#include \"" + m_filename + "\"\n";
  }

  /// \brief returns string with c++ definitions to be written before
  /// Clexulator.
  virtual std::string cpp_definitions_string(std::string const &_indent) const {
    return "";
  }

  /// \brief Clone the ParamPackMixIn
  std::unique_ptr<ParamPackMixIn> clone() const {
    return std::unique_ptr<ParamPackMixIn>(_clone());
  }

 private:
  /// \brief Clone the ParamPackMixIn (implementation)
  virtual ParamPackMixIn *_clone() const { return new ParamPackMixIn(*this); }

  std::string m_name;
  std::string m_filename;
  std::map<std::string, std::string> m_scalar_specializations;
};

/// \brief Key for indexing clexulator parameters
/// Contains pointer to implementation of the Key
class ClexParamKey {
 public:
  /// \brief Allow default construction
  ClexParamKey() {}

  /// \brief Construct ClexParamKey by cloning specialization of BaseKey
  ClexParamKey(ClexParamPack_impl::BaseKey const &_key)
      : m_key_ptr(_key.clone()) {}

  /// \brief Returns name of key
  std::string const &name() const { return m_key_ptr->name(); }

  /// \brief Returns pointer to Key implementation
  ClexParamPack_impl::BaseKey const *ptr() const {
    return m_key_ptr.unique().get();
  }

  /// \brief use parentheses operator to set identifiers.
  /// Permits syntax like
  /// \code
  ///  ClexParamPack * pp_ptr = &(my_clexulator.param_pack());
  ///  ClexParamKey pp_ptr->key("some_inherently_5dimensional_parameter");
  ///  // accesses slice of 5-dimensional data
  ///  Eigen::MatrixXd my_params = pp_ptr->read(key(i,
  ///  std::pair<Index,Index>(j,k)));
  /// \endcode
  template <typename... Args>
  ClexParamKey &operator()(Args const &...args) {
    m_key_ptr->set_identifiers(args...);
    return *this;
  }

 private:
  /// \brief  ptr to BaseKey class that hides implementation-specific access
  /// attributes
  notstd::cloneable_ptr<ClexParamPack_impl::BaseKey> m_key_ptr;
};

/// \brief Abstract base class for reading/writing clexulator parameters
/// Parameters are assumed to be stored as, or at least sliceable into, 2D or 1D
/// arrays in order to standardize read/write interface via Eigen types Arrays
/// are layed out along two dimensions to form (m x N) or (m x 1) matrices where
/// 'N' is the number of sites in the neighborhood of a unit cell, and 'm' is
/// the intrinsic dimension of a particular parameter. Macroscopic or
/// homogeneous properties of the crystal, such as evaluated crystal basis
/// functions or homogeneous strain, which are not associated with individual
/// sites, are represented as a (m x 1) array. The length of the 'm' dimension
/// of a parameter array is referred do as 'dim' and the length of the N'
/// dimension is referred to as 'size'
class ClexParamPack {
 public:
  typedef unsigned int size_type;

  /// \brief Abstract class must define virtual destructor
  virtual ~ClexParamPack() {}

  /// \brief Obtain registry of all keys for data blocks managed by this
  /// ClexParamPack
  std::map<std::string, ClexParamKey> const &keys() const { return m_keys; }

  /// \brief Obtain key for managed data block by name
  ClexParamKey const &key(std::string const &_name) const {
    auto it = keys().find(_name);
    if (it == keys().end()) {
      std::stringstream ss;
      for (it = keys().begin(); it != keys().end(); ++it)
        ss << "  " << it->first << "\n";
      throw std::runtime_error(
          "In ClexParamPack::key(), ClexParamPack does not contain parameters "
          "corresponding to name " +
          _name + ". Options are:\n" + ss.str());
    }
    return it->second;
  }

  /// \brief 'N' dimension of parameter array (either '1', or size of unit cell
  /// neighborhood) returns number of columns in parameter array
  virtual size_type size(ClexParamKey const &_key) const = 0;

  /// \brief 'm' dimension of parameter array
  /// returns number of rows in parameter array
  virtual size_type dim(ClexParamKey const &_key) const = 0;

  /// \brief Check evaluation mode for ClexParamPack parameter
  /// Choices are at least
  /// - "default" (Clexulator manages evaluation and writing of parameters)
  /// - "read" (User fixes values via ClexParamPack::write and Clexulator reads
  /// provided values)
  /// - "custom" (Clexulator utilizes user-provided function handles to evaluate
  /// parameters)
  virtual std::string eval_mode(ClexParamKey const &_key) const = 0;

  /// \brief Returns const reference to parameter array for parameter specified
  /// by \param _key NOTE: Take care in constructing pointers or references from
  /// the return value. If _key accepts additional index identifiers, (because
  /// the parameter has more than 2 dimensions) ClexParamPack::read() may be
  /// implemented using  an internally cached matrix holding the current data
  /// slice. As such, subsequent requests to ClexParamPack::read() for a
  /// different slice may alter the internally cached values while returning a
  /// reference for the same block of memory
  virtual Eigen::MatrixXd const &read(ClexParamKey const &_key) const = 0;

  /// \brief Returns const reference to element of 1D parameter array for
  /// parameter specified by \param _key
  virtual double const &read(ClexParamKey const &_key, size_type _i) const = 0;

  /// \brief Returns const reference to element of 2D parameter array for
  /// parameter specified by \param _key \param _i {row index} and \param _j
  /// {column index}
  virtual double const &read(ClexParamKey const &_key, size_type _i,
                             size_type _j) const = 0;

  /// \brief Sets evaluation mode for ClexParamPack parameter
  /// Choices are at least
  /// - "default" (Clexulator manages evaluation and writing of parameters)
  /// - "read" (User fixes values via ClexParamPack::write and Clexulator reads
  /// provided values)
  /// - "custom" (Clexulator utilizes user-provided function handles to evaluate
  /// parameters)
  virtual void set_eval_mode(ClexParamKey const &_key,
                             std::string const &_mode) = 0;

  /// \brief Write parameter array for parameter specified by \param _key
  virtual void write(ClexParamKey const &_key,
                     Eigen::Ref<const Eigen::MatrixXd> const &_val) = 0;

  /// \brief Write element to 1D parameter array for parameter specified by
  /// \param _key
  virtual void write(ClexParamKey const &_key, size_type _i, double val) = 0;

  /// \brief Write element to 2D parameter array for parameter specified by
  /// \param _key \param _i {row index} and \param_j {column index}
  virtual void write(ClexParamKey const &_key, size_type _i, size_type _j,
                     double val) = 0;

  /// \brief May be specialized to perform preprocessing before function
  /// evaluation
  void pre_eval() {}

  /// \brief May be specialized to perform postprocessing after function
  /// evaluation
  void post_eval() {}

 protected:
  std::map<std::string, ClexParamKey> m_keys;

 private:
  // possible implementation:
  // std::vector<Eigen::MatrixXd> m_data;
};

/** @} */
}  // namespace clexulator
}  // namespace CASM

#endif
