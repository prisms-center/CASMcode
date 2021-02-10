#ifndef CASM_DataFormatterDecl
#define CASM_DataFormatterDecl

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

class DataStream;

// --- DataFormatterTools ---

template <typename DataObject>
class BaseDatumFormatter;

template <typename _DataObject,
          typename _DatumFormatterType = BaseDatumFormatter<_DataObject>>
class DataFormatterDictionary;

template <typename _DataObject>
class DataFormatter;

class FormattedPrintable;

template <typename DataObject, typename DatumFormatterType>
struct DictionaryConverter;

template <typename T>
bool always_true(const T &) {
  return true;
};

// --- DataFormatterTools ---

template <typename ValueType, typename ArgType, typename DataObject>
class DataFormatterOperator;

template <typename DataObject>
class DatumFormatterAlias;

template <typename ValueType, typename DataObject>
class ConstantValueFormatter;

template <typename ValueType, typename DataObject>
class BaseValueFormatter;

template <typename ValueType, typename DataObject>
class GenericDatumFormatter;

template <typename Container, typename DataObject>
class Base1DDatumFormatter;

template <typename Container, typename DataObject>
class Generic1DDatumFormatter;

template <typename Container, typename DataObject>
class Base2DDatumFormatter;

template <typename Container, typename DataObject>
class Generic2DDatumFormatter;

/// \brief Dictionary of all Attribute types of DataObject
///
/// Default includes String, Boolean, Integer, Scalar, and VectorXd attributes
template <typename DataObject>
DataFormatterDictionary<DataObject> make_attribute_dictionary();

/// \brief Template alias for BaseValueFormatter returning std::string
///
/// \ingroup DataFormatterTypes
///
template <typename DataObject>
using StringAttribute = BaseValueFormatter<std::string, DataObject>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
using StringAttributeDictionary =
    DataFormatterDictionary<DataObject, StringAttribute<DataObject>>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
StringAttributeDictionary<DataObject> make_string_dictionary();

/// \brief Template alias for BaseValueFormatter returning bool
///
/// \ingroup DataFormatterTypes
///
template <typename DataObject>
using BooleanAttribute = BaseValueFormatter<bool, DataObject>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
using BooleanAttributeDictionary =
    DataFormatterDictionary<DataObject, BooleanAttribute<DataObject>>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
BooleanAttributeDictionary<DataObject> make_boolean_dictionary();

/// \brief Template alias for BaseValueFormatter returning Index
///
/// \ingroup DataFormatterTypes
///
template <typename DataObject>
using IntegerAttribute = BaseValueFormatter<Index, DataObject>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
using IntegerAttributeDictionary =
    DataFormatterDictionary<DataObject, IntegerAttribute<DataObject>>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
IntegerAttributeDictionary<DataObject> make_integer_dictionary();

/// \brief Template alias for BaseValueFormatter returning double
///
/// \ingroup DataFormatterTypes
///
template <typename DataObject>
using ScalarAttribute = BaseValueFormatter<double, DataObject>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
using ScalarAttributeDictionary =
    DataFormatterDictionary<DataObject, ScalarAttribute<DataObject>>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
ScalarAttributeDictionary<DataObject> make_scalar_dictionary();

/// \brief Template alias for BaseValueFormatter returning Eigen::VectorXi
///
/// \ingroup DataFormatterTypes
///
template <typename DataObject>
using VectorXiAttribute = Base1DDatumFormatter<Eigen::VectorXi, DataObject>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
using VectorXiAttributeDictionary =
    DataFormatterDictionary<DataObject, VectorXiAttribute<DataObject>>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
VectorXiAttributeDictionary<DataObject> make_vectorxi_dictionary();

/// \brief Template alias for BaseValueFormatter returning Eigen::VectorXd
///
/// \ingroup DataFormatterTypes
///
template <typename DataObject>
using VectorXdAttribute = Base1DDatumFormatter<Eigen::VectorXd, DataObject>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
using VectorXdAttributeDictionary =
    DataFormatterDictionary<DataObject, VectorXdAttribute<DataObject>>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
VectorXdAttributeDictionary<DataObject> make_vectorxd_dictionary();

/// \brief Template alias for BaseValueFormatter returning Eigen::MatrixXd
///
/// \ingroup DataFormatterTypes
///
template <typename DataObject>
using MatrixXdAttribute = Base2DDatumFormatter<Eigen::MatrixXd, DataObject>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
using MatrixXdAttributeDictionary =
    DataFormatterDictionary<DataObject, MatrixXdAttribute<DataObject>>;

/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
MatrixXdAttributeDictionary<DataObject> make_matrixxd_dictionary();

class jsonParser;
/// \brief Template to be specialized for constructing dictionaries for
/// particular DataObject
///
/// \ingroup DataFormatter
///
template <typename DataObject>
DataFormatterDictionary<DataObject, BaseValueFormatter<jsonParser, DataObject>>
make_json_dictionary();

}  // namespace CASM

#endif
