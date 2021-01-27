#ifndef CASM_app_query_QueryIO
#define CASM_app_query_QueryIO

namespace CASM {

class PermuteIterator;
template <typename ValueType>
struct QueryData;

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

template <typename DataObject>
struct Adapter<DataObject, QueryData<DataObject>> {
  DataObject const &operator()(QueryData<DataObject> const &adaptable) const;
};

}  // namespace adapter

/// Data structure to collect value (Supercell, Configuration, etc.) and data
/// for query formatting
template <typename ValueType>
struct QueryData {
  QueryData(PrimClex const &_primclex, ValueType const &_value,
            bool _include_equivalents = false, Index _equivalent_index = 0,
            PermuteIterator const *_permute_it = nullptr);

  /// Provides access to project data
  PrimClex const &primclex;

  /// The object under consideration (Supercell, Configuration, etc.)
  ValueType const &value;

  /// True if non-canonical equivalents are also being queried
  bool include_equivalents;

  /// Index incremented for each non-canonical equivalent
  Index equivalent_index;

  /// Permutation that transforms the canonical form to the current
  /// non-canonical equivalent value
  PermuteIterator const *permute_it;
};

namespace QueryIO {

/// Index incremented for each non-canonical equivalent
/// (QueryDataType::equivalent_index)
template <typename QueryDataType>
GenericDatumFormatter<Index, QueryDataType> equivalent_index();

/// Index of permute factor group operation in the supercell factor group
/// (QueryDataType::permute_it->factor_group_index())
template <typename QueryDataType>
GenericDatumFormatter<Index, QueryDataType> permute_scel_factor_group_op();

/// Index of permute factor group operation in the prim factor group
/// (QueryDataType::permute_it->prim_factor_group_index())
template <typename QueryDataType>
GenericDatumFormatter<Index, QueryDataType> permute_factor_group_op();

/// Description of permute factor group operation in the prim factor group
template <typename QueryDataType>
GenericDatumFormatter<std::string, QueryDataType>
permute_factor_group_op_desc();

/// Permute translation operation as (i, j, k) multiples of the prim lattice
/// vectors
template <typename QueryDataType>
Generic1DDatumFormatter<Eigen::Vector3l, QueryDataType> permute_translation();

}  // namespace QueryIO

}  // namespace CASM

#endif
