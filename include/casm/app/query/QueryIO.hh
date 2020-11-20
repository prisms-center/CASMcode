#ifndef CASM_app_query_QueryIO
#define CASM_app_query_QueryIO

namespace CASM {

  class PermuteIterator;
  template<typename ValueType>
  struct QueryData;

  namespace adapter {

    template <typename ToType, typename FromType> struct Adapter;

    template<typename DataObject>
    struct Adapter<DataObject, QueryData<DataObject>> {
      DataObject const &operator()(QueryData<DataObject> const &adaptable) const;
    };

  }

  /// Collect information during `enumerate_configurations` function for optional output
  template<typename ValueType>
  struct QueryData {

    QueryData(PrimClex const &_primclex,
              ValueType const &_value,
              bool _include_equivalents = false,
              Index _equivalent_index = 0,
              PermuteIterator const *_permute_it = nullptr);

    PrimClex const &primclex;
    ValueType const &value;
    bool include_equivalents;
    Index equivalent_index;
    PermuteIterator const *permute_it;
  };

  namespace QueryIO {

    template <typename QueryDataType>
    GenericDatumFormatter<Index, QueryDataType> equivalent_index();

    template <typename QueryDataType>
    GenericDatumFormatter<Index, QueryDataType> permute_scel_factor_group_op();

    template <typename QueryDataType>
    GenericDatumFormatter<Index, QueryDataType> permute_factor_group_op();

    template <typename QueryDataType>
    GenericDatumFormatter<std::string, QueryDataType> permute_factor_group_op_desc();

    template <typename QueryDataType>
    Generic1DDatumFormatter<Eigen::Vector3l, QueryDataType> permute_translation();

  }

}

#endif
