#ifndef CASM_app_query_QueryIO_impl
#define CASM_app_query_QueryIO_impl

#include "casm/app/query/QueryIO.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"

namespace CASM {

  namespace adapter {

    template<typename DataObject>
    DataObject const &Adapter<DataObject, QueryData<DataObject>>::operator()(
    QueryData<DataObject> const &adaptable) const {
      return adaptable.value;
    }

  }

  template<typename ValueType>
  QueryData<ValueType>::QueryData(PrimClex const &_primclex,
                                  ValueType const &_value,
                                  bool _include_equivalents,
                                  Index _equivalent_index,
                                  PermuteIterator const *_permute_it):
    primclex(_primclex),
    value(_value),
    include_equivalents(_include_equivalents),
    equivalent_index(_equivalent_index),
    permute_it(_permute_it) {}

  namespace QueryIO {

    template <typename QueryDataType>
    GenericDatumFormatter<Index, QueryDataType> equivalent_index() {
      return GenericDatumFormatter<Index, QueryDataType>(
               "equivalent_index",
               "Index that counts over distinct symmetrically equivalent objects",
      [](QueryDataType const & data) -> Index {
        return data.equivalent_index;
      });
    }

    template <typename QueryDataType>
    GenericDatumFormatter<Index, QueryDataType> permute_scel_factor_group_op() {
      return GenericDatumFormatter<Index, QueryDataType>(
               "permute_scel_factor_group_op",
               "Factor group operation (applied before translation operation) of the permutation "
               "that generated the equivalent object, as its index into the prim factor group",
      [](QueryDataType const & data) -> Index {
        return data.permute_it->factor_group_index();
      });
    }

    template <typename QueryDataType>
    GenericDatumFormatter<Index, QueryDataType> permute_factor_group_op() {
      return GenericDatumFormatter<Index, QueryDataType>(
               "permute_factor_group_op",
               "Factor group operation (applied before translation operation) of the permutation "
               "that generated the equivalent object, as its index into the prim factor group",
      [](QueryDataType const & data) -> Index {
        return data.permute_it->prim_factor_group_index();
      });
    }

    template <typename QueryDataType>
    GenericDatumFormatter<std::string, QueryDataType> permute_factor_group_op_desc() {
      return GenericDatumFormatter<std::string, QueryDataType>(
               "permute_factor_group_op_desc",
               "Description of prim factor group operation that generated equivalent object",
      [](QueryDataType const & data) -> std::string {

        auto const &factor_group_op = data.permute_it->factor_group()[data.permute_it->factor_group_index()];
        auto const &lattice = data.primclex.prim().lattice();
        std::stringstream ss;
        ss << "'" << brief_description(factor_group_op, lattice) << "'";
        return ss.str();
      });
    }

    template <typename QueryDataType>
    Generic1DDatumFormatter<Eigen::Vector3l, QueryDataType> permute_translation() {
      return Generic1DDatumFormatter<Eigen::Vector3l, QueryDataType>(
               "permute_translation",
               "Translation (applied after factor group operation) of the permutation that generated "
               "the equivalent object, as a linear combination of the lattice vectors",
      [](QueryDataType const & data) -> Eigen::Vector3l {
        return data.permute_it->sym_info().unitcell_index_converter()(data.permute_it->translation_index());
      });
    }

  }

}

#endif
