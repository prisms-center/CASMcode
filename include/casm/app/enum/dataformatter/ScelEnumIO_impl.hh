#ifndef CASM_enumerator_io_formatter_ScelEnumIO_impl
#define CASM_enumerator_io_formatter_ScelEnumIO_impl

#include "casm/app/enum/dataformatter/ScelEnumIO.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/database/ConfigDatabaseTools_impl.hh"
#include "casm/misc/TypeInfo.hh"

namespace CASM {

// Note:
// - The following formatters are all templated to expect a ScelEnumDataType
//   that follows the minimum interface defined by `ConfigEnumData`
// - Besides the functions included here, you can also use
//   `make_datum_formatter_adapter` to adapt DatumFormatters from ScelIO for
//   use with ScelEnumDataType.
//   Examples:
//   - `make_datum_formatter_adapter<ScelEnumDataType,
//   Configuration>(ScelIO::scelname())`
//   - `make_datum_formatter_adapter<ScelEnumDataType,
//   Supercell>(ScelIO::TransfMat())`
//   - Note that `make_datum_formatter_adapter` has overrides that can change
//   the formatter name
//     if that is appropriate.
//   - Using `make_datum_formatter_adapter` requires defining an
//   adapter::Adapter

namespace adapter {

template <typename EnumeratorType>
Supercell const &Adapter<Supercell, ScelEnumData<EnumeratorType>>::operator()(
    ScelEnumData<EnumeratorType> const &adaptable) const {
  return adaptable.supercell;
}

}  // namespace adapter

namespace ScelEnumIO {

template <typename ScelEnumDataType>
GenericDatumFormatter<std::string, ScelEnumDataType> name() {
  return GenericDatumFormatter<std::string, ScelEnumDataType>(
      "name", "Supercell name",
      [](ScelEnumDataType const &data) -> std::string {
        return data.supercell.name();
      });
}

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> selected() {
  return GenericDatumFormatter<bool, ScelEnumDataType>(
      "selected", "True if supercell exists in database, else false",
      [](ScelEnumDataType const &data) {
        if (data.is_excluded_by_filter) {
          return false;
        }
        return data.insert_result.first !=
               data.primclex.template db<Supercell>().end();
      });
}

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> is_new() {
  return GenericDatumFormatter<bool, ScelEnumDataType>(
      "is_new",
      "True if supercell did not exist before insertion. False "
      "if not checked or was existing before insertion.",
      [](ScelEnumDataType const &data) {
        if (data.is_excluded_by_filter) {
          return false;
        }
        return data.insert_result.second;
      });
}

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> is_existing() {
  return GenericDatumFormatter<bool, ScelEnumDataType>(
      "is_existing",
      "True if supercell did exist before insertion. False if "
      "not checked or not existing before insertion.",
      [](ScelEnumDataType const &data) {
        if (data.is_excluded_by_filter) {
          return false;
        }
        return !data.insert_result.second;
      });
}

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> is_excluded_by_filter() {
  return GenericDatumFormatter<bool, ScelEnumDataType>(
      "is_excluded_by_filter",
      "True if supercell was excluded by the filter, else false.",
      [](ScelEnumDataType const &data) { return data.is_excluded_by_filter; });
}

template <typename ScelEnumDataType>
Generic2DDatumFormatter<Eigen::MatrixXi, ScelEnumDataType>
transformation_matrix_to_super() {
  return Generic2DDatumFormatter<Eigen::MatrixXi, ScelEnumDataType>(
      "transformation_matrix_to_super",
      "Transformation matrix T, defining the supercell lattice vectors, S, in "
      "terms of the prim lattice vectors, P: `S = P * T`, where S and P are "
      "column vector matrices.",
      [](ScelEnumDataType const &data) -> Eigen::MatrixXi {
        return data.supercell.sym_info()
            .transformation_matrix_to_super()
            .template cast<int>();
      });
}

template <typename ScelEnumDataType>
Generic2DDatumFormatter<Eigen::MatrixXd, ScelEnumDataType>
supercell_lattice_column_matrix() {
  return Generic2DDatumFormatter<Eigen::MatrixXd, ScelEnumDataType>(
      "supercell_lattice_column_matrix",
      "Supercell lattice vectors, as columns of a 3x3 matrix, S = P * T, where "
      "P are the prim lattice vectors, as columns of a 3x3 matrix, and T is a "
      "3x3 integer transformation matrix.",
      [](ScelEnumDataType const &data) -> Eigen::MatrixXd {
        return data.supercell.sym_info().supercell_lattice().lat_column_mat();
      });
}

template <typename ScelEnumDataType>
Generic2DDatumFormatter<Eigen::MatrixXd, ScelEnumDataType>
supercell_lattice_row_vectors() {
  return Generic2DDatumFormatter<Eigen::MatrixXd, ScelEnumDataType>(
      "supercell_lattice_row_vectors",
      "Supercell lattice vectors, as rows of a 3x3 matrix.",
      [](ScelEnumDataType const &data) -> Eigen::MatrixXd {
        return data.supercell.sym_info()
            .supercell_lattice()
            .lat_column_mat()
            .transpose();
      });
}

}  // namespace ScelEnumIO

}  // namespace CASM

#endif
