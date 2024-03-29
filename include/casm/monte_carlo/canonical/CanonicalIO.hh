#ifndef CASM_CanonicalIO_HH
#define CASM_CanonicalIO_HH

namespace CASM {
class PrimClex;
class jsonParser;
template <typename T>
class DataFormatter;
}  // namespace CASM

namespace CASM {
namespace Monte {

class MonteCarlo;
typedef const MonteCarlo *ConstMonteCarloPtr;
class Canonical;
class CanonicalConditions;

/// \brief Make a LTE results formatter
DataFormatter<ConstMonteCarloPtr> make_results_formatter(const Canonical &mc);

/// \brief Store CanonicalConditions in JSON format
jsonParser &to_json(const CanonicalConditions &conditions, jsonParser &json);

/// \brief Read CanonicalConditions from JSON format
void from_json(CanonicalConditions &conditions, const PrimClex &primclex,
               const jsonParser &json, Canonical const &mc,
               bool incremental = false);

}  // namespace Monte
}  // namespace CASM

#endif
