#ifndef CASM_clex_CompositionAxes
#define CASM_clex_CompositionAxes

#include <map>
#include <set>
#include <string>

#include "casm/clex/CompositionConverter.hh"

namespace CASM {

/// Data structure to facilitate reading, writing, selecting composition axes in
/// a CASM project
struct CompositionAxes {
  CompositionAxes() {}

  /// \brief Iterate over list of CompositionConverter and insert each one as an
  /// enumerated axes set with a unique numerical name
  template <typename IterType>
  void insert_enumerated(IterType begin, IterType end);

  /// \brief Erase all enumerated axes and clear this->enumerated
  void erase_enumerated();

  /// \brief Set this->curr using key
  void select(std::string key);

  /// \brief True if curr_key is set
  bool has_current_axes() const { return !curr_key.empty(); }

  std::map<std::string, CompositionConverter> all_axes;
  std::set<std::string> enumerated;
  std::string curr_key;
  CompositionConverter curr;

  int err_code = 0;
  std::string err_message;
};

}  // namespace CASM

#endif
