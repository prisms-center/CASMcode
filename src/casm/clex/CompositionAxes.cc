#include "casm/clex/CompositionAxes_impl.hh"
#include "casm/clex/CompositionConverter.hh"

namespace CASM {

void CompositionAxes::erase_enumerated() {
  for (std::string const &el : enumerated) {
    all_axes.erase(el);
    if (curr_key == el) curr_key.clear();
  }
  enumerated.clear();
}

/// \brief Set this->curr using key
void CompositionAxes::select(std::string key) {
  auto it = all_axes.find(key);
  if (it == all_axes.end()) {
    std::stringstream ss;
    ss << "Warning: The composition axes " << key
       << " cannot be found among the posible composition axes.\n\n"
       << "Please use 'casm composition --select' to re-select your "
          "composition axes,\n"
       << "use 'casm composition --calc' to re-calc your standard axes,\n"
       << "or add custom composition axes manually.";

    err_message = ss.str();

    err_code = 2;
  } else {
    curr = it->second;
    curr_key = key;
  }
}

}  // namespace CASM
