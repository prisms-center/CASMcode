#ifndef CASM_json_optional
#define CASM_json_optional

namespace CASM {

template <typename T>
struct jsonConstructor;
class jsonParser;

template <typename T, typename... Args>
jsonParser &to_json(std::optional<T> const &value, jsonParser &json) {
  if (value.has_value()) {
    to_json(value.value(), json);
  } else {
    json.put_null();
  }
  return json;
}

template <typename T, typename... Args>
void from_json(std::optional<T> &value, jsonParser const &json,
               Args &&... args) {
  if (json.is_null()) {
    value.reset();
  } else {
    value = jsonConstructor<T>::from_json(json, std::forward<Args>(args)...);
  }
}

/// Helper struct for constructing std::optional<T> objects
template <typename T>
struct jsonConstructor<std::optional<T>> {
  /// \brief Default from_json is equivalent to \code
  /// CASM::from_json<ReturnType>(json) \endcode
  template <typename... Args>
  static std::optional<T> from_json(jsonParser const &json, Args &&... args) {
    if (json.is_null()) {
      return std::nullopt;
    } else {
      std::optional<T> value =
          jsonConstructor<T>::from_json(json, std::forward<Args>(args)...);
      return value;
    }
  }
};

}  // namespace CASM

#endif
