#ifndef SUPERLATTICEENUMERATORIO_HH
#define SUPERLATTICEENUMERATORIO_HH

namespace CASM {
  namespace xtal {
    class ScelEnumProps;
  }
  template<typename T> class InputParser;
  class jsonParser;

  jsonParser &to_json(const xtal::ScelEnumProps &props, jsonParser &json);

  /// Make a ScelEnumProps object from JSON input
  template <>
  struct jsonConstructor<xtal::ScelEnumProps> {
    static xtal::ScelEnumProps from_json(const jsonParser &json);
  };

  /// Make a ScelEnumProps object from JSON input
  void from_json(xtal::ScelEnumProps &props, const jsonParser &json);

  /// Make a ScelEnumProps object from JSON input
  void parse(InputParser<xtal::ScelEnumProps> &parser);
}

#endif