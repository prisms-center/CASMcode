#ifndef CASM_ClexDescription
#define CASM_ClexDescription

#include <string>
#include <iosfwd>

namespace CASM {

  class jsonParser;
  class DirectoryStructure;

  /// \brief Specifies a particular cluster expansion
  ///
  /// - Cluster expansions are given an identifying name as a shortcut
  /// - Cluster expansions are fully specified via:
  ///   - property: the property being expanded, for instance 'formation_energy'
  ///   - calctype: the type of calculations of the property the cluster
  ///     expansion is fit to
  ///   - ref: indicates a reference used to calculate the property being
  ///     expanded
  ///   - bset: the basis set used
  ///
  struct ClexDescription {
    ClexDescription() {}

    ClexDescription(std::string _name,
                    std::string _property,
                    std::string _calctype,
                    std::string _ref,
                    std::string _bset,
                    std::string _eci) :
      name(_name), property(_property), calctype(_calctype), ref(_ref), bset(_bset), eci(_eci) {}

    void print(std::ostream &sout, bool is_default, int indent = 0) const;

    std::string name;
    std::string property;
    std::string calctype;
    std::string ref;
    std::string bset;
    std::string eci;
  };

  /// \brief Compare using name strings: A.name < B.name
  bool operator<(const ClexDescription &A, const ClexDescription &B);

  jsonParser &to_json(const ClexDescription &desc, jsonParser &json);

  void from_json(ClexDescription &desc, const jsonParser &json);

  bool clex_exists(const DirectoryStructure &dir, const ClexDescription &desc);

}

#endif
